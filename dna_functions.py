from functools import cmp_to_key
from urllib.error import HTTPError
from Bio.Restriction.Restriction import RestrictionBatch, RestrictionType
from Bio.Seq import reverse_complement
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from pydantic_models import PCRSource, PrimerModel, RepositoryIdSource, RestrictionEnzymeDigestionSource, \
    GenbankSequence, SequenceEntity, StickyLigationSource
from pydna.parsers import parse as pydna_parse
from itertools import permutations, product
from pydna.primer import Primer
from pydna.amplify import Anneal
from pydna.amplicon import Amplicon
import requests
from bs4 import BeautifulSoup
import regex
from Bio.SeqFeature import SimpleLocation, Location, CompoundLocation


def sum_is_sticky(seq1: Dseq, seq2: Dseq, partial: bool=False) -> bool:
    """Return true if the 3' end of seq1 and 5' end of seq2 ends are sticky and compatible for ligation."""
    type_seq1, sticky_seq1 = seq1.three_prime_end()
    type_seq2, sticky_seq2 = seq2.five_prime_end()

    if not partial:
        return 'blunt' != type_seq2 and type_seq2 == type_seq1 and str(sticky_seq2) == str(reverse_complement(sticky_seq1))
    else:
        if type_seq1 != type_seq2 or type_seq2 == "blunt":
            return False
        elif type_seq2 == "5'":
            sticky_seq1 = str(reverse_complement(sticky_seq1))
        elif type_seq2 == "3'":
            sticky_seq2 = str(reverse_complement(sticky_seq2))
    
        ovhg_len = min(len(sticky_seq1), len(sticky_seq2))
        for i in range(1, ovhg_len+1):
            if sticky_seq1[-i:] == sticky_seq2[:i]:
                return True
        else:
            return False


def format_sequence_genbank(seq: Dseqrecord) -> SequenceEntity:

    if seq.name.lower() == 'exported':
        correct_name(seq)
    # In principle here we do not set the id from the id of the Dseqrecord
    overhang_crick_3prime, overhang_watson_3prime = both_overhangs_from_dseq(
        seq.seq)
    gb_seq = GenbankSequence(file_content=seq.format('genbank'),
                             overhang_crick_3prime=overhang_crick_3prime,
                             overhang_watson_3prime=overhang_watson_3prime)
    return SequenceEntity(sequence=gb_seq)


def read_primer_from_json(primer: PrimerModel) -> Primer:
    return Primer(
        primer.sequence,
        id=str(primer.id),
        name=primer.name,
    )


def read_dsrecord_from_json(seq: SequenceEntity) -> Dseqrecord:
    initial_dseqrecord: Dseqrecord = pydna_parse(seq.sequence.file_content)[0]
    if seq.sequence.overhang_watson_3prime == 0 and seq.sequence.overhang_crick_3prime == 0:
        out_dseq_record = initial_dseqrecord
    else:
        out_dseq_record = Dseqrecord(dseq_from_both_overhangs(
            str(initial_dseqrecord.seq),
            seq.sequence.overhang_crick_3prime,
            seq.sequence.overhang_watson_3prime))
        out_dseq_record.features = initial_dseqrecord.features
    # We set the id to the integer converted to integer (this is only
    # useful for assemblies)
    out_dseq_record.id = str(seq.id)
    return out_dseq_record


def dseq_from_both_overhangs(full_sequence: str, overhang_crick_3prime: int, overhang_watson_3prime: int) -> Dseq:
    full_sequence_rev = str(Dseq(full_sequence).reverse_complement())
    watson = full_sequence
    crick = full_sequence_rev

    # If necessary, we trim the left side
    if overhang_crick_3prime < 0:
        crick = crick[:overhang_crick_3prime]
    elif overhang_crick_3prime > 0:
        watson = watson[overhang_crick_3prime:]

    # If necessary, we trim the right side
    if overhang_watson_3prime < 0:
        watson = watson[:overhang_watson_3prime]
    elif overhang_watson_3prime > 0:
        crick = crick[overhang_watson_3prime:]

    return Dseq(watson, crick=crick, ovhg=overhang_crick_3prime)


def both_overhangs_from_dseq(dseq: Dseq):
    return dseq.ovhg, len(dseq.watson) - len(dseq.crick) + dseq.ovhg


def sort_by(list2sort, reference_list):
    argsort = sorted(range(len(reference_list)), key=reference_list.__getitem__)
    return [list2sort[i] for i in argsort]


def get_cutsite_order(dseqrecord: Dseqrecord, enzymes: RestrictionBatch) -> tuple[list[RestrictionType], list[int]]:
    """Return the cutsites in order."""
    cuts = enzymes.search(dseqrecord.seq, linear=dseqrecord.linear)
    positions = list()
    cutsites = list()
    for enzyme in cuts:
        for position in cuts[enzyme]:
            # batch.search returns 1-based indexes
            positions.append(position - 1)
            cutsites.append(str(enzyme))

    if len(positions) == 0:
        return cutsites, positions

    # Sort both lists by position
    cutsites = sort_by(cutsites, positions)
    positions = sorted(positions)

    # For digestion of linear sequences, the first one and last one are the molecule ends
    if dseqrecord.linear:
        cutsites.insert(0, '')
        cutsites.append('')
        positions.insert(0, 0)
        positions.append(len(dseqrecord))

    # For digestion of circular sequences, the first one and last one are the same
    if not dseqrecord.linear:
        cutsites.append(cutsites[0])
        positions.append(positions[0])

    return cutsites, positions


def get_invalid_enzyme_names(enzyme_names_list):
    rest_batch = RestrictionBatch()
    invalid_names = list()
    for name in enzyme_names_list:
        # Empty enzyme names are the natural edges of the molecule
        if name != '':
            try:
                rest_batch.format(name)
            except ValueError:
                invalid_names.append(name)
    return invalid_names


def get_restriction_enzyme_products_list(seq: Dseqrecord, source: RestrictionEnzymeDigestionSource) -> tuple[list[Dseqrecord], list[RestrictionEnzymeDigestionSource]]:

    enzyme_list = source.restriction_enzymes
    # If the output is known, we remove empty strings (representing pre-existing edges of the molecule)
    if len(source.fragment_boundaries) > 0:
        enzyme_list = [e for e in source.restriction_enzymes if e != '']
    # enzyme_list = sorted(enzyme_list)
    # TODO: error if enzyme does not exist
    enzymes = RestrictionBatch(first=enzyme_list)

    fragments: list[Dseqrecord] = seq.cut(enzymes)

    # TODO: potential inconsistent behaviour
    # cutsites_positions are the positions of the cut in the forward
    # strand, so they may not match fragment.pos
    # For example in this case cutsite would be 4,
    # but the fragment.pos would be 2
    # NNNN  NN
    # NN  NNNN
    # This may cause problem for overlapping sites, for now I just remove
    # the check
    cutsites, cutsites_positions = get_cutsite_order(seq, enzymes)

    sources = list()

    # We extract the fragment boundaries like this, because the dseqrecord shifts the
    # origin
    fragment_boundaries = [fragment.pos % len(seq) for fragment in seq.seq.cut(enzymes)]

    # Sort fragments by their position in the sequence
    fragments = sort_by(fragments, fragment_boundaries)

    for i, fragment in enumerate(fragments):
        newsource = source.model_copy()
        start = cutsites_positions[i]
        end = (start + len(fragment))
        newsource.fragment_boundaries = [start, end]

        # This check is removed, see TODO above
        # In a circular molecule with a single cut, the length of the fragment is
        # bigger than the length of the molecule (there are overhang of both sides)
        # so it is necessary to stack the sequence three times in order to check
        # if (3 * str(seq.seq))[start:end] != str(fragment.seq):
        #     extra = 'from cutsites: ' + (3 * str(seq.seq))[start:end] + '\nfrom fragment: ' + str(fragment.seq)
        #     raise ValueError(f'Something is wrong with cutsite processing:\n{extra}')
        newsource.restriction_enzymes = [cutsites[i], cutsites[i + 1]]
        sources.append(newsource)

        # Name the fragment:
        fragment.annotations['keywords'] = f'{seq.name}_{newsource.restriction_enzymes[0]}[{start}-{end}]{newsource.restriction_enzymes[1]}'

    return fragments, sources


def get_pcr_products_list(template: Dseqrecord, source: PCRSource, primers: list[Primer], minimal_annealing: int) -> tuple[list[Dseqrecord], list[PCRSource]]:
    anneal = Anneal(primers, template, limit=minimal_annealing)
    amplicon: Amplicon
    sources = list()

    for amplicon in anneal.products:
        new_source = source.model_copy()
        new_source.primers = [int(amplicon.forward_primer.id), int(amplicon.reverse_primer.id)]
        new_source.fragment_boundaries = [amplicon.forward_primer.position, amplicon.reverse_primer.position]
        new_source.primer_footprints = [amplicon.forward_primer._fp, amplicon.reverse_primer._fp]
        sources.append(new_source)

    return anneal.products, sources


def assembly_is_duplicate(assembly: StickyLigationSource) -> bool:
    """
    For linear assemblies we apply the constrain that first fragment is not inverted.
    For circular assemblies, we apply that constrain, plus that the smallest id comes first.
    """
    if assembly.circularised:
        return assembly.fragments_inverted[0] or min(assembly.input) != assembly.input[0]
    else:
        return assembly.fragments_inverted[0]


def get_assembly_list_from_sticky_ligation_source(seqs: list[Dseqrecord], source: StickyLigationSource) -> list[Dseqrecord]:

    assembly: list[Dseqrecord] = list()
    for index, id in enumerate(source.input):
        # Find the element in the list that has that id and add it to the assembly
        assembly.append(next(seq for seq in seqs if int(seq.id) == id))

        # Invert it necessary
        if source.fragments_inverted[index]:
            assembly[-1] = assembly[-1].reverse_complement()
    return assembly


def perform_assembly(assembly: tuple[Dseqrecord], circularise) -> Dseqrecord:
    out = sum((f for f in assembly), Dseqrecord(''))
    if circularise:
        out = out.looped()
    return out


def assembly_list_is_valid(assembly: tuple[Dseqrecord], circularise=False, partial=False) -> bool:
    for i in range(0, len(assembly) - 1):
        if not sum_is_sticky(assembly[i].seq, assembly[i + 1].seq, partial):
            return False

    if circularise:
        return assembly_list_can_be_circularised(assembly)

    return True


def assembly_list_can_be_circularised(assembly: tuple[Dseqrecord]):
    return sum_is_sticky(assembly[-1].seq, assembly[0].seq)


def get_sticky_ligation_products_list(seqs: list[Dseqrecord]) -> tuple[list[Dseqrecord], list[StickyLigationSource]]:

    sequence_ids = [s.id for s in seqs]

    # We generate all possible combinations of sequence ids, fragment_inverted and circularised
    # without duplicates (see assembly_is_duplicate).

    possible_sources = list()
    for _input in permutations(sequence_ids):
        for fragments_inverted in product([True, False], repeat=len(_input)):
            for circularised in [True, False]:
                possible_source = StickyLigationSource(input=_input, fragments_inverted=fragments_inverted,
                                                       circularised=circularised)
                if assembly_is_duplicate(possible_source):
                    continue
                possible_sources.append(possible_source)

    # Here we filter those that are compatible
    valid_sources = list()
    products = list()
    for source in possible_sources:
        assembly_list = get_assembly_list_from_sticky_ligation_source(seqs, source)
        if assembly_list_is_valid(assembly_list, source.circularised):

            # For now, if the assembly can be circularised, we make it happen
            # (we exclude linear assemblies that could be circularised)
            if not source.circularised and assembly_list_can_be_circularised(assembly_list):
                continue
            valid_sources.append(source)
            products.append(perform_assembly(assembly_list, source.circularised))

    return products, valid_sources


def get_sequences_from_gb_file_url(url: str) -> Dseqrecord:
    # TODO once pydna parse is fixed it should handle urls that point to non-gb files
    resp = requests.get(url)
    if resp.status_code != 200:
        raise HTTPError(url, 404, 'file requested from url not found', 'file requested from url not found', None)
    return pydna_parse(resp.text)


def request_from_addgene(source: RepositoryIdSource) -> tuple[list[Dseqrecord], list[RepositoryIdSource]]:
    # TODO here maybe it would be good to check that the addgeneID still returns the url requested.
    if 'url' in source.info:
        return [get_sequences_from_gb_file_url(source.info['url'])[0]], [source]

    url = f'https://www.addgene.org/{source.repository_id}/sequences/'
    resp = requests.get(url)
    if resp.status_code == 404:
        raise HTTPError(url, 404, 'wrong addgene id', 'wrong addgene id', None)
    soup = BeautifulSoup(resp.content, 'html.parser')
    print(soup)
    sequence_file_url_dict = dict()
    for _type in ['depositor-full', 'depositor-partial', 'addgene-full', 'addgene-partial']:
        sequence_file_url_dict[_type] = []
        if soup.find(id=_type) is not None:
            sequence_file_url_dict[_type] = [a.get('href') for a in soup.find(id=_type).findAll(class_='genbank-file-download')]

    # TODO provide addgene sequencing data supporting the sequence
    # We prefer to return addgene full if both available
    products = list()
    sources = list()
    for _type in ['addgene-full', 'depositor-full']:
        if len(sequence_file_url_dict[_type]) > 0:
            for seq_url in sequence_file_url_dict[_type]:
                new_source = source.model_copy()
                new_source.info = {'url': seq_url, 'type': _type}
                sources.append(new_source)
                # There should be only one sequence
                products.append(get_sequences_from_gb_file_url(seq_url)[0])
            break
    return products, sources


def correct_name(dseq: Dseqrecord):
    # Can set the name from keyword if locus is set to Exported
    if dseq.name.lower() == 'exported' and dseq.locus.lower() == 'exported' and 'keywords' in dseq.annotations:
        dseq.name = dseq.annotations['keywords'][0]


def create_location(start: int, end: int, strand: int, seq_len: int, is_circular: bool) -> Location:

    if not is_circular:
        return SimpleLocation(start, end, strand)

    coords = [start, end]
    if coords[0] < coords[1]:
        return SimpleLocation(*coords, strand)
    return SimpleLocation(coords[0], seq_len, strand) + SimpleLocation(0, coords[1], strand)


def location_sorter(x, y) -> int:
    """
    Sort by start, then length, then strand.

    It's a bit weird for origin-spanning features in circular DNA,
    as the start is always zero, ultimately the reason for this is
    so that the features are always in the same order even if sets
    are used in the function.
    """
    if x.start != y.start:
        return x.start - y.start
    elif x.end != y.end:
        return x.end - y.end
    return x.strand - y.strand


def get_all_regex_feature_edges(pattern: str, seq: str, is_circular: bool) -> list[tuple[int, int]]:

    subject = 2 * seq if is_circular else seq

    compiled_pattern = regex.compile(pattern, regex.IGNORECASE)
    compiled_pattern_rev = regex.compile('(?r)' + pattern, regex.IGNORECASE)

    matches = list(regex.finditer(compiled_pattern, subject, overlapped=True))
    matches += list(regex.finditer(compiled_pattern_rev, subject, overlapped=True))

    # In circular objects we remove the matches that span the sequence more than once: m.end() - m.start() <= len(seq)
    return list(set([(m.start(), m.end()) for m in matches if (m.end() - m.start() <= len(seq))]))


def find_sequence_regex(pattern: str, seq: str, is_circular: bool) -> list[Location]:

    feature_locations = list()

    # Below, we use +1 in "% (len(seq) + 1)" because start / end is a range, for instance for a
    # string of length 6 spanned entirely, start=0 end=6, so if you want to fold you need
    # to add 1 to the end.

    strand = 1
    feature_edges = get_all_regex_feature_edges(pattern, seq, is_circular)
    feature_locations += [create_location(start % len(seq), (end - 1) % len(seq) + 1, strand, len(seq), is_circular) for start, end in feature_edges]

    strand = -1
    feature_edges = get_all_regex_feature_edges(pattern, reverse_complement(seq), is_circular)
    feature_locations += [create_location(len(seq) - ((end - 1) % len(seq) + 1), len(seq) - (start % len(seq)), strand, len(seq), is_circular) for start, end in feature_edges]

    # We return a unique list, cannot use a set because Location is not hashable
    return sorted([x for i, x in enumerate(feature_locations) if x not in feature_locations[:i]], key=cmp_to_key(location_sorter))


def location_edges(location: Location):
    if isinstance(location, SimpleLocation):
        return location.start, location.end
    if isinstance(location, CompoundLocation):
        print(location.parts)
        return location.parts[0].start, location.parts[-1].end


def get_homologous_recombination_locations(template: Dseqrecord, insert: Dseqrecord, minimal_homology) -> list[Location]:
    """Return the locations of the possible homologous recombination sites."""
    template_seq = str(template.seq)
    insert_seq = str(insert.seq)
    regex_pattern = insert_seq[:minimal_homology] + '.*' + insert_seq[-minimal_homology:]
    locations = find_sequence_regex(regex_pattern, template_seq, template.circular)
    return locations


def perform_homologous_recombination(template: Dseqrecord, insert: Dseqrecord, location: Location):
    edges = location_edges(location)
    if template.circular:
        return (template[edges[1]:edges[0]] + insert).looped()
    return template[0:edges[0]] + insert + template[edges[1]:]
