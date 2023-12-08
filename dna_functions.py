from functools import cmp_to_key
from urllib.error import HTTPError
from Bio.Restriction.Restriction import RestrictionBatch
from Bio.Seq import reverse_complement
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from pydantic_models import PCRSource, PrimerModel, RepositoryIdSource, \
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
from pydna.utils import shift_location

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
    gb_seq = GenbankSequence(file_content=seq.format('genbank'),
                             overhang_crick_3prime=seq.seq.ovhg,
                             overhang_watson_3prime=seq.seq.watson_ovhg())
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
        out_dseq_record = Dseqrecord(Dseq.from_full_sequence_and_overhangs(
            str(initial_dseqrecord.seq),
            seq.sequence.overhang_crick_3prime,
            seq.sequence.overhang_watson_3prime
        ), features=initial_dseqrecord.features )
    # We set the id to the integer converted to integer (this is only
    # useful for assemblies)
    out_dseq_record.id = str(seq.id)
    return out_dseq_record



def get_invalid_enzyme_names(enzyme_names_list: list[str|None]) -> list[str]:
    rest_batch = RestrictionBatch()
    invalid_names = list()
    for name in enzyme_names_list:
        # Empty enzyme names are the natural edges of the molecule
        if name is not None:
            try:
                rest_batch.format(name)
            except ValueError:
                invalid_names.append(name)
    return invalid_names




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


def location_sorter(x, y) -> int:
    """
    Sort by start, then length, then strand.
    """
    if x.parts[0].start != y.parts[0].start:
        return x.parts[0].start - y.parts[0].start
    elif x.parts[-1].end != y.parts[-1].end:
        return x.parts[-1].end - y.parts[-1].end
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

    # Strand 1
    feature_edges = get_all_regex_feature_edges(pattern, seq, is_circular)
    # We use shift_location to format origin-spanning features in circular DNA
    feature_locations += [ shift_location(SimpleLocation(start, end, 1), 0, len(seq)) for start, end in feature_edges]

    # Strand -1
    feature_edges = get_all_regex_feature_edges(pattern, reverse_complement(seq), is_circular)
    feature_locations += [ shift_location(SimpleLocation(start, end, 1)._flip(len(seq)), 0, len(seq)) for start, end in feature_edges]

    # We return a unique list, cannot use a set because Location is not hashable
    return sorted([x for i, x in enumerate(feature_locations) if x not in feature_locations[:i]], key=cmp_to_key(location_sorter))


def location_edges(location: Location):
    if isinstance(location, SimpleLocation):
        return location.start, location.end
    if isinstance(location, CompoundLocation):
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

