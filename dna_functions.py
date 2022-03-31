
from Bio.Restriction.Restriction import RestrictionBatch, RestrictionType
from Bio.Seq import reverse_complement
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from pydantic_models import PCRSource, PrimerModel, RestrictionEnzymeDigestionSource,\
    GenbankSequence, SequenceEntity, StickyLigationSource
from pydna.parsers import parse as pydna_parse
from itertools import permutations, product
from pydna.primer import Primer
from pydna.amplify import Anneal
from pydna.amplicon import Amplicon


def sum_is_sticky(seq1: Dseq, seq2: Dseq) -> bool:
    """Return true if the 3' end of seq1 and 5' end of seq2 ends are sticky and compatible for ligation."""
    # TODO: add support for restriction sites that partially overlap
    type_seq1, sticky_seq1 = seq1.three_prime_end()
    type_seq2, sticky_seq2 = seq2.five_prime_end()
    return 'blunt' != type_seq2 and type_seq2 == type_seq1 and str(sticky_seq2) == str(reverse_complement(sticky_seq1))


def format_sequence_genbank(seq: Dseqrecord) -> SequenceEntity:

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


def get_restriction_enzyme_products_list(seq: Dseqrecord, source: RestrictionEnzymeDigestionSource) -> tuple[list[Dseqrecord], list[RestrictionEnzymeDigestionSource]]:

    # The output is known, so we remove empty strings
    enzyme_list = source.restriction_enzymes
    if len(source.fragment_boundaries) > 0:
        enzyme_list = [e for e in source.restriction_enzymes if e != '']
    # enzyme_list = sorted(enzyme_list)
    # TODO: error if enzyme does not exist
    enzymes = RestrictionBatch(first=enzyme_list)

    fragments: list[Dseqrecord] = seq.cut(enzymes)
    cutsites, cutsites_positions = get_cutsite_order(seq, enzymes)

    sources = list()

    # We extract the fragment boundaries like this, because the dseqrecord shifts the
    # origin
    fragment_boundaries = [fragment.pos % len(seq) for fragment in seq.seq.cut(enzymes)]

    # Sort fragments by their position in the sequence
    fragments = sort_by(fragments, fragment_boundaries)

    for i, fragment in enumerate(fragments):
        newsource = source.copy()
        start = cutsites_positions[i]
        end = (start + len(fragment))
        newsource.fragment_boundaries = [start, end]

        # In a circular molecule with a single cut, the length of the fragment is
        # bigger than the length of the molecule (there are overhang of both sides)
        # so it is necessary to stack the sequence three times in order to check

        if (3 * str(seq.seq))[start:end] != str(fragment.seq):
            extra = 'from cutsites: ' + (3 * str(seq.seq))[start:end] + '\nfrom fragment: ' + str(fragment.seq)
            print(fragment_boundaries)
            print([f.pos for f in seq.seq.cut(enzymes)])
            print(cutsites, cutsites_positions)
            print(seq.seq)
            print(extra)
            exit()
            raise ValueError(f'Something is wrong with cutsite processing:\n{extra}')
        newsource.restriction_enzymes = [cutsites[i], cutsites[i + 1]]
        sources.append(newsource)

    return fragments, sources


def get_pcr_products_list(template: Dseqrecord, source: PCRSource, primers: list[Primer]) -> tuple[list[Dseqrecord], list[PCRSource]]:
    anneal = Anneal(primers, template, limit=source.primer_annealing_settings.minimum_annealing)
    amplicon: Amplicon
    sources = list()

    for amplicon in anneal.products:
        new_source = source.copy()
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
    if not assembly.fragments_inverted[0]:
        if not assembly.circularised or (min(assembly.input) == assembly.input[0]):
            return True
    return False


def get_sticky_ligation_source_from_assembly_list(assembly: tuple[Dseqrecord]) -> StickyLigationSource:
    fragments_inverted = list()
    fragments_order = list()
    for i in range(len(assembly)):
        fragments_inverted.append(
            assembly[i].id.endswith('_rc'))
        fragments_order.append(
            int(assembly[i].id[:-3]) if fragments_inverted[i]
            else int(assembly[i].id)
        )
        assembly_summary = StickyLigationSource(
            input=fragments_order,
            fragments_inverted=fragments_inverted)
    return assembly_summary


def get_assembly_list_from_sticky_ligation_source(seqs: list[Dseqrecord], source: StickyLigationSource) -> list[Dseqrecord]:

    assembly: list[Dseqrecord] = list()
    for index, id in enumerate(source.input):
        # Find the element in the list that has that id and add it to the assembly
        assembly.append(next(seq for seq in seqs if int(seq.id) == id))

        # Invert it necessary
        if source.fragments_inverted[index]:
            assembly[-1].reverse_complement()
    return assembly


def perform_assembly(assembly: tuple[Dseqrecord], circularise) -> Dseqrecord:
    out = sum((f for f in assembly), Dseqrecord(''))
    if circularise:
        out = out.looped()
    return out


def assembly_is_valid(assembly: tuple[Dseqrecord], circularise=False) -> bool:
    for i in range(0, len(assembly) - 1):
        if not sum_is_sticky(assembly[i].seq, assembly[i + 1].seq):
            return False

    if circularise:
        return sum_is_sticky(assembly[-1].seq, assembly[0].seq)

    return True


def get_sticky_ligation_products_list(seqs: list[Dseqrecord]) -> tuple[list[Dseqrecord], list[StickyLigationSource]]:

    # TODO: include also partial ligations, it could also be made more performant by creating
    # a class with only the minimal information, rather than reverse-complementing everything,
    # but this is not a priority, I would say.
    possible_products = list()
    possible_assemblies: list[StickyLigationSource] = list()

    # We try all permutations of fragments
    for perm in permutations(seqs):
        # We generate all possible arrangements of fragments in this particular order
        arrangements = [(fragment, fragment.reverse_complement())
                        for fragment in perm]

        assembly: tuple[Dseqrecord]
        for assembly in product(*arrangements):
            if assembly_is_valid(assembly):

                source = get_sticky_ligation_source_from_assembly_list(
                    assembly)

                # Sometimes they can be circularised, for now, we circularise by
                # default if possible
                source.circularised = sum_is_sticky(assembly[-1].seq, assembly[0].seq)

                if assembly_is_duplicate(source):
                    continue

                possible_assemblies.append(source)
                possible_products.append(perform_assembly(assembly, source.circularised))

    return possible_products, possible_assemblies
