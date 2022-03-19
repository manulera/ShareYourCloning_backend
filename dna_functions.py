from Bio.Restriction.Restriction import CommOnly
from Bio.Seq import reverse_complement
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from typing import List
from pydantic_models import RestrictionEnzymeDigestionSource,\
    GenbankSequence, SequenceEntity, StickyLigationSource
from pydna.parsers import parse as pydna_parse
from itertools import permutations, product


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


def get_restriction_enzyme_products_list(seq: Dseqrecord, source: RestrictionEnzymeDigestionSource) -> tuple[List[Dseqrecord], List[int]]:

    # TODO: error if enzyme does not exist
    enzymes = [CommOnly.format(e) for e in source.restriction_enzymes]

    output_list: List[Dseqrecord] = seq.cut(enzymes)
    fragment_boundaries = [fragment.pos for fragment in seq.seq.cut(enzymes)]
    if not seq.circular:
        fragment_boundaries.append(len(seq))
    else:
        fragment_boundaries.append(fragment_boundaries[0])

    return output_list, fragment_boundaries


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


def get_assembly_list_from_sticky_ligation_source(seqs: List[Dseqrecord], source: StickyLigationSource) -> List[Dseqrecord]:

    assembly: List[Dseqrecord] = list()
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


def assembly_is_valid(assembly: tuple[Dseqrecord]) -> bool:
    for i in range(0, len(assembly) - 1):
        if not sum_is_sticky(assembly[i].seq, assembly[i + 1].seq):
            return False
    return True


def get_sticky_ligation_products_list(seqs: List[Dseqrecord]) -> tuple[List[Dseqrecord], List[StickyLigationSource]]:

    # TODO: include also partial ligations, it could also be made more performant by creating
    # a class with only the minimal information, rather than reverse-complementing everything,
    # but this is not a priority, I would say.
    possible_products = list()
    possible_assemblies: List[StickyLigationSource] = list()

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
