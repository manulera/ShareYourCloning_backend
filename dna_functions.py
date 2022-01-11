from Bio.Restriction.Restriction import CommOnly
from Bio.Seq import reverse_complement
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from typing import List
from pydantic_models import RestrictionEnzymeDigestionSource,\
    GenbankSequence, SequenceEntity, StickyLigationSource
from pydna.parsers import parse as pydna_parse
from itertools import permutations, product, chain
import numpy


def sum_is_sticky(seq1: Dseq, seq2: Dseq) -> bool:
    """Returns true if the 3' end of seq1 and 5' end of seq2 ends are sticky and compatible for ligation"""
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


def assemblies_are_circular_permutations(assembly_1: StickyLigationSource, assembly_2: StickyLigationSource):
    """Check if lists are circular permutations of one another"""
    # TODO either fix or remove this function. As of now it fails to detect the same circular assembly when all
    # fragments are inverted
    # We multiply the inputs by -1 to indicate the inversion, and then compare circular list
    inputs1 = assembly_1.input
    inverted1 = assembly_1.fragments_inverted
    values1 = [[-inputs1[i], inputs1[i]] if inverted1[i] else [inputs1[i], -inputs1[i]]
               for i in range(len(inputs1))]

    values1 = list(chain.from_iterable(values1))

    inputs2 = assembly_2.input
    inverted2 = assembly_2.fragments_inverted
    values2 = [[-inputs2[i], inputs2[i]] if inverted2[i] else [inputs2[i], -inputs2[i]]
               for i in range(len(inputs2))]
    values2 = list(chain.from_iterable(values2))

    list_len = len(values1)

    # We concatenate a list with itself to compare the circular permutation
    summed_list = values1 + values1
    for i in range(list_len):
        if values2 == summed_list[i:i + list_len]:
            return True
    return False


def eliminate_assembly_duplicates(assemblies_in: List[StickyLigationSource],
                                  products_in: List[Dseqrecord]) -> \
        tuple[List[StickyLigationSource], List[Dseqrecord]]:
    """Circular assemblies can be equivalent: [1,2,3] == [2,3,1] == [3,1,2]"""
    assemblies_out = list()
    products_out = list()
    while(len(assemblies_in)):
        assembly = assemblies_in.pop()
        product = products_in.pop()
        # For the linear assemblies, we apply the constrain that the
        # first fragment is not inverted
        # For the circular assemblies, we apply the constrain that the
        # smallest id comes first, and that is not inverted. This eliminates
        # all equivalent assemblies.
        if not assembly.fragments_inverted[0]:
            if not assembly.circularised or (min(assembly.input) == assembly.input[0]):
                assemblies_out.append(assembly)
                products_out.append(product)

    return assemblies_out, products_out


def get_assembly_summary_from_fragment_list(assembly: tuple[Dseqrecord]) -> StickyLigationSource:
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


def sum_assembly_fragments(assembly: tuple[Dseqrecord]) -> Dseqrecord:
    if len(assembly) == 1:
        return assembly[0]
    else:
        out = assembly[0]
        for f in assembly[1:]:
            out += f
        return out


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
            assembly_is_valid = True
            for i in range(0, len(assembly) - 1):
                if not sum_is_sticky(assembly[i].seq, assembly[i + 1].seq):
                    assembly_is_valid = False
                    break
            if assembly_is_valid:

                linear_ligation = sum_assembly_fragments(assembly)
                assembly_summary = get_assembly_summary_from_fragment_list(
                    assembly)
                # Sometimes they can be circularised, for now, we circularise by
                # default if possible
                if sum_is_sticky(linear_ligation.seq, linear_ligation.seq):
                    possible_products.append(linear_ligation.looped())
                    assembly_summary.circularised = True
                else:
                    possible_products.append(linear_ligation)
                    assembly_summary.circularised = False
                possible_assemblies.append(assembly_summary)
    possible_assemblies, possible_products = eliminate_assembly_duplicates(
        possible_assemblies, possible_products)
    return possible_products, possible_assemblies
