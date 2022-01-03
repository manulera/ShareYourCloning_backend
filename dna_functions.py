from Bio.Restriction.Restriction import CommOnly
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from typing import List
from pydantic_models import RestrictionEnzymeDigestionSource, GenbankSequence, SequenceEntity
from pydna.parsers import parse as pydna_parse


def format_sequence_genbank(seq: Dseqrecord) -> SequenceEntity:
    overhang_crick_3prime, overhang_watson_3prime = both_overhangs_from_dseq(
        seq.seq)
    gb_seq = GenbankSequence(file_content=seq.format('genbank'),
                             overhang_crick_3prime=overhang_crick_3prime,
                             overhang_watson_3prime=overhang_watson_3prime)
    return SequenceEntity(sequence=gb_seq)


def read_dsrecord_from_json(seq: SequenceEntity) -> Dseqrecord:
    initial_dseqrecord: Dseqrecord = pydna_parse(seq.sequence.file_content)[0]
    if seq.sequence.overhang_watson_3prime == 0 and seq.sequence.overhang_crick_3prime == 0:
        return initial_dseqrecord
    else:
        dseqrecord_with_overhang = Dseqrecord(dseq_from_both_overhangs(
            str(initial_dseqrecord.seq),
            seq.sequence.overhang_crick_3prime,
            seq.sequence.overhang_watson_3prime))
        dseqrecord_with_overhang.features = initial_dseqrecord.features
        return dseqrecord_with_overhang


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
