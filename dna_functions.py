from Bio.Restriction.Restriction import CommOnly
from Bio.SeqFeature import FeatureLocation
from pydna.dseqrecord import Dseqrecord
from typing import List, OrderedDict

from pydna.seqfeature import SeqFeature

from pydantic_models import RestrictionEnzymeDigestionSource


def get_restriction_enzyme_products_list(seq: Dseqrecord, source: RestrictionEnzymeDigestionSource) -> tuple[List[Dseqrecord], List[int]]:

    # TODO: error if enzyme does not exist
    enzymes = [CommOnly.format(e) for e in source.restriction_enzymes]

    output_list: List[Dseqrecord] = seq.cut(enzymes)
    fragment_boundaries = [fragment.pos for fragment in seq.seq.cut(enzymes)]
    if not seq.circular:
        fragment_boundaries.append(len(seq))
    else:
        fragment_boundaries.append(fragment_boundaries[0])
    # For now, to represent the overhangs of the enzyme cut, we just add a
    # feature to the Dseqrecord
    for fragment in output_list:
        five_prime_end = fragment.seq.five_prime_end()
        if five_prime_end[0] != 'blunt':
            feature_name = 'overhang_' + five_prime_end[0]
            start = 0
            end = len(five_prime_end[1])
            fragment.features.append(SeqFeature(
                location=FeatureLocation(start, end),
                type="misc_feature",
                qualifiers=OrderedDict({"label": feature_name}), strand=None))

        three_prime_end = fragment.seq.three_prime_end()
        if three_prime_end[0] != 'blunt':
            feature_name = 'overhang_' + three_prime_end[0]
            start = len(fragment) - len(three_prime_end[1])
            end = len(fragment)
            fragment.features.append(SeqFeature(
                location=FeatureLocation(start, end),
                type="misc_feature",
                qualifiers=OrderedDict({"label": feature_name}), strand=None))

    return output_list, fragment_boundaries
