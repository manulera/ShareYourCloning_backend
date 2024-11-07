from Bio.Data.IUPACData import ambiguous_dna_values as _ambiguous_dna_values
from Bio.Seq import reverse_complement
from pydna.dseqrecord import Dseqrecord as _Dseqrecord
import re
import itertools as _itertools

ambiguous_only_dna_values = {**_ambiguous_dna_values}
for normal_base in 'ACGT':
    del ambiguous_only_dna_values[normal_base]


def compute_regex_site(site: str) -> str:
    upper_site = site.upper()
    for k, v in ambiguous_only_dna_values.items():
        if len(v) > 1:
            upper_site = upper_site.replace(k, f"[{''.join(v)}]")

    # Make case insensitive
    upper_site = f'(?i){upper_site}'
    return upper_site


def dseqrecord_finditer(pattern: str, seq: _Dseqrecord) -> list[re.Match]:
    query = str(seq.seq) if not seq.circular else str(seq.seq) * 2
    matches = re.finditer(pattern, query)
    return (m for m in matches if m.start() <= len(seq))


# Taken from gateway_features.txt, shipped with ApE
raw_gateway_sites = {
    'attB1': 'CHWVTWTGTACAAAAAANNNG',
    'attB2': 'CHWVTWTGTACAAGAAANNNG',
    'attB3': 'CHWVTWTGTATAATAAANNNG',
    'attB4': 'CHWVTWTGTATAGAAAANNNG',
    'attB5': 'CHWVTWTGTATACAAAANNNG',
    'attL1': 'VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTGTACAAAAAANNNG',
    'attL2': 'VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTGTACAAGAAANNNG',
    'attL3': 'VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTGTATAATAAANNNG',
    'attL4': 'VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTGTATAGAAAANNNG',
    'attL5': 'VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTGTATACAAAANNNG',
    'attP1': 'VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTGTACAAAAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV',
    'attP2': 'VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTGTACAAGAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV',
    'attP3': 'VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTGTATAATAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV',
    'attP4': 'VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTGTATAGAAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV',
    'attP5': 'VAAWWAWKRWTTTWWTTYGACTGATAGTGACCTGTWCGTYGMAACAVATTGATRAGCAATKMTTTYYTATAWTGHCMASTWTGTATACAAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV',
    'attR1': 'CHWVTWTGTACAAAAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV',
    'attR2': 'CHWVTWTGTACAAGAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV',
    'attR3': 'CHWVTWTGTATAATAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV',
    'attR4': 'CHWVTWTGTATAGAAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV',
    'attR5': 'CHWVTWTGTATACAAAAGYWGARCGAGAARCGTAARRTGATATAAATATCAATATATTAAATTAGAYTTTGCATAAAAAACAGACTACATAATRCTGTAARACACAACATATBCAGTCV',
    'overlap_4': 'twtGTATAGAaaag',
    'overlap_3': 'twtGTATAATaaag',
    'overlap_2': 'twtGTACAAGaaag',
    'overlap_1': 'twtGTACAAAaaag',
}

gateway_sites = {
    k: {
        'forward_regex': compute_regex_site(v),
        'reverse_regex': compute_regex_site(reverse_complement(v)),
        'consensus_sequence': v,
    }
    for k, v in raw_gateway_sites.items()
}


def gateway_overlap(seqx: _Dseqrecord, seqy: _Dseqrecord, type=None):
    """Find gateway overlaps"""
    if type not in ['BP', 'LR']:
        raise ValueError(f'Invalid overlap type: {type}')

    out = list()
    # Iterate over the four possible att sites
    for num in range(1, 5):
        # Iterate over the two possible orientations
        # The sites have to be in the same orientation (fwd + fwd or rev + rev)
        for pattern in ['forward_regex', 'reverse_regex']:
            # The overlap regex is the same for all types
            overlap_regex = gateway_sites[f'overlap_{num}'][pattern]

            # Iterate over pairs B, P and P, B for BP and L, R and R, L for LR
            for site_x, site_y in zip(type, type[::-1]):
                site_x_regex = gateway_sites[f'att{site_x}{num}'][pattern]
                matches_x = list(dseqrecord_finditer(site_x_regex, seqx))
                if len(matches_x) == 0:
                    continue

                site_y_regex = gateway_sites[f'att{site_y}{num}'][pattern]
                matches_y = list(dseqrecord_finditer(site_y_regex, seqy))
                if len(matches_y) == 0:
                    continue

                for match_x, match_y in _itertools.product(matches_x, matches_y):
                    # Find the overlap sequence within each match
                    overlap_x = re.search(overlap_regex, match_x.group())
                    overlap_y = re.search(overlap_regex, match_y.group())

                    # Sanity check
                    assert (
                        overlap_x is not None and overlap_y is not None
                    ), 'Something went wrong, no overlap found within the matches'

                    out.append(
                        (
                            match_x.start() + overlap_x.start(),
                            match_y.start() + overlap_y.start(),
                            len(overlap_x.group()),
                        )
                    )

    return out
