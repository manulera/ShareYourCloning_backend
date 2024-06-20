from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design
from Bio.SeqFeature import SimpleLocation
from pydna.utils import location_boundaries


def homologous_recombination_primers(
    pcr_seq: Dseqrecord,
    pcr_loc: SimpleLocation,
    hr_seq: Dseqrecord,
    hr_loc: SimpleLocation,
    homology_length: int,
    minimal_hybridization_length: int,
    insert_forward: bool,
    target_tm: float,
) -> tuple[str, str]:

    if hr_seq.circular:
        # TODO: Conditions will have to be adapted for circular sequences, e.g. fwd_homology_start < 0
        # is not enough
        raise ValueError('Homology target sequence must be linear.')
    start, end = location_boundaries(pcr_loc)
    fragment2amplify = pcr_seq[start:end]
    amplicon = primer_design(fragment2amplify, limit=minimal_hybridization_length, target_tm=target_tm)

    if insert_forward:
        fwd_primer, rvs_primer = amplicon.primers()
    else:
        rvs_primer, fwd_primer = amplicon.primers()

    if fwd_primer is None or rvs_primer is None:
        raise ValueError('Primers could not be designed, try changing settings.')

    fwd_homology_start = hr_loc.start - homology_length
    if fwd_homology_start < 0:
        raise ValueError('Forward homology region is out of bounds.')
    fwd_homology = hr_seq[fwd_homology_start : hr_loc.start]

    rvs_homology_end = hr_loc.end + homology_length
    if rvs_homology_end > len(hr_seq):
        raise ValueError('Reverse homology region is out of bounds.')
    rvs_homology = hr_seq[hr_loc.end : rvs_homology_end].reverse_complement()

    return str(fwd_homology.seq).lower() + str(fwd_primer.seq), str(rvs_homology.seq).lower() + str(rvs_primer.seq)
