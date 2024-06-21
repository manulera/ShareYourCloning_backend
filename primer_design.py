from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design
from Bio.SeqFeature import SimpleLocation
from pydna.utils import locations_overlap, shift_location, location_boundaries


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

    fragment2amplify = pcr_loc.extract(pcr_seq)
    amplicon = primer_design(fragment2amplify, limit=minimal_hybridization_length, target_tm=target_tm)

    if insert_forward:
        fwd_primer, rvs_primer = amplicon.primers()
    else:
        rvs_primer, fwd_primer = amplicon.primers()

    if fwd_primer is None or rvs_primer is None:
        raise ValueError('Primers could not be designed, try changing settings.')

    hr_loc_start, hr_loc_end = location_boundaries(hr_loc)
    fwd_homology_start = hr_loc_start - homology_length
    rvs_homology_end = hr_loc_end + homology_length

    if not hr_seq.circular:
        if fwd_homology_start < 0:
            raise ValueError('Forward homology region is out of bounds.')
        if rvs_homology_end > len(hr_seq):
            raise ValueError('Reverse homology region is out of bounds.')

    # Convert to locations
    fwd_arm = shift_location(SimpleLocation(fwd_homology_start, hr_loc_start), 0, len(hr_seq))
    rvs_arm = shift_location(SimpleLocation(hr_loc_end, rvs_homology_end), 0, len(hr_seq))

    if locations_overlap(fwd_arm, rvs_arm, len(hr_seq)):
        raise ValueError('Homology arms overlap.')

    fwd_homology = fwd_arm.extract(hr_seq)
    rvs_homology = rvs_arm.extract(hr_seq).reverse_complement()

    return str(fwd_homology.seq).lower() + str(fwd_primer.seq), str(rvs_homology.seq).lower() + str(rvs_primer.seq)
