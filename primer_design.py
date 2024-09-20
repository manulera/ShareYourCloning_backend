from pydna.dseqrecord import Dseqrecord
from pydna.design import primer_design, assembly_fragments
from Bio.SeqFeature import SimpleLocation
from pydna.utils import locations_overlap, shift_location, location_boundaries
from pydna.amplicon import Amplicon
from pydantic_models import PrimerModel


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


def gibson_assembly_primers(
    templates: list[Dseqrecord],
    homology_length: int,
    minimal_hybridization_length: int,
    target_tm: float,
    circular: bool,
) -> list[PrimerModel]:

    initial_amplicons = [
        primer_design(template, limit=minimal_hybridization_length, target_tm=target_tm) for template in templates
    ]
    if circular:
        initial_amplicons.append(initial_amplicons[0])

    templates_with_no_primers = []
    for i, amplicon in enumerate(initial_amplicons):
        if None in amplicon.primers():
            # In circular cases
            index = i + 1 if i < len(templates) else 1
            templates_with_no_primers.append(index)
    templates_with_no_primers = list(sorted(set(templates_with_no_primers)))

    if templates_with_no_primers:
        raise ValueError(
            f'Primers could not be designed for template {", ".join(map(str, templates_with_no_primers))}, try changing settings.'
        )

    assembly_amplicons: list[Amplicon] = assembly_fragments(initial_amplicons, overlap=homology_length)

    all_primers = sum((list(amplicon.primers()) for amplicon in assembly_amplicons), [])
    if circular:
        all_primers[0] = all_primers[-2]
        all_primers = all_primers[:-2]

    for i in range(0, len(all_primers), 2):
        fwd, rvs = all_primers[i : i + 2]
        template = templates[i // 2]
        template_name = template.name if template.name != 'name' else f'seq_{template.id}'
        fwd.name = f'{template_name}_fwd'
        rvs.name = f'{template_name}_rvs'

    return [PrimerModel(id=0, name=primer.name, sequence=str(primer.seq)) for primer in all_primers]
