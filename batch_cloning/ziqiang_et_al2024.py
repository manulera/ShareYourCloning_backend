from Bio.Seq import reverse_complement
from fastapi import APIRouter
from fastapi.responses import FileResponse
from pathlib import Path


router = APIRouter()


def validate_protospacers(protospacers):
    """
    Ensure that the protospacers are compatible for the golden gate assembly, by checking that the central 4
    bps are different in all of them, and also different from the other assembly joints.
    """
    forbidden_joints = ['AACA', 'GCTT']

    for i, ps in enumerate(protospacers):
        # Must be only ACGT
        if not set(ps).issubset({'A', 'C', 'G', 'T'}):
            raise ValueError(f"Protospacer {i} {ps} contains invalid bases")
        if len(ps) != 20:
            raise ValueError(f"Protospacer {i} {ps} is not 20 bp long")
        if ps[8:12] in forbidden_joints:
            raise ValueError(
                f"Protospacer {i} has a forbidden joint {ps[8:12]}, which is used for the constant parts of the assembly"
            )
        # Find any other protospacers with the same joint sequence
        same_joint_indexes = [
            j + 1 for j, other_ps in enumerate(protospacers) if i != j and ps[8:12] == other_ps[8:12]
        ]
        if same_joint_indexes:
            raise ValueError(
                f"Protospacer {i + 1} has the same joint as protospacers {', '.join(map(str, same_joint_indexes))}"
            )


def design_primers(protospacers):

    fwd_prefix = 'taggtctcc'
    rvs_prefix = 'atggtctca'
    fwd_suffix = 'gttttagagctagaa'
    rvs_suffix = 'tgcaccagccgggaa'

    primers = list()
    for protospacer in protospacers:
        primers.append(rvs_prefix + reverse_complement(protospacer[:12]) + rvs_suffix)
        primers.append(fwd_prefix + protospacer[8:] + fwd_suffix)
    return primers


@router.get('/batch_cloning/ziqiang_et_al2024')
def ziqiang_et_al2024():
    return FileResponse(Path(__file__).parent / 'ziqiang_et_al2024.html')
