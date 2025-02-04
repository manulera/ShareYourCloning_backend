from fastapi import Query, HTTPException
from Bio.Restriction.Restriction_Dictionary import rest_dict


from ..dna_functions import (
    format_sequence_genbank,
    read_dsrecord_from_json,
)
from ..dna_utils import align_sanger_traces
from ..pydantic_models import (
    TextFileSequence,
    BaseCloningStrategy,
)
from ..get_router import get_router
from ..utils import api_version


router = get_router()


@router.get('/version')
async def get_version():
    return api_version()


@router.get('/restriction_enzyme_list', response_model=dict[str, list[str]])
async def get_restriction_enzyme_list():
    """Return the dictionary of restriction enzymes"""
    return {'enzyme_names': list(rest_dict.keys())}


@router.post(
    '/validate',
    summary='Validate a cloning strategy',
)
async def cloning_strategy_is_valid(
    cloning_strategy: BaseCloningStrategy,
) -> bool:
    """Validate a cloning strategy"""
    return True


@router.post('/rename_sequence', response_model=TextFileSequence)
async def rename_sequence(
    sequence: TextFileSequence,
    name: str = Query(..., description='The new name of the sequence.', pattern=r'^[^\s]+$'),
):
    """Rename a sequence"""
    dseqr = read_dsrecord_from_json(sequence)
    return format_sequence_genbank(dseqr, name)


@router.post('/align_sanger', response_model=list[str])
async def align_sanger(
    sequence: TextFileSequence,
    traces: list[str],
):
    """Align a list of sanger traces to a sequence"""

    dseqr = read_dsrecord_from_json(sequence)
    try:
        return align_sanger_traces(dseqr, traces)
    except RuntimeError as e:
        raise HTTPException(status_code=400, detail=str(e))
