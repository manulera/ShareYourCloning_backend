from fastapi import Query, HTTPException
from pydantic import create_model
from urllib.error import HTTPError
import os

from ..get_router import get_router
from ..pydantic_models import TextFileSequence, AnnotationSource
from ..dna_functions import (
    read_dsrecord_from_json,
    annotate_with_plannotate as _annotate_with_plannotate,
    format_sequence_genbank,
)
from ..gateway import find_gateway_sites

router = get_router()

PLANNOTATE_URL = os.environ['PLANNOTATE_URL'] if 'PLANNOTATE_URL' in os.environ else None
PLANNOTATE_TIMEOUT = int(os.environ['PLANNOTATE_TIMEOUT']) if 'PLANNOTATE_TIMEOUT' in os.environ else 20

# Handle trailing slash:
if PLANNOTATE_URL is not None and not PLANNOTATE_URL.endswith('/'):
    PLANNOTATE_URL += '/'


@router.post('/annotation/get_gateway_sites', response_model=dict[str, list[str]])
async def get_gateway_sites(
    sequence: TextFileSequence, greedy: bool = Query(False, description='Whether to use the greedy algorithm.')
) -> dict[str, list[str]]:
    """
    Get a dictionary with the names of the gateway sites present in the sequence and their locations as strings.
    """
    dseqr = read_dsrecord_from_json(sequence)
    sites_dict = find_gateway_sites(dseqr, greedy)
    for site in sites_dict:
        sites_dict[site] = [str(loc) for loc in sites_dict[site]]
    return sites_dict


if PLANNOTATE_URL is not None:

    @router.post(
        '/annotate/plannotate',
        summary='Annotate a sequence with Plannotate',
        response_model=create_model(
            'PlannotateResponse',
            sources=(list[AnnotationSource], ...),
            sequences=(list[TextFileSequence], ...),
        ),
    )
    async def annotate_with_plannotate(
        sequence: TextFileSequence,
        source: AnnotationSource,
    ):
        input_seqr = read_dsrecord_from_json(sequence)
        # Make a request submitting sequence as a file:
        try:
            seqr, annotations, version = await _annotate_with_plannotate(
                sequence.file_content, f'{sequence.id}.gb', PLANNOTATE_URL + 'annotate', PLANNOTATE_TIMEOUT
            )
        except HTTPError as e:
            raise HTTPException(e.code, e.msg) from e

        source.annotation_report = annotations
        source.annotation_tool = 'plannotate'
        source.annotation_tool_version = version
        seqr.name = input_seqr.name + '_annotated'

        return {'sources': [source], 'sequences': [format_sequence_genbank(seqr, source.output_name)]}
