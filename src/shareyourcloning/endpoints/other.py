from fastapi import Query
from Bio.Restriction.Restriction_Dictionary import rest_dict
import tempfile
import os
from Bio import SeqIO


from ..dna_functions import (
    format_sequence_genbank,
    read_dsrecord_from_json,
)
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


@router.post('/align', response_model=list[str])
async def align_sequences(
    sequences: list[str],
):
    """Align a list of sequences"""
    # Create temporary directory and file
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, 'sequences.fa')
        aln_path = os.path.join(tmpdir, 'aligned.fa')

        # Write sequences to FASTA file
        with open(fasta_path, 'w') as f:
            for i, seq in enumerate(sequences):
                f.write(f">seq{i+1}\n{seq}\n")

        # Run clustalo alignment
        os.system(f"clustalo -i {fasta_path} --force --outfmt=fa -o {aln_path}")

        # Read the file and return the sequences
        with open(aln_path, 'r') as f:
            records = list(SeqIO.parse(f, 'fasta'))

        return [str(record.seq) for record in records]
