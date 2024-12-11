"""
Some extra tests that are not covered by the main endpoint tests
"""

import shareyourcloning.ncbi_requests as ncbi_requests
import pytest


def test_empty_get_assembly_accession_from_sequence_accession():
    # For accessions that are not linked to assemblies
    assert [] == ncbi_requests.get_assembly_accession_from_sequence_accession('DQ208311.2')


@pytest.mark.xfail(reason='waiting on https://github.com/ncbi/datasets/issues/380#issuecomment-2231142816')
def test_get_assembly_accession_from_sequence_accession():
    # For accessions that are linked to assemblies
    assert ['GCF_000002945.2', 'GCF_000002945.1'] == ncbi_requests.get_assembly_accession_from_sequence_accession(
        'NC_003424.3'
    )
