"""
Some extra tests that are not covered by the main endpoint tests
"""

import ncbi_requests


def test_empty_get_assembly_accession_from_sequence_accession():
    # For accessions that are not linked to assemblies
    assert [] == ncbi_requests.get_assembly_accession_from_sequence_accession('DQ208311.2')
