import requests
from fastapi import HTTPException
from urllib.error import HTTPError
from pydna.genbank import Genbank
import os

NCBI_API_KEY = os.environ.get('NCBI_API_KEY')

if NCBI_API_KEY is None:
    headers = None
else:
    headers = {'api_key': NCBI_API_KEY}


# TODO: this does not return old assembly accessions, see https://github.com/ncbi/datasets/issues/380#issuecomment-2231142816
def get_assembly_accession_from_sequence_accession(sequence_accession: str) -> list[str]:
    """Get the assembly accession from a sequence accession"""

    url = f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/sequence_accession/{sequence_accession}/sequence_assemblies'
    resp = requests.get(url, headers=headers)
    data = resp.json()
    if 'accessions' in data:
        return data['accessions']
    else:
        return []


def get_sequence_accessions_from_assembly_accession(assembly_accession: str) -> list[str]:
    """Get the sequence accessions from an assembly accession"""
    url = f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{assembly_accession}/sequence_reports'
    resp = requests.get(url, headers=headers)
    data = resp.json()
    if 'reports' in data:
        return [report['refseq_accession'] for report in data['reports']] + [
            report['genbank_accession'] for report in data['reports']
        ]
    elif 'total_count' in data:
        raise HTTPException(400, f'No sequence accessions linked, see {url}')
    else:
        raise HTTPException(404, 'Wrong assembly accession number')


def get_annotation_from_locus_tag(locus_tag: str, assembly_accession: str) -> dict:
    url = f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{assembly_accession}/annotation_report?search_text={locus_tag}'
    resp = requests.get(url, headers=headers)
    if resp.status_code == 404:
        raise HTTPException(404, 'wrong accession number')
    data = resp.json()
    if 'reports' not in data:
        raise HTTPException(404, 'wrong locus_tag')

    matching_annotations = list(a['annotation'] for a in data['reports'] if a['annotation']['locus_tag'] == locus_tag)

    if len(matching_annotations) == 0:
        raise HTTPException(404, 'wrong locus_tag')
    elif len(matching_annotations) > 1:
        # Not sure if this can ever happen, but just in case
        raise HTTPException(400, 'multiple matches for locus_tag')

    return matching_annotations[0]


def get_genbank_sequence_subset(sequence_accession, start, end, strand):
    gb_strand = 1 if strand == 1 else 2
    gb = Genbank('example@gmail.com')
    try:
        return gb.nucleotide(sequence_accession, start, end, gb_strand)
    except HTTPError as e:
        if e.code == 400:
            raise HTTPException(404, 'wrong sequence accession')
        else:
            raise HTTPException(503, 'NCBI returned an error')


def validate_coordinates_pre_request(start, end, strand):
    # TODO: move this to the class
    if strand not in [1, -1]:
        raise HTTPException(422, 'strand must be 1 or -1')
    if start >= end:
        raise HTTPException(422, 'start must be less than end')
    if start < 1:
        raise HTTPException(422, 'start must be greater than 0')
