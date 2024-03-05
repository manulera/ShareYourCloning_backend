import requests
from json import JSONDecodeError
from fastapi import HTTPException
from urllib.error import HTTPError
from pydna.genbank import Genbank


def get_assembly_accession_from_sequence_accession(sequence_accession: str) -> str | None:
    """Get the assembly accession from a sequence accession"""
    url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=assembly&id={sequence_accession}&idtype=acc&retmode=json'
    resp = requests.get(url)
    try:
        data = resp.json()
    except JSONDecodeError:
        raise HTTPException(404, 'sequence accession not found')
    assembly_accession = None
    # sequence_accessions is associated with any assembly
    if 'error' in data:
        raise HTTPException(503, f'failed to access NCBI: {data["error"]}')
    if 'linksetdbs' in data['linksets'][0]:
        links = assembly_id = data['linksets'][0]['linksetdbs'][0]['links']
        if len(links) > 1:
            # raise not implemented
            return None
        assembly_id = links[0]
        url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&retmode=json&id={assembly_id}'
        resp = requests.get(url)
        try:
            data = resp.json()
        except JSONDecodeError:
            raise HTTPException(503, 'Error accessing assembly id')
        if 'error' in data:
            raise HTTPException(503, f'failed to access NCBI: {data["error"]}')
        assembly_accession = data['result'][assembly_id]['assemblyaccession']
    return assembly_accession


def get_annotation_from_locus_tag(locus_tag: str, assembly_accession: str) -> dict:
    url = f'https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/{assembly_accession}/annotation_report?search_text={locus_tag}'
    resp = requests.get(url)
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


def get_genbank_sequence_subset(sequence_accession, start, stop, strand):
    gb_strand = 1 if strand == 1 else 2
    gb = Genbank("example@gmail.com")
    try:
        return gb.nucleotide(sequence_accession, start, stop, gb_strand)
    except HTTPError as e:
        if e.code == 400:
            raise HTTPException(404, 'wrong sequence accession')
        else:
            raise HTTPException(503, 'NCBI returned an error')


def validate_coordinates_pre_request(start, stop, strand):
    # TODO: move this to the class
    if strand not in [1, -1]:
        raise HTTPException(422, 'strand must be 1 or -1')
    if start >= stop:
        raise HTTPException(422, 'start must be less than stop')
    if start < 1:
        raise HTTPException(422, 'start must be greater than 0')
