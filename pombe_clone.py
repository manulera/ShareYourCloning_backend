# from main import genome_coordinates
from pydantic_models import GenomeCoordinatesSource, TextFileSequence

# from httpx import AsyncClient, Response
# from ncbi_requests import get_annotation_from_query
import asyncio
import json


async def main(gene: str, assembly_accession: str, padding: int = 1000):
    # annotation = await get_annotation_from_query(gene, assembly_accession)
    # gene_range = annotation['genomic_regions'][0]['gene_range']['range'][0]
    # sequence_accession = annotation['genomic_regions'][0]['gene_range']['accession_version']
    # locus_tag = annotation.get('locus_tag', None)
    # gene_id = annotation.get('gene_id', None)
    # start = int(gene_range['begin'])
    # end = int(gene_range['end'])
    # orientation = 1 if gene_range['orientation'] == 'plus' else -1

    # source = GenomeCoordinatesSource(
    #     id=1,
    #     start=start - padding,
    #     end=end + padding,
    #     strand=orientation,
    #     assembly_accession=assembly_accession,
    #     sequence_accession=sequence_accession,
    #     locus_tag=locus_tag,
    #     gene_id=gene_id,
    #     output_name=gene,
    # )
    # locus = await genome_coordinates(source)
    # with open('locus.json', 'w') as f:
    #     locus['sequences'][0] = locus['sequences'][0].model_dump()
    #     locus['sources'][0] = locus['sources'][0].model_dump()
    #     f.write(json.dumps(locus, indent=2))

    # Read the file:
    with open('locus.json', 'r') as f:
        locus = json.load(f)

    locus_seq: TextFileSequence = TextFileSequence.model_validate(locus['sequences'][0])
    locus_seq.id = 2
    locus_source: GenomeCoordinatesSource = GenomeCoordinatesSource.model_validate(locus['sources'][0])
    locus_source.output = 2

    return


if __name__ == '__main__':
    gene = 'SPAPB1A10.09'
    assembly_accession = 'GCF_000002945.2'
    asyncio.run(main(gene, assembly_accession))
