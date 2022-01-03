from fastapi import FastAPI, UploadFile, File, Query
from pydna.dseqrecord import Dseqrecord
from pydantic import conlist, create_model
from pydna.parsers import parse as pydna_parse
from Bio.SeqIO import read as seqio_read
from pydna.genbank import Genbank
from dna_functions import get_restriction_enzyme_products_list
from pydantic_models import GenbankSequence, SequenceEntity, SequenceFileFormat, GenbankIdSource, RestrictionEnzymeDigestionSource, UploadedFileSource
from fastapi.middleware.cors import CORSMiddleware

# Instance of the API object
app = FastAPI()

# Allow CORS
origins = ["http://localhost:3000"]
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


def format_sequence_genbank(seq: Dseqrecord) -> SequenceEntity:
    gb_seq = GenbankSequence(file_content=seq.format('genbank'))
    return SequenceEntity(sequence=gb_seq)


def read_dsrecord_from_json(seq: SequenceEntity) -> Dseqrecord:
    return pydna_parse(seq.sequence.file_content)[0]


@ app.post('/read_from_file', response_model=create_model(
    'UploadedFileResponse',
    source=(UploadedFileSource, ...),
    sequence=(SequenceEntity, None)
))
async def read_from_file(file: UploadFile = File(...),
                         file_format: SequenceFileFormat = Query(
                             None,
                             description='Format of the sequence file. Unless specified, it will be guessed from the extension')
                         ):
    """Return a json sequence from a sequence file
    """

    if file_format is None:
        extension_dict = {
            'gbk': 'genbank',
            'gb': 'genbank',
            'dna': 'snapgene',
            'fasta': 'fasta'
        }

        # We guess the file type from the extension
        file_format = SequenceFileFormat(
            extension_dict[file.filename.split('.')[-1]])

    if file_format in ['fasta', 'genbank']:
        # Read the whole file to a string
        file_content = (await file.read()).decode()
        dseq = pydna_parse(file_content)[0]

    elif file_format == 'snapgene':
        seq = seqio_read(file.file, file_format)
        iscircular = 'topology' in seq.annotations.keys(
        ) and seq.annotations['topology'] == 'circular'
        dseq = Dseqrecord(seq, circular=iscircular)

    output_sequence = format_sequence_genbank(dseq)
    source = UploadedFileSource(
        file_format=file_format, file_name=file.filename, output_list=[output_sequence])

    return {'sequence': output_sequence, 'source': source}


@ app.post('/genebank_id', response_model=create_model(
    'GenbankIdResponse',
    source=(GenbankIdSource, ...),
    sequence=(SequenceEntity, None)
))
async def get_from_genebank_id(source: GenbankIdSource):
    gb = Genbank("example@gmail.com")
    seq = Dseqrecord(gb.nucleotide(source.genbank_id))
    output_sequence = format_sequence_genbank(seq)
    source.output_list = [output_sequence]
    return {'sequence': output_sequence, 'source': source}


@ app.post('/restriction', response_model=create_model(
    'RestrictionEnzymeDigestionResponse',
    source=(RestrictionEnzymeDigestionSource, ...),
    sequence=(SequenceEntity, None)
))
async def restriction(source: RestrictionEnzymeDigestionSource,
                      sequences: conlist(SequenceEntity, min_items=1, max_items=1)):
    dseq = read_dsrecord_from_json(sequences[0])
    products_dseq, fragment_boundaries = get_restriction_enzyme_products_list(
        dseq, source)

    source.output_list = [format_sequence_genbank(seq) for seq
                          in products_dseq]
    source.fragment_boundaries = fragment_boundaries

    if source.output_index is not None:
        output_sequence = source.output_list[source.output_index]
    else:
        output_sequence = None
    return {'sequence': output_sequence, 'source': source}
