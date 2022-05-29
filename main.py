from fastapi import FastAPI, UploadFile, File, Query, HTTPException
from pydna.dseqrecord import Dseqrecord
from pydantic import conlist, create_model
from pydna.parsers import parse as pydna_parse
from Bio.SeqIO import read as seqio_read
from pydna.genbank import Genbank
from dna_functions import assembly_list_is_valid, get_assembly_list_from_sticky_ligation_source, get_pcr_products_list, get_restriction_enzyme_products_list, \
    format_sequence_genbank, get_sticky_ligation_products_list, perform_assembly, \
    read_dsrecord_from_json, read_primer_from_json
from pydantic_models import PCRSource, PrimerAnnealingSettings, PrimerModel, SequenceEntity, SequenceFileFormat, \
    GenbankIdSource, RestrictionEnzymeDigestionSource, StickyLigationSource,\
    UploadedFileSource
from fastapi.middleware.cors import CORSMiddleware
from urllib.error import HTTPError, URLError
from fastapi.responses import HTMLResponse

# Instance of the API object
app = FastAPI()

# Allow CORS
# TODO put a wildcard on the shareyourcloning.netlify to
# allow for the draft websites to also work in netlify.
# TODO make this conditional to dev / prod using settings

origins = ["http://localhost:3000", "https://shareyourcloning.netlify.app"]
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# TODO limit the maximum size of submitted files


@app.get('/')
async def greeting():
    html_content = """
        <html>
            <head>
                <title>Welcome to ShareYourCloning API</title>
            </head>
            <body>
                <h1>Welcome to ShareYourCloning API</h1>
                <p>Visit /docs to see the documentation</p>
            </body>
        </html>
        """
    return HTMLResponse(content=html_content, status_code=200)


@ app.post('/read_from_file', response_model=create_model(
    'UploadedFileResponse',
    sources=(list[UploadedFileSource], ...),
    sequences=(list[SequenceEntity], ...)
))
async def read_from_file(file: UploadFile = File(...),
                         file_format: SequenceFileFormat = Query(
                             None,
                             description='Format of the sequence file. \
                                Unless specified, it will be guessed\
                                from the extension'),
                         index_in_file: int = Query(None, description='The index\
                             of the sequence in the file for multi-sequence files')
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
        extension = file.filename.split('.')[-1]
        if extension not in extension_dict:
            raise HTTPException(
                422, 'We could not guess the format of the file from its extension.\
                     Please provide file_format as a query parameter.')

        # We guess the file type from the extension
        file_format = SequenceFileFormat(extension_dict[extension])

    if file_format in ['fasta', 'genbank']:
        # Read the whole file to a string
        file_content = (await file.read()).decode()
        dseqs = pydna_parse(file_content)
        if len(dseqs) == 0:
            raise HTTPException(
                422, 'Pydna parser reader cannot process this file.')

    elif file_format == 'snapgene':
        try:
            seq = seqio_read(file.file, file_format)
        except ValueError:
            raise HTTPException(
                422, 'Biopython snapgene reader cannot process this file.')

        iscircular = 'topology' in seq.annotations.keys(
        ) and seq.annotations['topology'] == 'circular'
        dseqs = [Dseqrecord(seq, circular=iscircular)]

    # The common part
    parent_source = UploadedFileSource(file_format=file_format, file_name=file.filename)
    out_sources = list()
    for i in range(len(dseqs)):
        new_source = parent_source.copy()
        new_source.index_in_file = i
        out_sources.append(new_source)

    out_sequences = [format_sequence_genbank(s) for s in dseqs]

    return {'sequences': out_sequences, 'sources': out_sources}

# TODO: a bit inconsistent that here you don't put {source: {...}} in the request, but
# directly the object.


@ app.post('/genbank_id', response_model=create_model(
    'GenbankIdResponse',
    sources=(list[GenbankIdSource], ...),
    sequences=(list[SequenceEntity], ...)
))
async def get_from_genbank_id(source: GenbankIdSource):
    gb = Genbank("example@gmail.com")
    try:
        seq = Dseqrecord(gb.nucleotide(source.genbank_id))
    except HTTPError as exception:
        if exception.code == 500:
            raise HTTPException(
                503, f'GenBank returned: {exception} - GenBank might be down')
        elif exception.code == 400:
            raise HTTPException(
                404, f'GenBank returned: {exception} - Likely you inserted\
                     a wrong GenBank id')
    except URLError as exception:
        raise HTTPException(504, f'Unable to connect to GenBank: {exception}')

    output_sequence = format_sequence_genbank(seq)

    return {'sequences': [output_sequence], 'sources': [source]}


@ app.post('/restriction', response_model=create_model(
    'RestrictionEnzymeDigestionResponse',
    sources=(list[RestrictionEnzymeDigestionSource], ...),
    sequences=(list[SequenceEntity], ...)
))
async def restriction(source: RestrictionEnzymeDigestionSource,
                      sequences: conlist(SequenceEntity, min_items=1, max_items=1)):
    dseq = read_dsrecord_from_json(sequences[0])
    # TODO: return error if the id of the sequence does not correspond
    # TODO: issue warning if number of enzymes>2 or if one does not cut

    # If the request provides the fragment_boundaries, the program should return only one output.
    output_is_known = False
    if len(source.fragment_boundaries) > 0:
        if len(source.fragment_boundaries) != 2 or len(source.restriction_enzymes) != 2:
            raise HTTPException(
                400, 'If `fragment_boundaries` are provided, the length of `fragment_boundaries` and `restriction_enzymes` must be 2.')
        output_is_known = True

    fragments, out_sources = get_restriction_enzyme_products_list(
        dseq, source)

    out_sequences = [format_sequence_genbank(seq) for seq
                     in fragments]

    if len(out_sequences) == 0:
        raise HTTPException(
            400, 'The enzymes do not cut.')

    # If the user has provided boundaries, we verify that they are correct, and return only those as the output
    if output_is_known:
        for i, out_source in enumerate(out_sources):
            if out_source == source:
                return {'sequences': [out_sequences[i]], 'sources': [out_source]}
        # If we don't find it, there was a mistake
        raise HTTPException(
            400, 'The fragment boundaries / enzymes provided do not correspond to the ones predicted.')

    return {'sources': out_sources, 'sequences': out_sequences}


@ app.post('/sticky_ligation', response_model=create_model(
    'StickyLigationResponse',
    sources=(list[StickyLigationSource], ...),
    sequences=(list[SequenceEntity], ...)
))
async def sticky_ligation(source: StickyLigationSource,
                          sequences: conlist(SequenceEntity, min_items=1)):
    dseqs = [read_dsrecord_from_json(seq) for seq in sequences]
    if len(source.fragments_inverted) > 0:
        # TODO Error if the list has different order or the ids are wrong.
        # TODO check input for unique ids
        assembly_list = get_assembly_list_from_sticky_ligation_source(dseqs, source)
        if not assembly_list_is_valid(assembly_list, source.circularised):
            raise HTTPException(
                400, 'Fragments are not compatible for sticky ligation')
        output_sequence = format_sequence_genbank(perform_assembly(assembly_list, source.circularised))

        return {'sequences': [output_sequence], 'sources': [source]}

    products_dseq, out_sources = get_sticky_ligation_products_list(dseqs)
    if len(products_dseq) == 0:
        raise HTTPException(
            400, 'No combination of these fragments is compatible for sticky ligation')

    out_sequences = [format_sequence_genbank(seq) for seq
                     in products_dseq]

    return {'sequences': out_sequences, 'sources': out_sources}


@ app.post('/pcr', response_model=create_model(
    'PCRResponse',
    sources=(list[PCRSource], ...),
    sequences=(list[SequenceEntity], ...)
))
async def pcr(source: PCRSource,
              sequences: conlist(SequenceEntity, min_items=1, max_items=1),
              primers: conlist(PrimerModel, min_items=1, max_items=2),):
    dseq = read_dsrecord_from_json(sequences[0])
    primers_pydna = [read_primer_from_json(p) for p in primers]

    # If the footprints are set, the output should be known
    output_is_known = len(source.primer_footprints) > 0 or len(source.primers) > 0 or len(source.fragment_boundaries) > 0

    # Here we set the annealing settings to match the input. If we send the exact annealing,
    # there is no point in sending these settings as well.
    if output_is_known:
        source.primer_annealing_settings = PrimerAnnealingSettings(
            minimum_annealing=min(source.primer_footprints)
        )

    # TODO: return error if the id of the sequence does not correspond.
    # TODO: return error if the ids not present in the list.

    # Error if no pair is generated.
    products, out_sources = get_pcr_products_list(dseq, source, primers_pydna)
    if len(products) == 0:
        raise HTTPException(
            400, 'No pair of annealing primers was found.' + ('' if output_is_known else ' Try changing the annealing settings.')
        )

    out_sequences = [format_sequence_genbank(seq) for seq
                     in products]

    # If the user has provided boundaries, we verify that they are correct.
    if output_is_known:
        for i, out_source in enumerate(out_sources):
            if out_source == source:
                return {'sequences': [out_sequences[i]], 'sources': [out_source]}
        # If we don't find it, there was a mistake
        raise HTTPException(
            400, 'The annealing positions of the primers seem to be wrong.')

    return {'sources': out_sources, 'sequences': out_sequences}
