from fastapi import FastAPI, UploadFile, File, Query, HTTPException
from pydna.dseqrecord import Dseqrecord
from pydantic import conlist, create_model
from pydna.parsers import parse as pydna_parse
from Bio.SeqIO import read as seqio_read
from pydna.genbank import Genbank
from dna_functions import assembly_is_valid, get_assembly_list_from_sticky_ligation_source, get_pcr_products_list, get_restriction_enzyme_products_list, \
    format_sequence_genbank, get_sticky_ligation_products_list, perform_assembly, \
    read_dsrecord_from_json, read_primer_from_json
from pydantic_models import PCRSource, PrimerAnnealingSettings, PrimerModel, SequenceEntity, SequenceFileFormat, \
    GenbankIdSource, RestrictionEnzymeDigestionSource, StickyLigationSource,\
    UploadedFileSource
from fastapi.middleware.cors import CORSMiddleware
from urllib.error import HTTPError, URLError

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


@ app.post('/read_from_file', response_model=create_model(
    'UploadedFileResponse',
    source=(UploadedFileSource, ...),
    sequence=(SequenceEntity, None)
))
async def read_from_file(file: UploadFile = File(...),
                         file_format: SequenceFileFormat = Query(
                             None,
                             description='Format of the sequence file. \
                                Unless specified, it will be guessed\
                                from the extension')
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
    # TODO: return error if the id of the sequence does not correspond
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


@ app.post('/sticky_ligation', response_model=create_model(
    'StickyLigationResponse',
    source=(StickyLigationSource, ...),
    sequence=(SequenceEntity, None)
))
async def sticky_ligation(source: StickyLigationSource,
                          sequences: conlist(SequenceEntity, min_items=1)):
    """Return `output_list` if `source.fragments_inverted` is not set, and `output` otherwise."""
    dseqs = [read_dsrecord_from_json(seq) for seq in sequences]
    if len(source.fragments_inverted) > 0:
        # TODO Error if the list has different order or the ids are wrong.
        # TODO It is problematic that both output_index and fragments_inverted could be set.
        # TODO check input for unique ids
        assembly = get_assembly_list_from_sticky_ligation_source(dseqs, source)
        if not assembly_is_valid(assembly, source.circularised):
            raise HTTPException(
                400, 'Fragments are not compatible for sticky ligation')
        output_sequence = format_sequence_genbank(perform_assembly(assembly, source.circularised))

    else:

        products_dseq, assemblies = get_sticky_ligation_products_list(
            dseqs)
        if len(products_dseq) == 0:
            raise HTTPException(
                400, 'No combination of these fragments is compatible for sticky ligation')

        source.output_list = [format_sequence_genbank(seq) for seq
                              in products_dseq]

        if source.output_index is not None:
            output_sequence = source.output_list[source.output_index]
        else:
            output_sequence = None
    return {'sequence': output_sequence, 'source': source}


@ app.post('/pcr', response_model=create_model(
    'PCRResponse',
    source=(PCRSource, ...),
    sequence=(SequenceEntity, None)
))
async def pcr(source: PCRSource,
              sequences: conlist(SequenceEntity, min_items=1, max_items=1),
              primers: conlist(PrimerModel, min_items=1, max_items=2),):
    dseq = read_dsrecord_from_json(sequences[0])
    primers_pydna = [read_primer_from_json(p) for p in primers]

    # If the primer_pair is set the sender already knows what they want as the output
    output_is_known = source.primer_pair is not None

    # Here we set the annealing settings to match the input. If we send the exact annealing,
    # there is no point in sending these settings as well.
    if output_is_known:
        source.primer_annealing_settings = PrimerAnnealingSettings(
            minimum_annealing=min(source.primer_pair.forward.footprint_length,
                                  source.primer_pair.reverse.footprint_length)
        )

    # TODO: return error if the id of the sequence does not correspond.
    # TODO: return error if the ids not present in the list.

    # Error if no pair is generated.
    products, primer_pairs = get_pcr_products_list(dseq, source, primers_pydna)
    if len(products) == 0:
        raise HTTPException(
            400, 'No pair of annealing primers was found.' + ('' if output_is_known else ' Try changing the annealing settings.')
        )

    if output_is_known:
        dseq_output = next((products[i] for i, pair in enumerate(primer_pairs) if pair == source.primer_pair), None)
        if dseq_output is None:
            raise HTTPException(
                400, 'The annealing positions of the primers seem to be wrong.'
            )
        output_sequence = format_sequence_genbank(dseq_output)
    else:
        source.possible_primer_pairs = primer_pairs
        output_sequence = None
        source.output_list = [format_sequence_genbank(seq) for seq
                              in products]

    return {'sequence': output_sequence, 'source': source}