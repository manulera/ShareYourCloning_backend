import json
from fastapi import FastAPI, UploadFile, File, Query, HTTPException, Request, Body, APIRouter, Form
from typing import Annotated
from fastapi.responses import JSONResponse
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from pydantic import conlist, create_model
from pydna.parsers import parse as pydna_parse
from Bio.SeqIO import read as seqio_read
from pydna.genbank import Genbank
from dna_functions import (
    get_invalid_enzyme_names,
    format_sequence_genbank,
    read_dsrecord_from_json,
    request_from_addgene,
    oligonucleotide_hybridization_overhangs,
)
from pydantic_models import (
    PCRSource,
    PrimerModel,
    SequenceEntity,
    SequenceFileFormat,
    RepositoryIdSource,
    RestrictionEnzymeDigestionSource,
    LigationSource,
    UploadedFileSource,
    HomologousRecombinationSource,
    GibsonAssemblySource,
    RestrictionAndLigationSource,
    AssemblySource,
    GenomeCoordinatesSource,
    ManuallyTypedSource,
    OligoHybridizationSource,
    PolymeraseExtensionSource,
)
from fastapi.middleware.cors import CORSMiddleware
from Bio.Restriction.Restriction import RestrictionBatch
from urllib.error import HTTPError, URLError
from fastapi.responses import HTMLResponse
from Bio.Restriction.Restriction_Dictionary import rest_dict
from assembly2 import (
    Assembly,
    assemble,
    sticky_end_sub_strings,
    PCRAssembly,
    gibson_overlap,
    filter_linear_subassemblies,
    restriction_ligation_overlap,
    SingleFragmentAssembly,
    blunt_overlap,
)
import request_examples
import ncbi_requests
import os
from record_stub_route import RecordStubRoute

record_stubs = os.environ['RECORD_STUBS'] == '1' if 'RECORD_STUBS' in os.environ else False

# Instance of the API object
app = FastAPI()
if record_stubs:
    router = APIRouter(route_class=RecordStubRoute)
else:
    router = APIRouter()


# Allow CORS
# TODO put a wildcard on the shareyourcloning.netlify to
# allow for the draft websites to also work in netlify.
# TODO make this conditional to dev / prod using settings

origins = ['http://localhost:3000', 'https://shareyourcloning.netlify.app']
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=['*'],
    allow_headers=['*'],
)

# TODO limit the maximum size of submitted files


def format_known_assembly_response(
    source: AssemblySource, out_sources: list[AssemblySource], fragments: list[Dseqrecord]
):
    """Common function for assembly sources, when assembly is known"""
    # If a specific assembly is requested
    assembly_plan = source.get_assembly_plan()
    for s in out_sources:
        if s == source:
            return {
                'sequences': [format_sequence_genbank(assemble(fragments, assembly_plan, s.circular))],
                'sources': [s],
            }
    raise HTTPException(400, 'The provided assembly is not valid.')


# Workaround for internal server errors: https://github.com/tiangolo/fastapi/discussions/7847#discussioncomment-5144709
@app.exception_handler(500)
async def custom_http_exception_handler(request: Request, exc: Exception):

    response = JSONResponse(content={'message': 'internal server error'}, status_code=500)

    origin = request.headers.get('origin')

    if origin:
        # Have the middleware do the heavy lifting for us to parse
        # all the config, then update our response headers
        cors = CORSMiddleware(
            app=app, allow_origins=origins, allow_credentials=True, allow_methods=['*'], allow_headers=['*']
        )

        # Logic directly from Starlette's CORSMiddleware:
        # https://github.com/encode/starlette/blob/master/starlette/middleware/cors.py#L152

        response.headers.update(cors.simple_headers)
        has_cookie = 'cookie' in request.headers

        # If request includes any cookie headers, then we must respond
        # with the specific origin instead of '*'.
        if cors.allow_all_origins and has_cookie:
            response.headers['Access-Control-Allow-Origin'] = origin

        # If we only allow specific origins, then we have to mirror back
        # the Origin header in the response.
        elif not cors.allow_all_origins and cors.is_allowed_origin(origin=origin):
            response.headers['Access-Control-Allow-Origin'] = origin
            response.headers.add_vary_header('Origin')

    return response


@router.get('/')
async def greeting(request: Request):
    html_content = f"""
        <html>
            <head>
                <title>Welcome to ShareYourCloning API</title>
            </head>
            <body>
                <h1>Welcome to ShareYourCloning API</h1>
                <p>You can access the endpoints documentation <a href="{request.url._url}docs">here</a></p>
            </body>
        </html>
        """
    return HTMLResponse(content=html_content, status_code=200)


@router.post(
    '/read_from_file',
    response_model=create_model(
        'UploadedFileResponse', sources=(list[UploadedFileSource], ...), sequences=(list[SequenceEntity], ...)
    ),
)
async def read_from_file(
    file: UploadFile = File(...),
    file_format: SequenceFileFormat = Query(
        None,
        description='Format of the sequence file. Unless specified, it will be guessed from the extension',
    ),
    info_str: str = Form(''),
    index_in_file: int = Query(
        None,
        description='The index of the sequence in the file for multi-sequence files',
    ),
):
    """Return a json sequence from a sequence file"""

    if file_format is None:
        extension_dict = {'gbk': 'genbank', 'gb': 'genbank', 'dna': 'snapgene', 'fasta': 'fasta'}
        extension = file.filename.split('.')[-1]
        if extension not in extension_dict:
            raise HTTPException(
                422,
                'We could not guess the format of the file from its extension. Please provide file_format as a query parameter.',
            )

        # We guess the file type from the extension
        file_format = SequenceFileFormat(extension_dict[extension])

    if file_format in ['fasta', 'genbank']:
        # Read the whole file to a string
        file_content = (await file.read()).decode()
        dseqs = pydna_parse(file_content)
        if len(dseqs) == 0:
            raise HTTPException(422, 'Pydna parser reader cannot process this file.')

    elif file_format == 'snapgene':
        try:
            seq = seqio_read(file.file, file_format)
        except ValueError:
            raise HTTPException(422, 'Biopython snapgene reader cannot process this file.')

        iscircular = 'topology' in seq.annotations.keys() and seq.annotations['topology'] == 'circular'
        dseqs = [Dseqrecord(seq, circular=iscircular)]

    # The common part
    parent_source = UploadedFileSource(file_format=file_format, file_name=file.filename)
    out_sources = list()
    info = json.loads(info_str) if info_str else dict()
    for i in range(len(dseqs)):
        new_source = parent_source.model_copy()
        new_source.index_in_file = i
        new_source.info = info
        out_sources.append(new_source)

    out_sequences = [format_sequence_genbank(s) for s in dseqs]

    if index_in_file is not None:
        if index_in_file >= len(out_sources):
            raise HTTPException(404, 'The index_in_file is out of range.')
        return {'sequences': [out_sequences[index_in_file]], 'sources': [out_sources[index_in_file]]}

    return {'sequences': out_sequences, 'sources': out_sources}


# TODO: a bit inconsistent that here you don't put {source: {...}} in the request, but
# directly the object.


@router.post(
    '/repository_id',
    response_model=create_model(
        'RepositoryIdResponse', sources=(list[RepositoryIdSource], ...), sequences=(list[SequenceEntity], ...)
    ),
)
async def get_from_repository_id(source: RepositoryIdSource):

    try:
        if source.repository == 'genbank':
            # This only returns one
            gb = Genbank('example@gmail.com')
            seq = Dseqrecord(gb.nucleotide(source.repository_id))
            dseqs = [seq]
            sources = [source.model_copy()]
        elif source.repository == 'addgene':
            # This may return more than one? Not clear
            dseqs, sources = request_from_addgene(source)
            # Special addgene exception, they may have only partial sequences
            if len(dseqs) == 0:
                raise HTTPException(
                    404,
                    f'The requested plasmid does not exist, or does not have full sequences, see https://www.addgene.org/{source.repository_id}/sequences/',
                )
        else:
            raise HTTPException(404, 'wrong ID')

        output_sequences = [format_sequence_genbank(s) for s in dseqs]

        return {'sequences': output_sequences, 'sources': sources}

    except HTTPError as exception:
        if exception.code == 500:  # pragma: no cover
            raise HTTPException(
                503, f'{source.repository.value} returned: {exception} - {source.repository} might be down'
            )
        elif exception.code == 400 or exception.code == 404:
            raise HTTPException(
                404,
                f'{source.repository.value} returned: {exception} - Likely you inserted a wrong {source.repository} id',
            )
    except URLError as exception:
        raise HTTPException(504, f'Unable to connect to {source.repository}: {exception}')


@router.post(
    '/genome_coordinates',
    response_model=create_model(
        'GenomeRegionResponse', sources=(list[GenomeCoordinatesSource], ...), sequences=(list[SequenceEntity], ...)
    ),
)
async def genome_coordinates(
    source: Annotated[GenomeCoordinatesSource, Body(openapi_examples=request_examples.genome_region_examples)]
):

    # Validate that coordinates make sense
    ncbi_requests.validate_coordinates_pre_request(source.start, source.stop, source.strand)

    # Source includes a locus tag in annotated assembly
    if source.locus_tag is not None:

        if source.assembly_accession is None:
            raise HTTPException(422, 'assembly_accession is required if locus_tag is set')

        annotation = ncbi_requests.get_annotation_from_locus_tag(source.locus_tag, source.assembly_accession)
        gene_range = annotation['genomic_regions'][0]['gene_range']['range'][0]
        gene_strand = 1 if gene_range['orientation'] == 'plus' else -1

        # This field will not be present in all cases, but should be there in reference genomes
        if source.gene_id is not None:
            if 'gene_id' not in annotation:
                raise HTTPException(400, 'gene_id is set, but not found in the annotation')
            if source.gene_id != int(annotation['gene_id']):
                raise HTTPException(400, 'gene_id does not match the locus_tag')
        elif 'gene_id' in annotation:
            source.gene_id = int(annotation['gene_id'])

        # The gene should fall within the range (range might be bigger if bases were requested upstream or downstream)
        if (
            int(gene_range['begin']) < source.start
            or int(gene_range['end']) > source.stop
            or gene_strand != source.strand
        ):
            raise HTTPException(
                400,
                f'wrong coordinates, expected to fall within {source.start}, {source.stop} on strand: {source.strand}',
            )

    elif source.gene_id is not None:
        raise HTTPException(422, 'gene_id is set, but not locus_tag')

    # We get the assembly accession (if it exists), and if the user provided one we validate it
    assembly_accessions = ncbi_requests.get_assembly_accession_from_sequence_accession(source.sequence_accession)

    if source.assembly_accession is not None:
        if source.assembly_accession not in assembly_accessions:
            raise HTTPException(422, 'assembly_accession does not match the one from the sequence_accession')

    # TODO: this could also be not set if there is more than one assembly linked to the sequence
    if len(assembly_accessions):
        source.assembly_accession = assembly_accessions[0]

    seq = ncbi_requests.get_genbank_sequence_subset(
        source.sequence_accession, source.start, source.stop, source.strand
    )

    # NCBI does not complain for coordinates that fall out of the sequence, so we have to check here
    if len(seq) != source.stop - source.start + 1:
        raise HTTPException(400, 'coordinates fall outside the sequence')

    return {'sequences': [format_sequence_genbank(seq)], 'sources': [source.model_copy()]}


@router.post(
    '/manually_typed',
    response_model=create_model(
        'ManuallyTypedResponse', sources=(list[ManuallyTypedSource], ...), sequences=(list[SequenceEntity], ...)
    ),
)
async def manually_typed(source: ManuallyTypedSource):
    """Return the sequence from a manually typed sequence"""
    seq = Dseqrecord(source.user_input, circular=source.circular)
    return {'sequences': [format_sequence_genbank(seq)], 'sources': [source]}


@router.get('/restriction_enzyme_list', response_model=dict[str, list[str]])
async def get_restriction_enzyme_list():
    """Return the dictionary of restriction enzymes"""
    return {'enzyme_names': list(rest_dict.keys())}


@router.post(
    '/restriction',
    response_model=create_model(
        'RestrictionEnzymeDigestionResponse',
        sources=(list[RestrictionEnzymeDigestionSource], ...),
        sequences=(list[SequenceEntity], ...),
    ),
)
async def restriction(
    source: RestrictionEnzymeDigestionSource, sequences: conlist(SequenceEntity, min_length=1, max_length=1)
):

    # TODO: this could be moved to the class
    invalid_enzymes = get_invalid_enzyme_names(source.restriction_enzymes)
    if len(invalid_enzymes):
        raise HTTPException(404, 'These enzymes do not exist: ' + ', '.join(invalid_enzymes))
    enzymes = RestrictionBatch(first=[e for e in source.restriction_enzymes if e is not None])

    seqr = read_dsrecord_from_json(sequences[0])
    # TODO: return error if the id of the sequence does not correspond

    cutsites = seqr.seq.get_cutsites(*enzymes)
    cutsite_pairs = seqr.seq.get_cutsite_pairs(cutsites)
    sources = [RestrictionEnzymeDigestionSource.from_cutsites(*p, source.input, source.id) for p in cutsite_pairs]

    all_enzymes = set(enzyme for s in sources for enzyme in s.restriction_enzymes)
    enzymes_not_cutting = set(source.restriction_enzymes) - set(all_enzymes)
    if len(enzymes_not_cutting):
        raise HTTPException(400, 'These enzymes do not cut: ' + ', '.join(enzymes_not_cutting))

    # If the output is known
    if source.left_edge is not None or source.right_edge is not None:
        # TODO: this could be moved to the class
        if len(source.restriction_enzymes) != 2:
            raise HTTPException(
                400,
                'If either `left_edge` or `right_edge` is provided, the length of `restriction_enzymes` must be 2.',
            )

        for i, s in enumerate(sources):
            if s == source:
                return {'sequences': [format_sequence_genbank(seqr.apply_cut(*cutsite_pairs[i]))], 'sources': [s]}

        raise HTTPException(400, 'Invalid restriction enzyme pair.')

    products = [format_sequence_genbank(seqr.apply_cut(*p)) for p in cutsite_pairs]

    return {'sequences': products, 'sources': sources}


@router.post(
    '/ligation',
    response_model=create_model(
        'LigationResponse', sources=(list[LigationSource], ...), sequences=(list[SequenceEntity], ...)
    ),
)
async def ligation(
    source: LigationSource,
    sequences: conlist(SequenceEntity, min_length=1),
    blunt: bool = Query(False, description='Use blunt ligation instead of sticky ends.'),
    allow_partial_overlap: bool = Query(True, description='Allow for partially overlapping sticky ends.'),
    circular_only: bool = Query(False, description='Only return circular assemblies.'),
):

    # Fragments in the same order as in source.input
    fragments = [
        next((read_dsrecord_from_json(seq) for seq in sequences if seq.id == id), None) for id in source.input
    ]

    # Lambda function for code clarity
    def create_source(a, is_circular):
        return LigationSource.from_assembly(assembly=a, circular=is_circular, id=source.id, input=source.input)

    # If the assembly is known, the blunt parameter is ignored, and we set the algorithm type from the assembly
    # (blunt ligations have features without length)
    if source.assembly is not None:
        asm = source.get_assembly_plan()
        blunt = len(asm[0][2]) == 0

    algo = blunt_overlap if blunt else sticky_end_sub_strings
    out_sources = []
    if len(fragments) > 1:
        asm = Assembly(
            fragments, algorithm=algo, limit=allow_partial_overlap, use_all_fragments=True, use_fragment_order=False
        )
        circular_assemblies = asm.get_circular_assemblies()
        out_sources += [create_source(a, True) for a in circular_assemblies]
        if not circular_only:
            out_sources += [
                create_source(a, False)
                for a in filter_linear_subassemblies(asm.get_linear_assemblies(), circular_assemblies, fragments)
            ]
    else:
        asm = SingleFragmentAssembly(fragments, algorithm=algo, limit=allow_partial_overlap)
        out_sources += [create_source(a, True) for a in asm.get_circular_assemblies()]
        # Not possible to have insertion assemblies in this case

    # If a specific assembly is requested
    if source.assembly is not None:
        return format_known_assembly_response(source, out_sources, fragments)

    if len(out_sources) == 0:
        raise HTTPException(400, 'No ligations were found.')

    out_sequences = [
        format_sequence_genbank(assemble(fragments, s.get_assembly_plan(), s.circular)) for s in out_sources
    ]

    return {'sources': out_sources, 'sequences': out_sequences}


@router.post(
    '/pcr',
    response_model=create_model('PCRResponse', sources=(list[PCRSource], ...), sequences=(list[SequenceEntity], ...)),
)
async def pcr(
    source: PCRSource,
    sequences: conlist(SequenceEntity, min_length=1, max_length=1),
    primers: conlist(PrimerModel, min_length=1, max_length=2),
    minimal_annealing: int = Query(20, description='The minimal annealing length for each primer.'),
    allowed_mismatches: int = Query(0, description='The number of mismatches allowed'),
):

    dseq = read_dsrecord_from_json(sequences[0])
    forward_primer = next((Dseqrecord(Dseq(p.sequence)) for p in primers if p.id == source.forward_primer), None)
    reverse_primer = next((Dseqrecord(Dseq(p.sequence)) for p in primers if p.id == source.reverse_primer), None)
    if forward_primer is None or reverse_primer is None:
        raise HTTPException(404, 'Invalid primer id.')

    # TODO: This may have to be re-written if we allow mismatches
    # If an assembly is provided, we ignore minimal_annealing
    # What happens if annealing is zero? That would mean
    # mismatch in the 3' of the primer, which maybe should
    # not be allowed.
    if source.assembly is not None:
        minimal_annealing = source.minimal_overlap()
        # Only the ones that match are included in the output assembly
        # location, so the submitted assembly should be returned without
        # allowed mistmatches
        # TODO: tests for this
        allowed_mismatches = 0
    fragments = [forward_primer, dseq, reverse_primer]

    asm = PCRAssembly(fragments, limit=minimal_annealing, mismatches=allowed_mismatches)
    try:
        possible_assemblies = asm.get_linear_assemblies()
    except ValueError as e:
        raise HTTPException(400, *e.args)

    out_sources = [
        PCRSource.from_assembly(
            id=source.id,
            input=source.input,
            assembly=a,
            forward_primer=source.forward_primer,
            reverse_primer=source.reverse_primer,
        )
        for a in possible_assemblies
    ]

    # If a specific assembly is requested
    if source.assembly is not None:
        return format_known_assembly_response(source, out_sources, fragments)

    if len(possible_assemblies) == 0:
        raise HTTPException(400, 'No pair of annealing primers was found. Try changing the annealing settings.')

    out_sequences = [
        format_sequence_genbank(assemble(fragments, a, s.circular)) for s, a in zip(out_sources, possible_assemblies)
    ]

    return {'sources': out_sources, 'sequences': out_sequences}


@router.post(
    '/oligonucleotide_hybridization',
    response_model=create_model(
        'OligoHybridizationResponse',
        sources=(list[OligoHybridizationSource], ...),
        sequences=(list[SequenceEntity], ...),
    ),
    openapi_extra={
        'requestBody': {
            'content': {'application/json': {'examples': request_examples.oligonucleotide_hybridization_examples}}
        }
    },
)
async def oligonucleotide_hybridization(
    source: OligoHybridizationSource,
    primers: conlist(PrimerModel, min_length=1, max_length=2),
    minimal_annealing: int = Query(20, description='The minimal annealing length for each primer.'),
):
    watson_seq = next((p.sequence for p in primers if p.id == source.forward_oligo), None)
    crick_seq = next((p.sequence for p in primers if p.id == source.reverse_oligo), None)

    if watson_seq is None or crick_seq is None:
        raise HTTPException(404, 'Invalid oligo id.')

    # The overhang is provided
    if source.overhang_crick_3prime is not None:
        ovhg_watson = len(watson_seq) - len(crick_seq) + source.overhang_crick_3prime
        minimal_annealing = len(watson_seq)
        if source.overhang_crick_3prime < 0:
            minimal_annealing += source.overhang_crick_3prime
        if ovhg_watson > 0:
            minimal_annealing -= ovhg_watson

    possible_overhangs = oligonucleotide_hybridization_overhangs(watson_seq, crick_seq, minimal_annealing)

    if len(possible_overhangs) == 0:
        raise HTTPException(400, 'No pair of annealing oligos was found. Try changing the annealing settings.')

    if source.overhang_crick_3prime is not None:
        if source.overhang_crick_3prime not in possible_overhangs:
            raise HTTPException(400, 'The provided overhang is not compatible with the primers.')

        return {
            'sources': [source],
            'sequences': [
                format_sequence_genbank(Dseqrecord(Dseq(watson_seq, crick_seq, source.overhang_crick_3prime)))
            ],
        }

    out_sources = list()
    out_sequences = list()
    for overhang in possible_overhangs:
        new_source = source.model_copy()
        new_source.overhang_crick_3prime = overhang
        out_sources.append(new_source)
        out_sequences.append(
            format_sequence_genbank(Dseqrecord(Dseq(watson_seq, crick_seq, source.overhang_crick_3prime)))
        )

    return {'sources': out_sources, 'sequences': out_sequences}


@router.post(
    '/polymerase_extension',
    response_model=create_model(
        'PolymeraseExtensionResponse',
        sources=(list[PolymeraseExtensionSource], ...),
        sequences=(list[SequenceEntity], ...),
    ),
)
async def polymerase_extension(
    source: PolymeraseExtensionSource,
    sequences: conlist(SequenceEntity, min_length=1, max_length=1),
):
    """Return the sequence from a polymerase extension reaction"""

    if source.input[0] != sequences[0].id:
        raise HTTPException(404, 'The provided input does not match the sequence id.')

    dseq = read_dsrecord_from_json(sequences[0])

    if dseq.circular:
        raise HTTPException(400, 'The sequence must be linear.')

    if dseq.seq.ovhg == dseq.seq.watson_ovhg() == 0:
        raise HTTPException(400, 'The sequence must have an overhang.')

    out_sequence = Dseqrecord(dseq.seq.fill_in(), features=dseq.features)

    return {'sequences': [format_sequence_genbank(out_sequence)], 'sources': [source]}


@router.post(
    '/homologous_recombination',
    response_model=create_model(
        'HomologousRecombinationResponse',
        sources=(list[HomologousRecombinationSource], ...),
        sequences=(list[SequenceEntity], ...),
    ),
)
async def homologous_recombination(
    source: HomologousRecombinationSource,
    sequences: conlist(SequenceEntity, min_length=2, max_length=2),
    minimal_homology: int = Query(40, description='The minimum homology between the template and the insert.'),
):
    # source.input contains the ids of the sequences in the order template, insert
    template, insert = [
        next((read_dsrecord_from_json(seq) for seq in sequences if seq.id == id), None) for id in source.input
    ]

    if template.circular:
        raise HTTPException(400, 'The template and the insert must be linear.')

    # If an assembly is provided, we ignore minimal_homology
    if source.assembly is not None:
        minimal_homology = source.minimal_overlap()

    asm = Assembly((template, insert), limit=minimal_homology, use_all_fragments=True)

    # The condition is that the first and last fragments are the template
    possible_assemblies = [a for a in asm.get_insertion_assemblies() if a[0][0] == 1]
    out_sources = [
        HomologousRecombinationSource.from_assembly(id=source.id, input=source.input, assembly=a, circular=False)
        for a in possible_assemblies
    ]

    # If a specific assembly is requested
    if source.assembly is not None:
        return format_known_assembly_response(source, out_sources, [template, insert])

    if len(possible_assemblies) == 0:
        raise HTTPException(400, 'No homologous recombination was found.')

    out_sequences = [format_sequence_genbank(assemble([template, insert], a, False)) for a in possible_assemblies]

    return {'sources': out_sources, 'sequences': out_sequences}


@router.post(
    '/gibson_assembly',
    response_model=create_model(
        'GibsonAssemblyResponse', sources=(list[GibsonAssemblySource], ...), sequences=(list[SequenceEntity], ...)
    ),
)
async def gibson_assembly(
    source: GibsonAssemblySource,
    sequences: conlist(SequenceEntity, min_length=1),
    minimal_homology: int = Query(
        40, description='The minimum homology between consecutive fragments in the assembly.'
    ),
    circular_only: bool = Query(False, description='Only return circular assemblies.'),
):

    fragments = [read_dsrecord_from_json(seq) for seq in sequences]

    # Lambda function for code clarity
    def create_source(a, is_circular):
        return GibsonAssemblySource.from_assembly(assembly=a, circular=is_circular, id=source.id, input=source.input)

    out_sources = []
    if len(fragments) > 1:
        asm = Assembly(
            fragments,
            algorithm=gibson_overlap,
            use_fragment_order=False,
            use_all_fragments=True,
            limit=minimal_homology,
        )
        circular_assemblies = asm.get_circular_assemblies()
        out_sources += [create_source(a, True) for a in circular_assemblies]
        if not circular_only:
            out_sources += [
                create_source(a, False)
                for a in filter_linear_subassemblies(asm.get_linear_assemblies(), circular_assemblies, fragments)
            ]
    else:
        asm = SingleFragmentAssembly(fragments, algorithm=gibson_overlap, limit=minimal_homology)
        out_sources += [create_source(a, True) for a in asm.get_circular_assemblies()]
        # Not possible to have insertion assemblies with gibson

    # If a specific assembly is requested
    if source.assembly is not None:
        return format_known_assembly_response(source, out_sources, fragments)

    if len(out_sources) == 0:
        raise HTTPException(400, 'No terminal homology was found.')

    out_sequences = [
        format_sequence_genbank(assemble(fragments, s.get_assembly_plan(), s.circular)) for s in out_sources
    ]

    return {'sources': out_sources, 'sequences': out_sequences}


@router.post(
    '/restriction_and_ligation',
    response_model=create_model(
        'RestrictionAndLigationResponse',
        sources=(list[RestrictionAndLigationSource], ...),
        sequences=(list[SequenceEntity], ...),
    ),
    summary='Restriction and ligation in a single step. Can also be used for Golden Gate assembly.',
)
async def restriction_and_ligation(
    source: RestrictionAndLigationSource,
    sequences: conlist(SequenceEntity, min_length=1),
    allow_partial_overlap: bool = Query(False, description='Allow for partially overlapping sticky ends.'),
    circular_only: bool = Query(False, description='Only return circular assemblies.'),
):

    fragments = [read_dsrecord_from_json(seq) for seq in sequences]
    invalid_enzymes = get_invalid_enzyme_names(source.restriction_enzymes)
    if len(invalid_enzymes):
        raise HTTPException(404, 'These enzymes do not exist: ' + ', '.join(invalid_enzymes))
    enzymes = RestrictionBatch(first=[e for e in source.restriction_enzymes if e is not None])

    # Arguments that are common when instantiating an Assembly pydantic object
    common_args = {'id': source.id, 'input': source.input, 'restriction_enzymes': source.restriction_enzymes}

    # Lambda function for code clarity
    def create_source(a, is_circular):
        return RestrictionAndLigationSource.from_assembly(assembly=a, circular=is_circular, **common_args)

    # Algorithm used by assembly class
    def algo(x, y, _l):
        return restriction_ligation_overlap(x, y, enzymes, allow_partial_overlap)

    out_sources = []
    if len(fragments) > 1:
        asm = Assembly(fragments, algorithm=algo, use_fragment_order=False, use_all_fragments=True)
        circular_assemblies = asm.get_circular_assemblies()
        out_sources += [create_source(a, True) for a in circular_assemblies]
        if not circular_only:
            out_sources += [
                create_source(a, False)
                for a in filter_linear_subassemblies(asm.get_linear_assemblies(), circular_assemblies, fragments)
            ]
    else:
        asm = SingleFragmentAssembly(fragments, algorithm=algo)
        circular_assemblies = asm.get_circular_assemblies()
        out_sources += [create_source(a, True) for a in circular_assemblies]
        if not circular_only:
            out_sources += [create_source(a, False) for a in asm.get_insertion_assemblies()]

    # If a specific assembly is requested
    if source.assembly is not None:
        return format_known_assembly_response(source, out_sources, fragments)

    if len(out_sources) == 0:
        raise HTTPException(400, 'No compatible restriction-ligation was found.')

    out_sequences = [
        format_sequence_genbank(assemble(fragments, s.get_assembly_plan(), s.circular)) for s in out_sources
    ]

    return {'sources': out_sources, 'sequences': out_sequences}


app.include_router(router)
