from fastapi import FastAPI, UploadFile, File, Query, HTTPException, Request
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from pydantic import conlist, create_model
from pydna.parsers import parse as pydna_parse
from Bio.SeqIO import read as seqio_read
from pydna.genbank import Genbank
from dna_functions import get_invalid_enzyme_names, format_sequence_genbank, \
    read_dsrecord_from_json, request_from_addgene
from pydantic_models import PCRSource, PrimerModel, SequenceEntity, SequenceFileFormat, \
    RepositoryIdSource, RestrictionEnzymeDigestionSource, StickyLigationSource, \
    UploadedFileSource, HomologousRecombinationSource, GibsonAssemblySource
from fastapi.middleware.cors import CORSMiddleware
from Bio.Restriction.Restriction import RestrictionBatch
from urllib.error import HTTPError, URLError
from fastapi.responses import HTMLResponse
from Bio.Restriction.Restriction_Dictionary import rest_dict
from assembly2 import Assembly, assemble, sticky_end_sub_strings, is_sublist, assembly2str, PCRAssembly, gibson_overlap
# Instance of the API object
app = FastAPI()

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


@app.get('/')
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
        new_source = parent_source.model_copy()
        new_source.index_in_file = i
        out_sources.append(new_source)

    out_sequences = [format_sequence_genbank(s) for s in dseqs]

    return {'sequences': out_sequences, 'sources': out_sources}

# TODO: a bit inconsistent that here you don't put {source: {...}} in the request, but
# directly the object.


@ app.post('/repository_id', response_model=create_model(
    'RepositoryIdResponse',
    sources=(list[RepositoryIdSource], ...),
    sequences=(list[SequenceEntity], ...)
))
async def get_from_repository_id(source: RepositoryIdSource):

    try:
        if source.repository == 'genbank':
            # This only returns one
            gb = Genbank("example@gmail.com")
            seq = Dseqrecord(gb.nucleotide(source.repository_id))
            dseqs = [seq]
            sources = [source.model_copy()]
        elif source.repository == 'addgene':
            # This may return more than one? Not clear
            dseqs, sources = request_from_addgene(source)
            # Special addgene exception, they may have only partial sequences
            if len(dseqs) == 0:
                raise HTTPException(
                    404, f'The requested plasmid does not exist, or does not have full sequences, see https://www.addgene.org/{source.repository_id}/sequences/')
        else:
            raise HTTPException(400, 'wrong repository name')

        output_sequences = [format_sequence_genbank(s) for s in dseqs]

        return {'sequences': output_sequences, 'sources': sources}

    except HTTPError as exception:
        if exception.code == 500:
            raise HTTPException(
                503, f'{source.repository} returned: {exception} - {source.repository} might be down')
        elif exception.code == 400 or exception.code == 404:
            raise HTTPException(
                404, f'{source.repository} returned: {exception} - Likely you inserted\
                     a wrong {source.repository} id')
    except URLError as exception:
        raise HTTPException(504, f'Unable to connect to {source.repository}: {exception}')


@ app.get('/restriction_enzyme_list', response_model=dict[str, list[str]])
async def get_restriction_enzyme_list():
    """Return the dictionary of restriction enzymes
    """
    return {'enzyme_names': list(rest_dict.keys())}


@ app.post('/restriction', response_model=create_model(
    'RestrictionEnzymeDigestionResponse',
    sources=(list[RestrictionEnzymeDigestionSource], ...),
    sequences=(list[SequenceEntity], ...)
))
async def restriction(source: RestrictionEnzymeDigestionSource,
                      sequences: conlist(SequenceEntity, min_length=1, max_length=1)):

    # TODO: this could be moved to the class
    invalid_enzymes = get_invalid_enzyme_names(source.restriction_enzymes)
    if len(invalid_enzymes):
        raise HTTPException(404, 'These enzymes do not exist: ' + ', '.join(invalid_enzymes))

    seqr = read_dsrecord_from_json(sequences[0])
    # TODO: return error if the id of the sequence does not correspond

    enzymes = RestrictionBatch(first=[e for e in source.restriction_enzymes if e is not None])
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
            raise HTTPException(400, 'If either `left_edge` or `right_edge` is provided, the length of `restriction_enzymes` must be 2.')

        for i, s in enumerate(sources):
            if s == source:
                return {'sequences': [format_sequence_genbank(seqr.apply_cut(*cutsite_pairs[i]))], 'sources': [s]}

        raise HTTPException(400, 'Invalid restriction enzyme pair.')

    products = [format_sequence_genbank(seqr.apply_cut(*p)) for p in cutsite_pairs]

    return {'sequences': products, 'sources': sources}



@ app.post('/sticky_ligation', response_model=create_model(
    'StickyLigationResponse',
    sources=(list[StickyLigationSource], ...),
    sequences=(list[SequenceEntity], ...)
))
async def sticky_ligation(source: StickyLigationSource,
                          sequences: conlist(SequenceEntity, min_length=1),
                          allow_partial_overlap: bool = Query(True, description='Allow for partially overlapping sticky ends.')):

    # Fragments in the same order as in source.input
    fragments = [next((read_dsrecord_from_json(seq) for seq in sequences if seq.id == id), None) for id in source.input]
    if any(f is None for f in fragments):
        raise HTTPException(400, f'Invalid fragment id in input')

    asm = Assembly(fragments, algorithm=sticky_end_sub_strings, limit=allow_partial_overlap, use_all_fragments=True, use_fragment_order=False)
    circular_assemblies = asm.get_circular_assemblies()

    linear_assemblies = asm.get_linear_assemblies()
    # Remove linear assemblies which are sub-assemblies of circular assemblies
    linear_assemblies = [a for a in linear_assemblies if not any(is_sublist(a, c, True) for c in circular_assemblies)]
    possible_assemblies = circular_assemblies + linear_assemblies

    out_sources = [StickyLigationSource.from_assembly(id= source.id, input=source.input, assembly=a, circular=(a[0][0] == a[-1][1])) for a in possible_assemblies]

    # If a specific assembly is requested
    if source.assembly is not None:
        for i, s in enumerate(out_sources):
            if s == source:
                # TODO: remove enumeration by better parsing
                return {'sequences': [format_sequence_genbank(assemble(fragments, possible_assemblies[i], s.circular))], 'sources': [s]}
        raise HTTPException(400, 'The provided assembly is not valid.')

    if len(possible_assemblies) == 0:
        raise HTTPException(400, 'No ligations were found.')

    out_sequences = [format_sequence_genbank(assemble(fragments, a, s.circular)) for s, a in zip(out_sources, possible_assemblies)]

    return {'sources': out_sources, 'sequences': out_sequences}


@ app.post('/pcr', response_model=create_model(
    'PCRResponse',
    sources=(list[PCRSource], ...),
    sequences=(list[SequenceEntity], ...)
))
async def pcr(source: PCRSource,
              sequences: conlist(SequenceEntity, min_length=1, max_length=1),
              primers: conlist(PrimerModel, min_length=1, max_length=2),
              minimal_annealing: int = Query(20, description='The minimal annealing length for each primer.')):

    dseq = read_dsrecord_from_json(sequences[0])
    forward_primer = next((Dseqrecord(Dseq(p.sequence)) for p in primers if p.id == source.forward_primer), None)
    reverse_primer = next((Dseqrecord(Dseq(p.sequence)) for p in primers if p.id == source.reverse_primer), None)
    if forward_primer is None or reverse_primer is None:
        raise HTTPException(404, 'Invalid primer id.')

    # TODO: This will have to be re-written if we allow mismatches
    # If an assembly is provided, we ignore minimal_annealing
    if source.assembly is not None:
        minimal_annealing = source.minimal_overlap()
    fragments = [forward_primer, dseq, reverse_primer]

    asm = PCRAssembly(fragments, limit=minimal_annealing)
    try:
        print(asm.G.edges)
        possible_assemblies = asm.get_linear_assemblies()
    except ValueError as e:
        raise HTTPException(400, *e.args)

    out_sources = [PCRSource.from_assembly(id= source.id, input=source.input, assembly=a, forward_primer=source.forward_primer, reverse_primer=source.reverse_primer) for a in possible_assemblies]

    # If a specific assembly is requested
    if source.assembly is not None:
        for i, s in enumerate(out_sources):
            if s == source:
                # TODO: remove enumeration by better parsing
                return {'sequences': [format_sequence_genbank(assemble(fragments, possible_assemblies[i], s.circular))], 'sources': [s]}
        raise HTTPException(400, 'The provided assembly is not valid.')

    if len(possible_assemblies) == 0:
        raise HTTPException(400, 'No pair of annealing primers was found. Try changing the annealing settings.')

    out_sequences = [format_sequence_genbank(assemble(fragments, a, s.circular)) for s, a in zip(out_sources, possible_assemblies)]

    return {'sources': out_sources, 'sequences': out_sequences}


@ app.post('/homologous_recombination', response_model=create_model(
    'HomologousRecombinationResponse',
    sources=(list[HomologousRecombinationSource], ...),
    sequences=(list[SequenceEntity], ...)
))
async def homologous_recombination(
    source: HomologousRecombinationSource,
    sequences: conlist(SequenceEntity, min_length=2, max_length=2),
    minimal_homology: int = Query(40, description='The minimum homology between the template and the insert.')
):
    template = read_dsrecord_from_json(sequences[0])
    insert = read_dsrecord_from_json(sequences[1])

    # If an assembly is provided, we ignore minimal_homology
    if source.assembly is not None:
        minimal_homology = source.minimal_overlap()

    asm = Assembly((template, insert, template), limit=minimal_homology, use_all_fragments=True)

    # The condition is that the first and last fragments are the template
    possible_assemblies = [a for a in asm.get_linear_assemblies() if a[0][0] == 1 and a[-1][1] == 3]

    # Replace the index of last fragment (3) by 1, since it is repeated
    possible_assemblies = [(a[0], (2, 1, a[1][2], a[1][3]), ) for a in possible_assemblies]
    out_sources = [HomologousRecombinationSource.from_assembly(id= source.id, input=source.input, assembly=a, circular=False) for a in possible_assemblies]


    # If a specific assembly is requested
    if source.assembly is not None:
        for i, s in enumerate(out_sources):
            if s == source:
                return {'sequences': [format_sequence_genbank(assemble([template, insert], possible_assemblies[i], False))], 'sources': [s]}
        raise HTTPException(400, 'The provided assembly is not valid.')

    if len(possible_assemblies) == 0:
        raise HTTPException(400, 'No homologous recombination was found.')

    out_sequences = [format_sequence_genbank(assemble([template, insert], a, False)) for a in possible_assemblies]

    return {'sources': out_sources, 'sequences': out_sequences}

@ app.post('/gibson_assembly', response_model=create_model(
    'GibsonAssemblyResponse',
    sources=(list[GibsonAssemblySource], ...),
    sequences=(list[SequenceEntity], ...)
))
async def gibson_assembly(source: GibsonAssemblySource,
                          sequences: conlist(SequenceEntity, min_length=1),
                          minimal_homology: int = Query(40, description='The minimum homology between consecutive fragments in the assembly.'),):

    fragments = [read_dsrecord_from_json(seq) for seq in sequences]

    asm = Assembly(fragments, algorithm=gibson_overlap, limit=minimal_homology, use_all_fragments=True, use_fragment_order=False)
    circular_assemblies = asm.get_circular_assemblies()

    linear_assemblies = asm.get_linear_assemblies()
    # Remove linear assemblies which are sub-assemblies of circular assemblies
    linear_assemblies = [a for a in linear_assemblies if not any(is_sublist(a, c, True) for c in circular_assemblies)]
    possible_assemblies = circular_assemblies + linear_assemblies

    out_sources = [GibsonAssemblySource.from_assembly(id= source.id, input=source.input, assembly=a, circular=(a[0][0] == a[-1][1])) for a in possible_assemblies]

    # If a specific assembly is requested
    if source.assembly is not None:
        for i, s in enumerate(out_sources):
            if s == source:
                return {'sequences': [format_sequence_genbank(assemble(fragments, possible_assemblies[i], s.circular))], 'sources': [s]}
        raise HTTPException(400, 'The provided assembly is not valid.')

    if len(possible_assemblies) == 0:
        raise HTTPException(400, 'No terminal homology was found.')

    out_sequences = [format_sequence_genbank(assemble(fragments, a, s.circular)) for s, a in zip(out_sources, possible_assemblies)]

    return {'sources': out_sources, 'sequences': out_sequences}

