from dna_functions import format_sequence_genbank, read_dsrecord_from_json
from main import app
from fastapi.testclient import TestClient
from pydna.parsers import parse as pydna_parse
from Bio.Restriction.Restriction import CommOnly
from Bio.SeqFeature import SimpleLocation
from pydantic_models import (
    RepositoryIdSource,
    PCRSource,
    PrimerModel,
    RestrictionEnzymeDigestionSource,
    TextFileSequence,
    LigationSource,
    UploadedFileSource,
    HomologousRecombinationSource,
    GibsonAssemblySource,
    RestrictionAndLigationSource,
    ManuallyTypedSource,
    GenomeCoordinatesSource,
    OligoHybridizationSource,
    PolymeraseExtensionSource,
    CRISPRSource,
    AddGeneIdSource,
    RestrictionSequenceCut,
    BaseCloningStrategy,
    SimpleSequenceLocation as PydanticSimpleLocation,
    BenchlingUrlSource,
    EuroscarfSource,
    SnapGenePlasmidSource,
)
from pydna.dseqrecord import Dseqrecord
import unittest
from pydna.dseq import Dseq
import request_examples
import copy
import json
import tempfile
import pytest
from Bio.Seq import reverse_complement
import os


# Custom decorator to run code before and after a specific test, with protection against interruptions
def run_before_after(before_func, after_func):
    def decorator(func):
        def wrapper(*args, **kwargs):
            before_func()  # Run the before function
            try:
                return func(*args, **kwargs)  # Run the actual test
            finally:
                after_func()  # Always run the after function, even if there's an exception

        return wrapper

    return decorator


client = TestClient(app)


class VersionTest(unittest.TestCase):
    def test_version_empty(self):
        response = client.get('/version')
        self.assertEqual(response.status_code, 200)
        resp = response.json()
        self.assertIsNone(resp['version'])
        self.assertIsNone(resp['commit_sha'])

    def create_version_files():
        with open('./version.txt', 'w') as f:
            f.write('1.2.3')
        with open('./commit_sha.txt', 'w') as f:
            f.write('1234567890')

    def delete_version_files():
        os.remove('./version.txt')
        os.remove('./commit_sha.txt')

    @run_before_after(create_version_files, delete_version_files)
    def test_version_file(self):
        response = client.get('/version')
        self.assertEqual(response.status_code, 200)
        resp = response.json()
        self.assertEqual(resp['version'], '1.2.3')
        self.assertEqual(resp['commit_sha'], '1234567890')


class ReadFileTest(unittest.TestCase):
    def test_read_files(self):
        """Test that uploading files with single and multiple sequences works."""

        examples = [
            {'file': './examples/sequences/pFA6a-hphMX6.gb', 'format': 'genbank', 'nb_sequences': 1},
            {'file': './examples/sequences/dummy_EcoRI.fasta', 'format': 'fasta', 'nb_sequences': 1},
            {'file': './examples/sequences/dummy_multi_fasta.fasta', 'format': 'fasta', 'nb_sequences': 2},
            {
                'file': './examples/sequences/addgene-plasmid-39296-sequence-49545.dna',
                'format': 'snapgene',
                'nb_sequences': 1,
            },
            {'file': './examples/sequences/ase1.embl', 'format': 'embl', 'nb_sequences': 1},
            # Ape files as of 2024-10-30 did not have a properly formatted LOCUS line
            {
                'file': './test_files/example.ape',
                'format': 'genbank',
                'nb_sequences': 1,
                'warning': True,
                'circular': False,
            },
            # Euroscarf files as of 2024-10-30 did not have a properly formatted LOCUS line
            {
                'file': './test_files/pKT128_euroscarf.gb',
                'format': 'genbank',
                'nb_sequences': 1,
                'warning': True,
                'circular': True,
            },
        ]

        for example in examples:
            with open(example['file'], 'rb') as f:
                response = client.post('/read_from_file', files={'file': f})
            self.assertEqual(response.status_code, 200)
            payload = response.json()

            resulting_sequences = [
                read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
            ]
            sources = [UploadedFileSource.model_validate(s) for s in payload['sources']]

            self.assertEqual(len(sources), example['nb_sequences'])
            self.assertEqual(len(resulting_sequences), example['nb_sequences'])
            for seq in resulting_sequences:
                self.assertGreater(len(seq), 3)
            for source in sources:
                self.assertEqual(source.sequence_file_format, example['format'])

            if 'warning' in example and example['warning']:
                self.assertIn('x-warning', response.headers)
            else:
                self.assertNotIn('x-warning', response.headers)
            if 'circular' in example:
                for seq in resulting_sequences:
                    self.assertEqual(seq.circular, example['circular'])

        # Test naming
        with open(example['file'], 'rb') as f:
            response = client.post('/read_from_file?output_name=hello', files={'file': f})
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        dseqr = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertEqual(dseqr.name, 'hello')

    def test_errors_read_files(self):
        # Create a temp empty file
        with tempfile.NamedTemporaryFile() as temp_empty_file:

            examples = [
                {
                    'file': './test_endpoints.py',
                    'format': None,
                    'error_message': 'could not guess',
                },
                {
                    'file': './test_endpoints.py',
                    'format': 'snapgene',
                    'error_message': 'Biopython cannot process',
                },
                {
                    'file': './test_endpoints.py',
                    'format': 'genbank',
                    'error_message': 'Biopython cannot process',
                },
                {
                    'file': './test_endpoints.py',
                    'format': 'embl',
                    'error_message': 'Biopython cannot process',
                },
                {
                    'file': './test_endpoints.py',
                    'format': 'fasta',
                    'error_message': 'Biopython cannot process',
                },
                {
                    'file': temp_empty_file.name,
                    'format': 'genbank',
                    'error_message': 'Biopython cannot process',
                },
                {
                    'file': temp_empty_file.name,
                    'format': 'embl',
                    'error_message': 'Biopython cannot process',
                },
                {
                    'file': temp_empty_file.name,
                    'format': 'fasta',
                    'error_message': 'Biopython cannot process',
                },
            ]
            for example in examples:
                if example['format'] is not None:
                    response = client.post(
                        f'/read_from_file?sequence_file_format={example["format"]}',
                        files={'file': temp_empty_file},
                    )
                else:
                    response = client.post('/read_from_file', files={'file': temp_empty_file})

                self.assertNotEqual(response.status_code, 200)
                self.assertIn(example['error_message'], response.json()['detail'])

    def test_file_index_known(self):
        """Test that if the index in file is specified it works."""

        with open('./examples/sequences/dummy_multi_fasta.fasta', 'rb') as f:
            response = client.post('/read_from_file?index_in_file=1', files={'file': f})

        self.assertEqual(response.status_code, 200)
        payload = response.json()

        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        sources = [UploadedFileSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(sources), 1)
        self.assertEqual(len(resulting_sequences), 1)

    def test_circularize_fasta_sequence(self):
        """Test that the circularize parameter works when reading from file"""
        with open('./examples/sequences/dummy_EcoRI.fasta', 'rb') as f:
            response = client.post('/read_from_file?circularize=True', files={'file': f})

        self.assertEqual(response.status_code, 200)
        payload = response.json()

        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        sources = [UploadedFileSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(sources), 1)
        self.assertEqual(len(resulting_sequences), 1)
        self.assertTrue(sources[0].circularize)
        self.assertTrue(resulting_sequences[0].circular)

        # If the file format is different it raises an error
        with open('./examples/sequences/dummy_EcoRI.fasta', 'rb') as f:
            response = client.post('/read_from_file?circularize=True&sequence_file_format=genbank', files={'file': f})

        self.assertEqual(response.status_code, 400)

        # Also when guessing from a gb file
        with open('./examples/sequences/pFA6a-hphMX6.gb', 'rb') as f:
            response = client.post('/read_from_file?circularize=True', files={'file': f})

        self.assertEqual(response.status_code, 400)


class GenBankTest(unittest.TestCase):

    # TODO these tests will not work off-line, so the case where connection cannot be established should be handled in some way
    def test_request_gene(self):
        """Test whether the gene is requested from GenBank"""
        source = RepositoryIdSource(
            id=1,
            repository_name='genbank',
            repository_id='NM_001018957.2',
        )
        response = client.post('/repository_id/genbank', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        sequence = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertIn('Ase1', sequence.description)

    def test_request_wrong_id(self):
        """Test a wrong Genbank id"""
        source = RepositoryIdSource(
            id=1,
            repository_name='genbank',
            repository_id='wrong_id',
        )
        response = client.post('/repository_id/genbank', json=source.model_dump())
        self.assertEqual(response.status_code, 404)

    def test_redirect(self):
        """The repository_id endpoint should redirect based on repository_name value"""
        source = RepositoryIdSource(
            id=1,
            repository_name='genbank',
            repository_id='NM_001018957.2',
        )
        response = client.post('/repository_id', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        sequence: Dseqrecord = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertIn('Ase1', sequence.description)

    def test_rename(self):
        """If passing output_name, it renames the output"""
        source = RepositoryIdSource(
            id=1,
            repository_name='genbank',
            repository_id='NM_001018957.2',
            output_name='hello',
        )
        response = client.post('/repository_id/genbank', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        sequence = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertEqual(sequence.name, 'hello')

    def test_long_sequence(self):
        """Test that a long sequence raises an error"""
        source = RepositoryIdSource(
            id=1,
            repository_name='genbank',
            repository_id='CU329670.1',
        )
        response = client.post('/repository_id/genbank', json=source.model_dump())
        self.assertEqual(response.status_code, 400)
        self.assertIn('sequence is too long', response.json()['detail'])


class AddGeneTest(unittest.TestCase):
    def test_request_plasmid(self):
        """Test whether the gene is requested from AddGene and returns the right info"""
        examples = [
            {
                'id': '39282',
                'url': 'https://media.addgene.org/snapgene-media/v2.0.0/sequences/240599/4936a6ae-6b4d-4d24-b7ac-2339fad5755d/addgene-plasmid-39282-sequence-240599.gbk',
                'type': 'addgene-full',
                'name': 'pFA6a-kanMX6-P81nmt1',
            },
            {
                'id': '39289',
                'url': 'https://media.addgene.org/snapgene-media/v2.0.0/sequences/49640/4fa9f18b-d5ca-4ac6-a50c-76a6cd15cbab/addgene-plasmid-39289-sequence-49640.gbk',
                'type': 'depositor-full',
                'name': 'pFA6a-kanMX6-P3nmt1-GFP',
            },
        ]
        for example in examples:
            source = AddGeneIdSource(
                id=1,
                repository_name='addgene',
                repository_id=example['id'],
            )

            response = client.post('/repository_id/addgene', json=source.model_dump())
            self.assertEqual(response.status_code, 200)
            payload = response.json()
            resulting_sequences = [
                read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
            ]
            sources = [AddGeneIdSource.model_validate(s) for s in payload['sources']]

            self.assertEqual(len(resulting_sequences), 1)
            self.assertEqual(len(sources), 1)
            self.assertEqual(sources[0].addgene_sequence_type, example['type'])
            self.assertEqual(sources[0].sequence_file_url, example['url'])
            self.assertEqual(resulting_sequences[0].name, example['name'])

            # We get the same response when making the response with the url
            response2 = client.post('/repository_id/addgene', json=payload['sources'][0])
            self.assertEqual(response.json(), response2.json())

    @pytest.mark.xfail(reason='This file was removed from AddGene, not sure what the ideal behavior should be')
    def test_old_url(self):
        """Works for an AddGene url that has now been replaced by a newer one"""
        source = AddGeneIdSource(
            id=1,
            repository_name='addgene',
            repository_id='65109',
            addgene_sequence_type='addgene-full',
            sequence_file_url='https://media.addgene.org/snapgene-media/v1.7.9-0-g88a3305/sequences/110162/c1c98803-c8ba-44a6-95b8-d6a94097e36f/addgene-plasmid-65109-sequence-110162.gbk',
        )
        response = client.post('/repository_id/addgene', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        sources = [AddGeneIdSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)
        self.assertEqual(sources[0].addgene_sequence_type, 'addgene-full')
        self.assertEqual(
            sources[0].sequence_file_url,
            'https://media.addgene.org/snapgene-media/v1.7.9-0-g88a3305/sequences/110162/c1c98803-c8ba-44a6-95b8-d6a94097e36f/addgene-plasmid-65109-sequence-110162.gbk',
        )
        self.assertEqual(resulting_sequences[0].name, 'pYTK002')

    def test_missing_sequences(self):
        # Non-existing id
        source = AddGeneIdSource(
            id=1,
            repository_name='addgene',
            repository_id='DUMMYTEST',
        )

        response = client.post('/repository_id/addgene', json=source.model_dump())
        self.assertEqual(response.status_code, 404)
        self.assertIn('wrong addgene id', response.json()['detail'])

        # Id that has no full-sequences
        source = AddGeneIdSource(
            id=1,
            repository_name='addgene',
            repository_id='39291',
        )
        response = client.post('/repository_id/addgene', json=source.model_dump())
        self.assertEqual(response.status_code, 404)
        self.assertIn('The requested plasmid does not have full sequences', response.json()['detail'])

        # url does not exist
        source = AddGeneIdSource(
            id=1,
            repository_name='addgene',
            repository_id='39282',
            sequence_file_url='https://media.addgene.org/snapgene-media/wrongggggggg.gbk',
        )
        response = client.post('/repository_id/addgene', json=source.model_dump())
        self.assertEqual(response.status_code, 404)

        # TODO url exists but does not match id, or does not match

    def test_redirect(self):
        """Test repository_id endpoint should redirect based on repository_name value"""
        source = AddGeneIdSource(
            id=1,
            repository_name='addgene',
            repository_id='39282',
        )
        response = client.post('/repository_id', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        sequence: Dseqrecord = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertIn('synthetic circular DNA', sequence.description)


class LigationTest(unittest.TestCase):
    def test_ligation(self):
        """Test whether assembly is executed when the order is provided"""
        # Load dummy sequence

        initial_sequence = pydna_parse('examples/sequences/dummy_EcoRI.fasta')[0]

        # Restriction cut
        enzyme = CommOnly.format('EcoRI')
        output_list: list[Dseqrecord] = initial_sequence.cut([enzyme])

        # Convert to json to use as input
        model_seqs = [format_sequence_genbank(seq) for seq in output_list]
        model_seqs[0].id = 1
        model_seqs[1].id = 2
        json_seqs = [seq.model_dump() for seq in model_seqs]

        # Assign ids to define deterministic assembly
        source = LigationSource.from_assembly(
            id=0,
            assembly=[(1, 2, SimpleLocation(7, 11), SimpleLocation(0, 4))],
            circular=False,
            fragments=model_seqs,
        )
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post('/ligation', json=data)
        payload = response.json()

        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        sources = [LigationSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)

        # Check that the assembly is correct
        self.assertEqual(resulting_sequences[0].seq, initial_sequence.seq)
        self.assertEqual(sources[0], source)

        # Check that the inverse assembly will not pass
        source = LigationSource.from_assembly(
            id=0,
            assembly=[(1, 2, SimpleLocation(8, 11), SimpleLocation(1, 4))],
            circular=False,
            fragments=model_seqs,
        )
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post('/ligation', json=data)
        self.assertEqual(response.status_code, 400)
        data = response.json()
        self.assertEqual(data['detail'], 'The provided assembly is not valid.')

        # Check that the circular assembly does not pass either
        source = LigationSource.from_assembly(
            id=0,
            assembly=[(1, 2, SimpleLocation(7, 11), SimpleLocation(0, 4))],
            circular=True,
            fragments=model_seqs,
        )
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post('/ligation', json=data)
        self.assertEqual(response.status_code, 400)
        data = response.json()
        self.assertEqual(data['detail'], 'The provided assembly is not valid.')

    def test_linear_assembly_no_order(self):
        """Test that when order is not provided, no duplicate sequences are returned as options."""

        # Load Ase1 sequence
        initial_sequence = pydna_parse('examples/sequences/dummy_EcoRI.fasta')[0]

        # Restriction cut
        enzyme = CommOnly.format('EcoRI')
        output_list: list[Dseqrecord] = initial_sequence.cut([enzyme])

        # Convert to json to use as input
        json_seqs = [format_sequence_genbank(seq) for seq in output_list]
        json_seqs[0].id = 1
        json_seqs[1].id = 2
        json_seqs = [seq.model_dump() for seq in json_seqs]

        # We don't set the fragments_inverted, so we will get all possibilities (in this case only one)
        source = LigationSource(id=0)
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post('/ligation', json=data)
        payload = response.json()

        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        sources = [LigationSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)

        # Check that the assembly is correct
        self.assertEqual(resulting_sequences[0].seq, initial_sequence.seq)

    def test_circular_assembly_no_order(self):
        initial_sequence = pydna_parse('examples/sequences/pFA6a-hphMX6.gb')[0]

        # Restriction cut
        output_list: list[Dseqrecord] = initial_sequence.cut([CommOnly.format('AscI'), CommOnly.format('SacI')])

        # Convert to json to use as input
        json_seqs = [format_sequence_genbank(seq) for seq in output_list]
        json_seqs[0].id = 1
        json_seqs[1].id = 2
        json_seqs = [seq.model_dump() for seq in json_seqs]

        source = LigationSource(id=0)
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post('/ligation', json=data)
        payload = response.json()
        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        sources = [LigationSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)

        # Check that the assembly is correct, the sequences cannot be compared with the equal operator since their origin is different
        self.assertEqual(resulting_sequences[0].seq.seguid(), initial_sequence.seq.seguid())

    def test_circularization(self):
        enzyme = CommOnly.format('EcoRI')
        fragment = Dseqrecord('AGAATTC', circular=True).cut(enzyme)[0]
        json_seqs = [format_sequence_genbank(fragment)]
        json_seqs[0].id = 1
        json_seqs = [seq.model_dump() for seq in json_seqs]
        source = LigationSource(id=0)
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post('/ligation', json=data)
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        self.assertEqual(resulting_sequences[0].seq, fragment.looped().seq)

    def test_blunt_ligation(self):
        seqs = [Dseqrecord('ATCC', circular=False), Dseqrecord('TAAT', circular=False)]
        json_seqs = [format_sequence_genbank(seq) for seq in seqs]
        json_seqs[0].id = 1
        json_seqs[1].id = 2
        json_seqs = [seq.model_dump() for seq in json_seqs]
        source = LigationSource(id=0)
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post('/ligation', json=data, params={'blunt': True})
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        self.assertEqual(len(resulting_sequences), 2)

        # We submit one of the resulting sources, to check that it does the
        # blunt ligation without adding the request parameter
        source = payload['sources'][0]
        response = client.post('/ligation', json={'source': source, 'sequences': json_seqs})
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        resulting_sequences2 = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        self.assertEqual(len(resulting_sequences2), 1)
        self.assertEqual(resulting_sequences2[0].seq, resulting_sequences[0].seq)

    def test_blunt_circularization(self):
        seqs = [Dseqrecord('ATCC', circular=False)]
        json_seqs = [format_sequence_genbank(seq) for seq in seqs]
        json_seqs[0].id = 1
        json_seqs = [seq.model_dump() for seq in json_seqs]
        source = LigationSource(id=0)
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post('/ligation', json=data, params={'blunt': True})
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(resulting_sequences[0].seq, seqs[0].looped().seq)

    def test_mixed_ligation(self):
        seqs = [
            Dseqrecord(Dseq.from_full_sequence_and_overhangs('ACCGTA', -3, 0)),
            Dseqrecord(Dseq.from_full_sequence_and_overhangs('AAGACC', 0, -3)),
        ]
        json_seqs = [format_sequence_genbank(seq) for seq in seqs]
        json_seqs[0].id = 1
        json_seqs[1].id = 2
        json_seqs = [seq.model_dump() for seq in json_seqs]
        source = LigationSource(id=0)
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post('/ligation', json=data, params={'blunt': True})
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(resulting_sequences[0].seguid(), (seqs[1] + seqs[0]).looped().seguid())

    def test_too_many_assemblies(self):
        # Too many paths
        seqs = [Dseqrecord('ATCC', circular=False)] * 20
        json_seqs = [format_sequence_genbank(seq) for seq in seqs]
        json_seqs = [seq.model_dump() for seq in json_seqs]
        for i in range(len(json_seqs)):
            json_seqs[i]['id'] = i
        source = LigationSource(id=0)
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post('/ligation', json=data, params={'blunt': True})
        self.assertEqual(response.status_code, 400)
        self.assertTrue('Too many possible paths' in response.json()['detail'])


class RestrictionTest(unittest.TestCase):
    def test_enzyme_doesnt_exist(self):
        dseq = Dseqrecord('AAAAAAGAATTCTTTTTT', circular=False)
        json_seq = format_sequence_genbank(dseq)
        json_seq.id = 1

        # One enzyme
        source = RestrictionEnzymeDigestionSource(id=0)
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data, params={'restriction_enzymes': ['helloworld']})
        self.assertEqual(response.status_code, 404)
        self.assertTrue(response.json()['detail'] == 'These enzymes do not exist: helloworld')

        # More than one
        source = RestrictionEnzymeDigestionSource(id=0)
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data, params={'restriction_enzymes': ['helloworld', 'byebye']})
        self.assertEqual(response.status_code, 404)
        self.assertIn('byebye', response.json()['detail'])
        self.assertIn('helloworld', response.json()['detail'])

    def test_enzymes_dont_cut(self):
        dseq = Dseqrecord('AAAAAAGAATTCTTTTTT', circular=False)
        json_seq = format_sequence_genbank(dseq)
        json_seq.id = 1

        # See if we get the right responses
        source = RestrictionEnzymeDigestionSource(id=0)
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data, params={'restriction_enzymes': ['FbaI']})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json()['detail'], 'These enzymes do not cut: FbaI')

        # Returns when one does not cut
        source = RestrictionEnzymeDigestionSource(id=0)
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data, params={'restriction_enzymes': ['EcoRI', 'BamHI']})
        self.assertEqual(response.status_code, 400)
        self.assertTrue(response.json()['detail'] == 'These enzymes do not cut: BamHI')

        # Returns all that don't cut
        source = RestrictionEnzymeDigestionSource(id=0)
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data, params={'restriction_enzymes': ['BamHI', 'FbaI']})
        self.assertEqual(response.status_code, 400)
        self.assertIn('BamHI', response.json()['detail'])
        self.assertIn('FbaI', response.json()['detail'])

    def test_linear_single_restriction(self):

        dseq = Dseqrecord('AAAAAAGAATTCTTTTTT', circular=False)
        json_seq = format_sequence_genbank(dseq)
        json_seq.id = 1

        # See if we get the right responses
        source = RestrictionEnzymeDigestionSource(id=0)
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data, params={'restriction_enzymes': ['EcoRI']})
        payload = response.json()

        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        sources = [RestrictionEnzymeDigestionSource.model_validate(s) for s in payload['sources']]

        # The right sequences are returned in the right order
        self.assertEqual(resulting_sequences[0].seq.watson.upper(), 'AAAAAAG')
        self.assertEqual(resulting_sequences[1].seq.watson.upper(), 'AATTCTTTTTT')

        # The edges of the fragments are correct:
        self.assertEqual(sources[0].left_edge, None)
        self.assertEqual(sources[0].right_edge.to_cutsite_tuple()[0], (7, -4))

        self.assertEqual(sources[1].left_edge.to_cutsite_tuple()[0], (7, -4))
        self.assertEqual(sources[1].right_edge, None)

        # The enzyme names are correctly returned:
        self.assertEqual(sources[0].get_enzymes(), ['EcoRI'])
        self.assertEqual(sources[1].get_enzymes(), ['EcoRI'])

        # Now we specify the output
        source = RestrictionEnzymeDigestionSource(
            id=0, left_edge=None, right_edge=RestrictionSequenceCut.from_cutsite_tuple(((7, -4), 'EcoRI'))
        )
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data)
        payload = response.json()

        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        sources = [RestrictionEnzymeDigestionSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)
        self.assertEqual(sources[0].left_edge, None)
        self.assertEqual(sources[0].right_edge.to_cutsite_tuple()[0], (7, -4))
        self.assertEqual(str(sources[0].right_edge.to_cutsite_tuple()[1]), 'EcoRI')
        self.assertEqual(sources[0].get_enzymes(), ['EcoRI'])

    def test_circular_single_restriction(self):

        dseq = Dseqrecord('AAAAAAGAATTCTTTTTT', circular=True)
        json_seq = format_sequence_genbank(dseq)
        json_seq.id = 1

        # See if we get the right responses
        source = RestrictionEnzymeDigestionSource(id=0)
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data, params={'restriction_enzymes': ['EcoRI']})
        payload = response.json()

        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        sources = [RestrictionEnzymeDigestionSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)

        # The right sequences are returned in the right order
        self.assertEqual(resulting_sequences[0].seq.watson.upper(), 'AATTCTTTTTTAAAAAAG')

        # The edges of the fragments are correct:
        self.assertEqual(sources[0].left_edge.to_cutsite_tuple()[0], (7, -4))
        self.assertEqual(sources[0].right_edge.to_cutsite_tuple()[0], (7, -4))

        # The enzyme names are correctly returned:
        self.assertEqual(sources[0].get_enzymes(), ['EcoRI'])

        # When the cutting site spans the origin
        sequences = ['AATTCTTTTTTG', 'ATTCTTTTTTGA']
        cut_positions = [(0, -4), (11, -4)]

        for s, pos in zip(sequences, cut_positions):
            dseq = Dseqrecord(s, circular=True)
            json_seq = format_sequence_genbank(dseq)
            json_seq.id = 1

            # See if we get the right responses
            source = RestrictionEnzymeDigestionSource(
                id=0,
            )
            data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
            response = client.post('/restriction', json=data, params={'restriction_enzymes': ['EcoRI']})
            payload = response.json()

            resulting_sequences = [
                read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
            ]
            sources = [RestrictionEnzymeDigestionSource.model_validate(s) for s in payload['sources']]

            self.assertEqual(len(resulting_sequences), 1)
            self.assertEqual(len(sources), 1)

            # The right sequences are returned in the right order
            self.assertEqual(resulting_sequences[0].seq.watson.upper(), 'AATTCTTTTTTG')

            # The edges of the fragments are correct:
            self.assertEqual(sources[0].left_edge.to_cutsite_tuple()[0], pos)
            self.assertEqual(sources[0].right_edge.to_cutsite_tuple()[0], pos)

            # The enzyme names are correctly returned:
            self.assertEqual(sources[0].get_enzymes(), ['EcoRI'])

    def test_linear_multiple_restriction(self):

        dseq = Dseqrecord('AAAGGATCCAAAAGATATCAAAAA', circular=False)
        json_seq = format_sequence_genbank(dseq)
        json_seq.id = 1

        # We try EcoRV, which gives blunt ends
        source = RestrictionEnzymeDigestionSource(id=0)
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data, params={'restriction_enzymes': ['BamHI', 'EcoRV']})
        payload = response.json()

        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        sources = [RestrictionEnzymeDigestionSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 3)
        self.assertEqual(len(sources), 3)

        # The right sequences are returned in the right order
        fragment_sequences = 'AAAG^GATCCAAAAGAT^ATCAAAAA'.split('^')
        for i, s in enumerate(fragment_sequences):
            self.assertEqual(resulting_sequences[i].seq.watson.upper(), s)

        # The edges of the fragments are correct:
        # AAAG^GATC_CAAAAGAT^ATCAAAAA
        edges = [
            (None, RestrictionSequenceCut(cut_watson=4, overhang=-4, restriction_enzyme='BamHI')),
            (
                RestrictionSequenceCut(cut_watson=4, overhang=-4, restriction_enzyme='BamHI'),
                RestrictionSequenceCut(cut_watson=16, overhang=0, restriction_enzyme='EcoRV'),
            ),
            (RestrictionSequenceCut(cut_watson=16, overhang=0, restriction_enzyme='EcoRV'), None),
        ]
        for i, e in enumerate(edges):
            self.assertEqual(sources[i].left_edge, e[0])
            self.assertEqual(sources[i].right_edge, e[1])

        # The enzyme names are correct
        restriction_enzymes = [['BamHI'], ['BamHI', 'EcoRV'], ['EcoRV']]
        for i, e in enumerate(restriction_enzymes):
            self.assertEqual(sources[i].get_enzymes(), e)

        # Submitting the known fragments

        for i in range(len(edges)):
            source = RestrictionEnzymeDigestionSource(id=0, left_edge=edges[i][0], right_edge=edges[i][1])
            data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
            response = client.post('/restriction', json=data)
            payload = response.json()

            resulting_sequences = [
                read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
            ]
            sources = [RestrictionEnzymeDigestionSource.model_validate(s) for s in payload['sources']]

            self.assertEqual(len(resulting_sequences), 1)
            self.assertEqual(len(sources), 1)
            self.assertEqual(resulting_sequences[0].seq.watson, fragment_sequences[i])

    def test_circular_multiple_restriction(self):

        dseq = Dseqrecord('AAAGGATCCAAAAGATATCAAAAA', circular=True)
        json_seq = format_sequence_genbank(dseq)
        json_seq.id = 1

        # We try EcoRV, which gives blunt ends
        source = RestrictionEnzymeDigestionSource(id=0)
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data, params={'restriction_enzymes': ['BamHI', 'EcoRV']})
        payload = response.json()

        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        sources = [RestrictionEnzymeDigestionSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 2)
        self.assertEqual(len(sources), 2)

        fragment_sequences = 'GATCCAAAAGAT^ATCAAAAAAAAG'.split('^')
        # The right sequences are returned in the right order
        for i, s in enumerate(fragment_sequences):
            self.assertEqual(resulting_sequences[i].seq.watson.upper(), s)

        # The edges of the fragments are correct:
        # AAAG^GATC_CAAAAGAT^ATCAAAAA
        # We use 8 + len because it's an origin-spanning sequence
        edges = [
            (
                RestrictionSequenceCut(cut_watson=4, overhang=-4, restriction_enzyme='BamHI'),
                RestrictionSequenceCut(cut_watson=16, overhang=0, restriction_enzyme='EcoRV'),
            ),
            (
                RestrictionSequenceCut(cut_watson=16, overhang=0, restriction_enzyme='EcoRV'),
                RestrictionSequenceCut(cut_watson=4, overhang=-4, restriction_enzyme='BamHI'),
            ),
        ]

        for i, e in enumerate(edges):
            self.assertEqual(sources[i].left_edge, e[0])
            self.assertEqual(sources[i].right_edge, e[1])

        # The enzyme names are correct
        restriction_enzymes = [['BamHI', 'EcoRV'], ['EcoRV', 'BamHI']]
        for i, e in enumerate(restriction_enzymes):
            self.assertEqual(sources[i].get_enzymes(), e)

        # Submitting the known fragments
        for i in range(len(edges)):
            source = RestrictionEnzymeDigestionSource(
                id=0,
                left_edge=edges[i][0],
                right_edge=edges[i][1],
            )
            data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
            response = client.post('/restriction', json=data)
            payload = response.json()

            resulting_sequences = [
                read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
            ]
            sources = [RestrictionEnzymeDigestionSource.model_validate(s) for s in payload['sources']]

            self.assertEqual(len(resulting_sequences), 1)
            self.assertEqual(len(sources), 1)
            self.assertEqual(resulting_sequences[0].seq.watson, fragment_sequences[i])


class PCRTest(unittest.TestCase):
    def test_pcr(self):

        template = Dseqrecord(Dseq('TTTTACGTACGTAAAAAAGCGCGCGCTTTTT'))

        json_seq = format_sequence_genbank(template)
        json_seq.id = 1

        submitted_source = PCRSource(id=0)

        primer_fwd = PrimerModel(sequence='ACGTACGT', id=2, name='forward')

        primer_rvs = PrimerModel(sequence='GCGCGCGC', id=3, name='reverse')

        data = {
            'source': submitted_source.model_dump(),
            'sequences': [json_seq.model_dump()],
            'primers': [primer_fwd.model_dump(), primer_rvs.model_dump()],
        }
        response = client.post('/pcr', json=data, params={'minimal_annealing': 8})
        payload = response.json()

        sources = [PCRSource.model_validate(s) for s in payload['sources']]
        sequences = [read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']]

        # Only one possible output
        self.assertEqual(len(sources), 1)
        self.assertEqual(len(sequences), 1)
        dseq1 = sequences[0]
        source1 = sources[0]

        # The sequence matches what we expect
        self.assertEqual(str(dseq1.seq), 'ACGTACGTAAAAAAGCGCGCGC')

        # Now we submit the deterministic PCR (we already know which fragment we want)

        data = {
            'source': source1.model_dump(),
            'sequences': [json_seq.model_dump()],
            'primers': [primer_fwd.model_dump(), primer_rvs.model_dump()],
        }
        response = client.post('/pcr', json=data)
        payload = response.json()

        sources = [PCRSource.model_validate(s) for s in payload['sources']]
        sequences = [read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']]

        # Only one possible output
        self.assertEqual(len(sources), 1)
        self.assertEqual(len(sequences), 1)
        dseq2 = sequences[0]
        source2 = sources[0]
        self.assertEqual(source1, source2)
        self.assertEqual(str(dseq2.seq), 'ACGTACGTAAAAAAGCGCGCGC')

    def test_same_primer_twice(self):
        """
        Special case where you want only one sequence to be returned, since otherwise
        the same sequence is returned as a forward and reverse complement.
        """
        template = Dseqrecord(Dseq('TTTTACGTACGTAAAAAAACGTACGTTTTTT'))

        json_seq = format_sequence_genbank(template)
        json_seq.id = 1

        submitted_source = PCRSource(id=0)

        primer_fwd = PrimerModel(sequence='ACGTACGT', id=2, name='forward')

        data = {
            'source': submitted_source.model_dump(),
            'sequences': [json_seq.model_dump()],
            'primers': [primer_fwd.model_dump(), primer_fwd.model_dump()],
        }

        response = client.post('/pcr', json=data, params={'minimal_annealing': 8})
        payload = response.json()
        self.assertEqual(response.status_code, 200)

        sources = [PCRSource.model_validate(s) for s in payload['sources']]
        sequences = [read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']]
        self.assertEqual(len(sources), 1)
        self.assertEqual(len(sequences), 1)

    def test_wrong_primers(self):

        template = Dseqrecord(Dseq('TTTTACGTACGTAAAAAAGCGCGCGCTTTTT'))
        json_seq = format_sequence_genbank(template)
        json_seq.id = 1

        submitted_source = PCRSource(id=0, input=[1])

        primer_fwd = PrimerModel(sequence='CCCCCCCC', id=2, name='forward')

        primer_rvs = PrimerModel(sequence='GCGCGCGC', id=3, name='reverse')

        # Without specifying the pair
        data = {
            'source': submitted_source.model_dump(),
            'sequences': [json_seq.model_dump()],
            'primers': [primer_fwd.model_dump(), primer_rvs.model_dump()],
        }
        response = client.post('/pcr', json=data, params={'minimal_annealing': 8})
        payload = response.json()
        self.assertEqual(response.status_code, 400)
        self.assertEqual(
            payload['detail'], 'No pair of annealing primers was found. Try changing the annealing settings.'
        )

        # We submit the right pair of primers, but the wrong annealing information

        primer_fwd = PrimerModel(sequence='ACGTACGT', id=2, name='forward')

        submitted_source = PCRSource.from_assembly(
            id=0,
            assembly=[
                (1, 2, SimpleLocation(0, 8), SimpleLocation(4, 12)),
                (2, -3, SimpleLocation(18, 26), SimpleLocation(0, 8)),
            ],
            circular=False,
            fragments=[primer_fwd, json_seq, primer_rvs],
        )

        data = {
            'source': submitted_source.model_dump(),
            'sequences': [json_seq.model_dump()],
            'primers': [primer_fwd.model_dump(), primer_rvs.model_dump()],
        }
        response = client.post('/pcr', json=data)
        payload = response.json()
        self.assertEqual(response.status_code, 400)

        # This is the wrong annealing info
        submitted_source = PCRSource.from_assembly(
            id=0,
            assembly=[
                (2, -3, SimpleLocation(18, 26), SimpleLocation(0, 8)),
                (1, 2, SimpleLocation(0, 8), SimpleLocation(4, 12)),
            ],
            circular=False,
            fragments=[primer_fwd, json_seq, primer_rvs],
        )

        data = {
            'source': submitted_source.model_dump(),
            'sequences': [json_seq.model_dump()],
            'primers': [primer_fwd.model_dump(), primer_rvs.model_dump()],
        }
        response = client.post('/pcr', json=data)
        payload = response.json()
        self.assertEqual(response.status_code, 400)
        self.assertEqual(payload['detail'], 'The provided assembly is not valid.')

    def test_wrong_stoichiometry(self):
        # If not 2 primers per sequence, bad request
        template = Dseqrecord(Dseq('ACGTACGTGCGCGCGC'))

        json_seq = format_sequence_genbank(template)
        submitted_source = PCRSource(id=0)

        primer_fwd = PrimerModel(sequence='ACGTACGTG', id=2, name='forward')

        data = {
            'source': submitted_source.model_dump(),
            'sequences': [json_seq.model_dump()],
            'primers': [primer_fwd.model_dump()],
        }
        response = client.post('/pcr', json=data, params={'minimal_annealing': 8})
        self.assertEqual(response.status_code, 400)

    def test_too_many_assemblies(self):
        # Too many assemblies
        template = Dseqrecord(Dseq('A' * 200))
        json_template = format_sequence_genbank(template)
        json_template.id = 1

        source = PCRSource(id=0)
        primers = [
            PrimerModel(sequence='A' * 20, id=1, name=f'primer{1}'),
            PrimerModel(sequence='T' * 20, id=2, name=f'primer{2}'),
        ]

        data = {
            'source': source.model_dump(),
            'sequences': [json_template.model_dump()],
            'primers': [primer.model_dump() for primer in primers],
        }
        response = client.post('/pcr', json=data, params={'minimal_annealing': 8})
        self.assertEqual(response.status_code, 400)
        self.assertTrue('Too many assemblies' in response.json()['detail'])


class OligoHybridizationTest(unittest.TestCase):

    def test_hybridization(self):
        default_example = request_examples.oligonucleotide_hybridization_examples['default']['value']
        watson_sequence = default_example['primers'][0]['sequence']
        crick_sequence = default_example['primers'][1]['sequence']
        response = client.post('/oligonucleotide_hybridization', json=default_example)
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['sources']), 1)
        self.assertEqual(len(payload['sequences']), 1)

        source = OligoHybridizationSource.model_validate(payload['sources'][0])
        sequence = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))

        self.assertEqual(source.overhang_crick_3prime, sequence.seq.ovhg)
        self.assertEqual(sequence.seq.watson, watson_sequence.upper())
        self.assertEqual(sequence.seq.crick, crick_sequence.upper())

        # The minimal_annealing parameter is ignored if the overhangs are specified
        modified_example = copy.deepcopy(default_example)
        modified_example['source']['overhang_crick_3prime'] = source.overhang_crick_3prime

        response = client.post(
            '/oligonucleotide_hybridization', json=modified_example, params={'minimal_annealing': 60}
        )
        payload2 = response.json()
        self.assertEqual(response.status_code, 200)
        self.assertEqual(payload, payload2)

    def test_no_valid_result(self):
        default_example = request_examples.oligonucleotide_hybridization_examples['default']['value']
        response = client.post(
            '/oligonucleotide_hybridization', json=default_example, params={'minimal_annealing': 60}
        )
        self.assertEqual(response.status_code, 400)

    def test_invalid_oligo_id(self):
        invalid_oligo_example = copy.deepcopy(
            request_examples.oligonucleotide_hybridization_examples['default']['value']
        )
        invalid_oligo_example['source']['forward_oligo'] = 5
        response = client.post('/oligonucleotide_hybridization', json=invalid_oligo_example)
        self.assertEqual(response.status_code, 404)

    def test_hybridization_in_the_middle(self):
        default_example = request_examples.oligonucleotide_hybridization_examples['default']['value']
        default_example['primers'][0]['sequence'] = 'aaaACGGCAGCCCGTaaa'
        default_example['primers'][1]['sequence'] = 'cccTGCCGTCGGGCAccc'[::-1]
        response = client.post(
            '/oligonucleotide_hybridization', json=default_example, params={'minimal_annealing': 12}
        )
        self.assertEqual(response.status_code, 400)
        payload = response.json()
        self.assertIn('mismatch', payload['detail'])


class HomologousRecombinationTest(unittest.TestCase):
    def test_homologous_recombination(self):
        template = Dseqrecord('TTTTacgatAAtgctccCCCC', circular=False)
        json_template = format_sequence_genbank(template)
        json_template.id = 1

        insert = Dseqrecord('acgatCCCtgctcc', circular=False)
        json_insert = format_sequence_genbank(insert)
        json_insert.id = 2
        # One enzyme
        source = HomologousRecombinationSource(id=0)
        data = {'source': source.model_dump(), 'sequences': [json_template.model_dump(), json_insert.model_dump()]}
        response = client.post('/homologous_recombination', params={'minimal_homology': 5}, json=data)
        self.assertEqual(response.status_code, 200)

        payload = response.json()
        sequences = [read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']]
        self.assertEqual(len(sequences), 1)
        self.assertEqual(str(sequences[0].seq), 'TTTTacgatCCCtgctccCCCC'.upper())

        response = client.post(
            '/homologous_recombination',
            json={
                'source': payload['sources'][0],
                'sequences': [json_template.model_dump(), json_insert.model_dump()],
            },
        )
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        sequences = [read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']]
        self.assertEqual(len(sequences), 1)
        self.assertEqual(str(sequences[0].seq), 'TTTTacgatCCCtgctccCCCC'.upper())

    def test_multiple_insertions(self):
        homology = 'ATGCAAACAGTAATGATGGATGACATTCAAAGCACTGATT'
        template = Dseqrecord(f'aaaaaa{homology}aattggaa{homology}tttttttt', circular=False)
        insert = Dseqrecord(f'{homology}acaa{homology}', circular=False)

        json_template = format_sequence_genbank(template)
        json_template.id = 1
        json_insert = format_sequence_genbank(insert)
        json_insert.id = 2

        source = HomologousRecombinationSource(id=0)
        data = {'source': source.model_dump(), 'sequences': [json_template.model_dump(), json_insert.model_dump()]}

        response = client.post('/homologous_recombination', json=data)
        self.assertEqual(response.status_code, 200)
        payload = response.json()

        sequences = [read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']]
        self.assertEqual(len(sequences), 3)

    def test_too_many_assemblies(self):
        template = Dseqrecord(Dseq('A' * 200))
        insert = Dseqrecord(Dseq('A' * 40))
        json_template = format_sequence_genbank(template)
        json_template.id = 1
        json_insert = format_sequence_genbank(insert)
        json_insert.id = 2

        source = HomologousRecombinationSource(id=0)
        data = {'source': source.model_dump(), 'sequences': [json_template.model_dump(), json_insert.model_dump()]}

        response = client.post('/homologous_recombination', json=data)
        self.assertEqual(response.status_code, 400)
        self.assertTrue('Too many assemblies' in response.json()['detail'])


class GibsonAssemblyTest(unittest.TestCase):
    def test_gibson_assembly(self):
        # A circular one
        fragments = [
            Dseqrecord('TTTTacgatAAtgctccCCCC', circular=False),
            Dseqrecord('CCCCtcatGGGG', circular=False),
            Dseqrecord('GGGGatataTTTT', circular=False),
        ]

        json_fragments = [format_sequence_genbank(f) for f in fragments]
        for i, f in enumerate(json_fragments):
            f.id = i + 1

        source = GibsonAssemblySource(id=0)

        data = {'source': source.model_dump(), 'sequences': [f.model_dump() for f in json_fragments]}
        response = client.post('/gibson_assembly', json=data, params={'minimal_homology': 4})
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        sequences = [read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']]

        self.assertEqual(len(sequences), 2)
        self.assertEqual(str(sequences[0].seq), 'TTTTacgatAAtgctccCCCCtcatGGGGatata'.upper())
        self.assertEqual(str(sequences[1].seq), 'TTTTacgatAAtgctccCCCCatgaGGGGatata'.upper())

        # Circularisation works
        f1 = Dseqrecord('AGAGACCaaaAGAGACC')
        json_fragments = [format_sequence_genbank(f1)]
        for i, f in enumerate(json_fragments):
            f.id = i + 1

        source = GibsonAssemblySource(id=0)
        data = {'source': source.model_dump(), 'sequences': [f.model_dump() for f in json_fragments]}
        response = client.post('/gibson_assembly', json=data, params={'minimal_homology': 7})

        self.assertEqual(response.status_code, 200)
        payload = response.json()
        sequences = [read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']]

        self.assertEqual(len(sequences), 1)
        self.assertEqual(str(sequences[0].seq), 'AGAGACCaaa'.upper())

    def test_circular_constrain(self):
        # A circular one
        fragments = [Dseqrecord('TTTTacgatAAtgctccCCCC', circular=False), Dseqrecord('CCCCtcatGGGG', circular=False)]
        json_fragments = [format_sequence_genbank(f) for f in fragments]
        for i, f in enumerate(json_fragments):
            f.id = i + 1

        source = GibsonAssemblySource(id=0)

        data = {'source': source.model_dump(), 'sequences': [f.model_dump() for f in json_fragments]}
        response = client.post('/gibson_assembly', json=data, params={'minimal_homology': 4, 'circular_only': True})
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json()['detail'], 'No circular assembly with at least 4 bps of homology was found.')

    def test_too_many_assemblies(self):
        homology = 'ATGCAAACAGTAATGATGGATGACATTCAAAGCACTGATT'
        fragments = [Dseqrecord(f'{homology}acgt{homology}') for _ in range(10)]
        json_fragments = [format_sequence_genbank(f) for f in fragments]
        for i, f in enumerate(json_fragments):
            f.id = i + 1

        source = GibsonAssemblySource(id=0)
        data = {'source': source.model_dump(), 'sequences': [f.model_dump() for f in json_fragments]}
        response = client.post('/gibson_assembly', json=data)
        self.assertEqual(response.status_code, 400)
        self.assertTrue('Too many possible paths' in response.json()['detail'])


class RestrictionAndLigationTest(unittest.TestCase):
    def test_restriction_and_ligation(self):
        fragments = [Dseqrecord('AAAGAATTCAAA'), Dseqrecord('CCCCGAATTCCCC')]
        json_fragments = [format_sequence_genbank(f) for f in fragments]
        for i, f in enumerate(json_fragments):
            f.id = i + 1

        source = RestrictionAndLigationSource(
            id=0,
            restriction_enzymes=['EcoRI'],
        )

        data = {'source': source.model_dump(), 'sequences': [f.model_dump() for f in json_fragments]}
        response = client.post('/restriction_and_ligation', json=data, params={'minimal_homology': 4})

        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['sequences']), 4)
        self.assertEqual(len(payload['sources']), 4)

        # Test with blunt ends

        fragments = [Dseqrecord('cccAGCGCTcgcAGCTtat'), Dseqrecord('aaaAGCGCTggaAGCTctt')]
        json_fragments = [format_sequence_genbank(f) for f in fragments]
        for i, f in enumerate(json_fragments):
            f.id = i + 1

        source = RestrictionAndLigationSource(
            id=0,
            restriction_enzymes=['AfeI', 'AluI'],
        )
        data = {'source': source.model_dump(), 'sequences': [f.model_dump() for f in json_fragments]}
        response = client.post('/restriction_and_ligation', json=data, params={'circular_only': True})

        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['sequences']), 2)
        self.assertEqual(len(payload['sources']), 2)

        seqs = [read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']]
        self.assertEqual(str(seqs[0].seq), 'GCTcgcAGGCTggaAG'.upper())
        self.assertEqual(str(seqs[1].seq), 'GCTcgcAGCTtccAGC'.upper())

    def test_golden_gate(self):
        fragments = [
            Dseqrecord('GGTCTCAattaAAAAAttaaAGAGACC'),
            Dseqrecord('GGTCTCAttaaCCCCCatatAGAGACC'),
            Dseqrecord('GGTCTCAatatGGGGGccggAGAGACC'),
            Dseqrecord('TTTTattaAGAGACCTTTTTGGTCTCAccggTTTT', circular=True),
        ]

        json_fragments = [format_sequence_genbank(f) for f in fragments]
        for i, f in enumerate(json_fragments):
            f.id = i + 1

        source = RestrictionAndLigationSource(
            id=0,
            restriction_enzymes=['BsaI'],
        )

        data = {'source': source.model_dump(), 'sequences': [f.model_dump() for f in json_fragments]}
        response = client.post('/restriction_and_ligation', json=data, params={'circular_only': True})

        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['sequences']), 1)
        self.assertEqual(len(payload['sources']), 1)

    def test_single_input(self):
        fragments = [Dseqrecord('AAAGAATTCAAAGAATTCAAAA')]
        json_fragments = [format_sequence_genbank(f) for f in fragments]
        for i, f in enumerate(json_fragments):
            f.id = i + 1

        source = RestrictionAndLigationSource(
            id=0,
            restriction_enzymes=['EcoRI'],
        )

        data = {'source': source.model_dump(), 'sequences': [f.model_dump() for f in json_fragments]}
        response = client.post('/restriction_and_ligation', json=data)

        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['sequences']), 2)
        self.assertEqual(len(payload['sources']), 2)

        sequences = [read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']]

        self.assertEqual(str(sequences[0].seq), 'AATTCAAAG')
        self.assertEqual(str(sequences[1].seq), 'AAAGAATTCAAAA')

        sources = [RestrictionAndLigationSource.model_validate(s) for s in payload['sources']]
        # Submitting the known fragments
        data = {'source': sources[0].model_dump(), 'sequences': [f.model_dump() for f in json_fragments]}
        response = client.post('/restriction_and_ligation', json=data)
        payload = response.json()
        self.assertEqual(response.status_code, 200)
        self.assertEqual(len(payload['sequences']), 1)
        self.assertEqual(len(payload['sources']), 1)
        self.assertEqual(payload['sources'][0], sources[0].model_dump())
        self.assertEqual(
            read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0])), sequences[0]
        )

    def test_errors(self):
        fragments = [Dseqrecord('AAAGAATTCAAAGAATTCAAAA')]
        json_fragments = [format_sequence_genbank(f) for f in fragments]
        for i, f in enumerate(json_fragments):
            f.id = i + 1

        # Enzyme that does not exist
        source = RestrictionAndLigationSource(
            id=0,
            restriction_enzymes=['dummy'],
        )

        data = {'source': source.model_dump(), 'sequences': [f.model_dump() for f in json_fragments]}
        response = client.post('/restriction_and_ligation', json=data)
        self.assertEqual(response.status_code, 404)

        # enzyme that does not cut
        source = RestrictionAndLigationSource(
            id=0,
            restriction_enzymes=['HindIII'],
        )

        data = {'source': source.model_dump(), 'sequences': [f.model_dump() for f in json_fragments]}
        response = client.post('/restriction_and_ligation', json=data)

    def test_too_many_assemblies(self):
        fragments = [
            Dseqrecord('aaGCGGCCGCaaGCGGCCGC', circular=True),
            Dseqrecord('aaGCGGCCGCaaGCGGCCGC', circular=True),
            Dseqrecord('aaGCGGCCGCaaGCGGCCGCaaGCGGCCGCaaGCGGCCGC', circular=True),
        ]

        source = RestrictionAndLigationSource(
            id=0,
            restriction_enzymes=['NotI'],
        )

        json_fragments = [format_sequence_genbank(f) for f in fragments]
        for i, f in enumerate(json_fragments):
            f.id = i + 1

        data = {'source': source.model_dump(), 'sequences': [f.model_dump() for f in json_fragments]}
        response = client.post('/restriction_and_ligation', json=data)
        self.assertEqual(response.status_code, 400)
        self.assertTrue('Too many assemblies' in response.json()['detail'])


class ManuallyTypedTest(unittest.TestCase):
    def test_manually_typed(self):
        """Test the manually_typed endpoint"""

        # Test linear (default)
        source = ManuallyTypedSource(
            id=0,
            user_input='ATGC',
        )

        response = client.post('/manually_typed', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        sources = [ManuallyTypedSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)
        self.assertEqual(str(resulting_sequences[0].seq), 'ATGC')
        self.assertEqual(sources[0], source)

        # Test circular
        source.circular = True
        response = client.post('/manually_typed', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        sources = [ManuallyTypedSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)
        self.assertEqual(str(resulting_sequences[0].seq), 'ATGC')
        self.assertEqual(resulting_sequences[0].seq.circular, True)
        self.assertEqual(sources[0], source)

        # Test overhangs
        source = ManuallyTypedSource(
            id=0,
            user_input='ATGC',
            overhang_crick_3prime=1,
            overhang_watson_3prime=2,
        )

        response = client.post('/manually_typed', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(resulting_sequences[0].seq, Dseq.from_full_sequence_and_overhangs('ATGC', 1, 2))

        # Test that if overhangs are set, it cannot be circular
        wrong_source = source.model_dump()
        wrong_source['circular'] = True

        response = client.post('/manually_typed', json=wrong_source)
        self.assertEqual(response.status_code, 422)

    # Test that it fails if not acgt or empty
    def test_manually_typed_fail(self):

        response = client.post('/manually_typed', json={'user_input': 'ATGZ'})
        self.assertEqual(response.status_code, 422)

        response = client.post('/manually_typed', json={'user_input': ''})
        self.assertEqual(response.status_code, 422)


class GenomeRegionTest(unittest.TestCase):
    def assertStatusCode(self, response_code: int, expected: int, msg: str = ''):
        if response_code == 503:
            self.skipTest('NCBI not available')
        else:
            self.assertEqual(response_code, expected, msg)

    @pytest.mark.flaky(reruns=3, reruns_delay=2)
    def test_examples(self):
        for example_name in request_examples.genome_region_examples:
            example = request_examples.genome_region_examples[example_name]
            response = client.post('/genome_coordinates', json=example['value'])
            msg = f'Error in example {example_name}'
            self.assertStatusCode(response.status_code, 200, msg)
            payload = response.json()
            try:
                sources = [GenomeCoordinatesSource.model_validate(s) for s in payload['sources']]
            except Exception:
                self.fail(f'Cannot parse the sources for example {example_name}')
            response_source = sources[0]
            request_source = GenomeCoordinatesSource.model_validate(example['value'])
            if example_name == 'full' or example_name == 'full_with_genbank_accession':
                self.assertEqual(response_source, request_source, msg)
            elif example_name == 'id_omitted':
                request_source.gene_id = 2543372
                self.assertEqual(response_source, request_source, msg)
            elif example_name == 'assembly_accession_omitted':
                self.assertEqual(response_source, request_source, msg)
            elif example_name == 'viral_sequence':
                self.assertEqual(response_source, request_source, msg)

    @pytest.mark.flaky(reruns=3, reruns_delay=2)
    def test_exceptions(self):
        # Load first example
        correct_source = GenomeCoordinatesSource.model_validate(
            request_examples.genome_region_examples['full']['value']
        )

        # Ommit assembly accession
        s = correct_source.model_copy()
        s.assembly_accession = None
        response = client.post('/genome_coordinates', json=s.model_dump())
        self.assertStatusCode(response.status_code, 422)

        # Ommit locus_tag keeping gene id is now supported
        # Before it was not, but see https://github.com/ncbi/datasets/issues/397
        s = correct_source.model_copy()
        s.locus_tag = None
        response = client.post('/genome_coordinates', json=s.model_dump())
        self.assertStatusCode(response.status_code, 200)

        # Wrong gene_id (not matching that of the locus_tag)
        s = correct_source.model_copy()
        s.gene_id = 123
        response = client.post('/genome_coordinates', json=s.model_dump())
        self.assertStatusCode(response.status_code, 400)

        # Wrong assembly accession
        s = correct_source.model_copy()
        s.assembly_accession = 'blah'
        response = client.post('/genome_coordinates', json=s.model_dump())
        self.assertStatusCode(response.status_code, 404)

        # Wrong locus_tag
        s = correct_source.model_copy()
        s.locus_tag = 'blah'
        response = client.post('/genome_coordinates', json=s.model_dump())
        self.assertStatusCode(response.status_code, 404)

        # Wrong coordinates
        s = correct_source.model_copy()
        s.start = 1
        s.end = 10
        response = client.post('/genome_coordinates', json=s.model_dump())
        self.assertStatusCode(response.status_code, 400)

        # Wrong assembly accession
        s = correct_source.model_copy()
        s.locus_tag = None
        s.gene_id = None
        s.assembly_accession = 'blah'
        response = client.post('/genome_coordinates', json=s.model_dump())
        self.assertStatusCode(response.status_code, 404)
        self.assertIn('Wrong assembly accession', response.json()['detail'])

        # Assembly accession not linked to any sequence record
        s = correct_source.model_copy()
        s.locus_tag = None
        s.gene_id = None
        s.assembly_accession = 'GCF_000146045.1'
        response = client.post('/genome_coordinates', json=s.model_dump())
        self.assertStatusCode(response.status_code, 400)
        self.assertIn('No sequence accessions linked', response.json()['detail'])

        # Assembly accession not linked to that sequence accession
        s = correct_source.model_copy()
        s.locus_tag = None
        s.gene_id = None
        s.assembly_accession = 'GCF_000146045.2'
        response = client.post('/genome_coordinates', json=s.model_dump())
        self.assertStatusCode(response.status_code, 400)
        self.assertIn('not contained in assembly accession', response.json()['detail'])

        # Wrong sequence accession
        s = correct_source.model_copy()
        s.locus_tag = None
        s.gene_id = None
        s.assembly_accession = None
        s.sequence_accession = 'blah'
        response = client.post('/genome_coordinates', json=s.model_dump())
        self.assertStatusCode(response.status_code, 404)

        # Coordinates malformatted
        viral_source = GenomeCoordinatesSource.model_validate(
            request_examples.genome_region_examples['viral_sequence']['value']
        )
        viral_source.start = 10
        viral_source.end = 1
        response = client.post('/genome_coordinates', json=viral_source.model_dump())
        self.assertStatusCode(response.status_code, 422)

        viral_source.start = 0
        viral_source.end = 20
        response = client.post('/genome_coordinates', json=viral_source.model_dump())
        self.assertStatusCode(response.status_code, 422)

        viral_source.start = 1
        viral_source.end = 20
        viral_source.strand = 0
        response = client.post('/genome_coordinates', json=viral_source.model_dump())
        self.assertStatusCode(response.status_code, 422)

        # Coordinates outside of the sequence
        viral_source.start = 1
        # the length is 2151
        viral_source.end = 2152
        viral_source.strand = 1
        response = client.post('/genome_coordinates', json=viral_source.model_dump())
        self.assertStatusCode(response.status_code, 400)

        # Coordinates too long
        viral_source.start = 1
        viral_source.end = 100004
        response = client.post('/genome_coordinates', json=viral_source.model_dump())
        self.assertStatusCode(response.status_code, 400)


class PolymeraseExtensionTest(unittest.TestCase):

    def test_polymerase_extension(self):
        dseq = Dseq.from_full_sequence_and_overhangs('ACGTT', 1, 1)
        template = Dseqrecord(dseq, circular=False)
        json_template = format_sequence_genbank(template)
        json_template.id = 1

        source = PolymeraseExtensionSource(
            id=2,
        )

        data = {'source': source.model_dump(), 'sequences': [json_template.model_dump()]}
        response = client.post('/polymerase_extension', json=data)
        payload = response.json()
        self.assertEqual(response.status_code, 200)
        sequences = [read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']]

        self.assertEqual(len(sequences), 1)
        self.assertEqual(sequences[0].seq, dseq.fill_in())
        response_source = PolymeraseExtensionSource.model_validate(payload['sources'][0])
        self.assertEqual(response_source, source)

    def test_exceptions(self):

        # Sequence without overhangs
        cases = [
            (Dseqrecord('ACGTT'), 400, 1),  # No overhangs
            (Dseqrecord('ACGTT', circular=True), 400, 1),  # circular
        ]

        for template, status_code, id in cases:
            json_template = format_sequence_genbank(template)
            json_template.id = id

            source = PolymeraseExtensionSource(
                id=2,
            )

            data = {'source': source.model_dump(), 'sequences': [json_template.model_dump()]}
            response = client.post('/polymerase_extension', json=data)
            self.assertEqual(response.status_code, status_code)


class CrisprTest(unittest.TestCase):

    def test_crispr(self):

        template = Dseqrecord('aaccggttcaatgcaaacagtaatgatggatgacattcaaagcac')
        json_template = format_sequence_genbank(template)
        json_template.id = 1

        insert = Dseqrecord('aaccggttAAAAAAAAAttcaaagcac')
        json_insert = format_sequence_genbank(insert)
        json_insert.id = 2

        guide = PrimerModel(sequence='ttcaatgcaaacagtaatga', id=3, name='guide_1')
        source = CRISPRSource(
            id=0,
            guides=[3],
        )
        data = {
            'source': source.model_dump(),
            'sequences': [json_template.model_dump(), json_insert.model_dump()],
            'guides': [guide.model_dump()],
        }
        params = {'minimal_homology': 8}

        response = client.post('/crispr', json=data, params=params)

        payload = response.json()
        self.assertEqual(response.status_code, 200)
        resulting_sequences = [
            read_dsrecord_from_json(TextFileSequence.model_validate(s)) for s in payload['sequences']
        ]
        sources = [CRISPRSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)

        self.assertEqual(str(resulting_sequences[0].seq), 'aaccggttAAAAAAAAAttcaaagcac'.upper())

    def test_errors(self):

        # Wrong guide
        template = Dseqrecord('aaccggttcaatgcaaacagtaatgatggatgacattcaaagcac')
        json_template = format_sequence_genbank(template)
        json_template.id = 1

        insert = Dseqrecord('aaccggttAAAAAAAAAttcaaagcac')
        json_insert = format_sequence_genbank(insert)
        json_insert.id = 2

        guide = PrimerModel(sequence='AAAAAAAA', id=3, name='guide_1')

        source = CRISPRSource(id=0, guides=[3])

        data = {
            'source': source.model_dump(),
            'sequences': [json_template.model_dump(), json_insert.model_dump()],
            'guides': [guide.model_dump()],
        }

        params = {'minimal_homology': 8}

        response = client.post('/crispr', json=data, params=params)

        payload = response.json()
        self.assertEqual(response.status_code, 400)

        self.assertIn('Could not find Cas9 cutsite in the target sequence using the guide: guide_1', payload['detail'])

        # Wrong and right guide
        guide2 = PrimerModel(sequence='ttcaatgcaaacagtaatga', id=4, name='guide_2')
        data['guides'] = [guide.model_dump(), guide2.model_dump()]

        response = client.post('/crispr', json=data, params=params)
        payload = response.json()
        self.assertEqual(response.status_code, 400)

        self.assertIn('Could not find Cas9 cutsite in the target sequence using the guide: guide_1', payload['detail'])

        # Homology too short
        data['guides'] = [guide2.model_dump()]
        response = client.post('/crispr', json=data)
        payload = response.json()
        self.assertEqual(response.status_code, 400)

        # cut outside of homology

        template2 = Dseqrecord('aaccggttcaatgcaaacagtaatgatggatgacattcaaagcacaaaAAAGAGAGACAGGTTTTGAGtgg')
        json_template2 = format_sequence_genbank(template2)
        json_template2.id = 1

        guide_outside = PrimerModel(sequence='AAAGAGAGACAGGTTTTGAG', id=3, name='guide_outside')
        data = {
            'source': source.model_dump(),
            'sequences': [json_template2.model_dump(), json_insert.model_dump()],
            'guides': [guide_outside.model_dump()],
        }

        response = client.post('/crispr', json=data, params=params)
        payload = response.json()
        self.assertEqual(response.status_code, 400)
        self.assertEqual(
            payload['detail'],
            'A Cas9 cutsite was found, and a homologous recombination region, but they do not overlap.',
        )

    def test_too_many_assemblies(self):
        template = Dseqrecord(30 * 'A' + 'aattcaatgcaaacagtaatgatggatgaca' + 30 * 'A')
        json_template = format_sequence_genbank(template)
        json_template.id = 1

        insert = Dseqrecord(30 * 'A' + 'AAAAAAAAA' + 30 * 'A')
        json_insert = format_sequence_genbank(insert)
        json_insert.id = 2

        guide = PrimerModel(sequence='ttcaatgcaaacagtaatga', id=3, name='guide_1')
        source = CRISPRSource(
            id=0,
            guides=[3],
        )
        data = {
            'source': source.model_dump(),
            'sequences': [json_template.model_dump(), json_insert.model_dump()],
            'guides': [guide.model_dump()],
        }
        params = {'minimal_homology': 8}

        response = client.post('/crispr', json=data, params=params)

        self.assertEqual(response.status_code, 400)
        self.assertIn('Too many assemblies', response.json()['detail'])


class ValidatorTest(unittest.TestCase):
    def test_validator(self):
        with open('test_files/homologous_recombination.json') as ins:
            # Read it to json
            data = json.load(ins)
        BaseCloningStrategy.model_validate(data)


class PrimerDesignTest(unittest.TestCase):

    def test_homologous_recombination(self):
        pcr_seq = format_sequence_genbank(Dseqrecord('GAAATGGAACAGTGCCAGAAATTTTT'))
        pcr_seq.id = 1
        pcr_loc = PydanticSimpleLocation(start=1, end=23)
        hr_seq = format_sequence_genbank(Dseqrecord('AAACGTTT'))
        hr_seq.id = 2
        hr_loc_replace = PydanticSimpleLocation(start=3, end=5)

        homology_length = 3
        minimal_hybridization_length = 10
        insert_forward = True
        target_tm = 30

        # First we replace the CG
        data = {
            'pcr_template': {
                'sequence': pcr_seq.model_dump(),
                'location': pcr_loc.model_dump(),
            },
            'homologous_recombination_target': {
                'sequence': hr_seq.model_dump(),
                'location': hr_loc_replace.model_dump(),
            },
        }
        params = {
            'homology_length': homology_length,
            'minimal_hybridization_length': minimal_hybridization_length,
            'insert_forward': insert_forward,
            'target_tm': target_tm,
        }
        response = client.post('/primer_design/homologous_recombination', json=data, params=params)
        payload = response.json()

        self.assertEqual(response.status_code, 200)
        self.assertEqual(payload['primers'][0]['sequence'], 'aaaAAATGGAACAG')
        self.assertEqual(payload['primers'][1]['sequence'], 'aaaAATTTCTGGC')

        # Raise valuerror
        params['homology_length'] = 10
        response = client.post('/primer_design/homologous_recombination', json=data, params=params)
        self.assertEqual(response.status_code, 400)
        payload = response.json()
        self.assertEqual(payload['detail'], 'Forward homology region is out of bounds.')

        # Test an insertion with spacers and reversed insert
        params['homology_length'] = 3
        data['homologous_recombination_target']['forward_orientation'] = False
        data['homologous_recombination_target']['location'] = PydanticSimpleLocation(start=3, end=3).model_dump()
        data['spacers'] = ['attt', 'cggg']
        response = client.post('/primer_design/homologous_recombination', json=data, params=params)
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(payload['primers'][0]['sequence'], 'aaaatttAAATGGAACAG')
        self.assertEqual(payload['primers'][1]['sequence'], 'acgcccgAATTTCTGGC')

        # Raise error if the number of spacers is incorrect
        data['spacers'] = ['attt', 'cggg', 'tttt']
        response = client.post('/primer_design/homologous_recombination', json=data, params=params)
        self.assertEqual(response.status_code, 422)
        self.assertIn('The number of spacers must be', response.json()['detail'])

        # Raise error if the spacer is not DNA
        data['spacers'] = ['zzz', 'cggg']
        response = client.post('/primer_design/homologous_recombination', json=data, params=params)
        self.assertEqual(response.status_code, 422)
        self.assertIn('Spacer can only contain ACGT bases', response.json()['detail'])

    def test_gibson_assembly(self):
        # Test case for gibson_assembly_primers endpoint
        templates = [
            Dseqrecord('AAACAGTAATACGTTCCTTTTTTATGATGATGGATGACATTCAAAGCACTGATTCTAT'),
            Dseqrecord('GTTTACAACGGCAATGAACGTTCCTTTTTTATGATATGCCCAGCTTCATGAAATGGAA'),
            Dseqrecord('AAGGACAACGTTCCTTTTTTATGATATATATGGCACAGTATGATCAAAAGTTAAGTAC'),
        ]

        queries = []
        for i, template in enumerate(templates):
            json_seq = format_sequence_genbank(template)
            json_seq.id = i
            queries.append(
                {
                    'sequence': json_seq.model_dump(),
                    'location': PydanticSimpleLocation(start=0, end=len(template)).model_dump(),
                    'forward_orientation': True,
                }
            )

        params = {'homology_length': 20, 'minimal_hybridization_length': 15, 'target_tm': 55, 'circular': True}

        response = client.post(
            '/primer_design/gibson_assembly', json={'pcr_templates': queries, 'spacers': None}, params=params
        )

        self.assertEqual(response.status_code, 200)

        payload = response.json()
        self.assertEqual(len(payload['primers']), 6)  # 2 primers per template
        for i, p in enumerate(payload['primers']):
            p = PrimerModel.model_validate(p)
            # Check that the name is correct
            if i % 2 == 0:
                self.assertEqual(p.name, f'seq_{i//2}_fwd')
            else:
                self.assertEqual(p.name, f'seq_{i//2}_rvs')

        # Primer naming also work for named sequences
        for i, t in enumerate(templates):
            t.name = f'template_{i}'

        queries = []
        for i, template in enumerate(templates):
            json_seq = format_sequence_genbank(template)
            json_seq.id = i
            queries.append(
                {
                    'sequence': json_seq.model_dump(),
                    'location': PydanticSimpleLocation(start=0, end=len(template)).model_dump(),
                    'forward_orientation': True,
                }
            )

        params = {'homology_length': 20, 'minimal_hybridization_length': 15, 'target_tm': 55, 'circular': True}

        response = client.post(
            '/primer_design/gibson_assembly', json={'pcr_templates': queries, 'spacers': None}, params=params
        )

        self.assertEqual(response.status_code, 200)

        payload = response.json()
        self.assertEqual(len(payload['primers']), 6)  # 2 primers per template
        for i, p in enumerate(payload['primers']):
            p = PrimerModel.model_validate(p)
            # Check that the name is correct
            if i % 2 == 0:
                self.assertEqual(p.name, f'template_{i//2}_fwd')
            else:
                self.assertEqual(p.name, f'template_{i//2}_rvs')

        # Test error case with invalid parameters
        params['minimal_hybridization_length'] = 100  # Too long
        response = client.post(
            '/primer_design/gibson_assembly', json={'pcr_templates': queries, 'spacers': None}, params=params
        )
        self.assertEqual(response.status_code, 400)
        self.assertIn('Primers could not be designed', response.json()['detail'])

        # Test case with spacers
        params['homology_length'] = 20
        params['circular'] = True
        params['minimal_hybridization_length'] = 15
        spacers = ['aaaa', 'tttt', 'cccc']
        response = client.post(
            '/primer_design/gibson_assembly', json={'pcr_templates': queries, 'spacers': spacers}, params=params
        )
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['primers']), 6)  # 2 primers per template
        primers = [PrimerModel.model_validate(p) for p in payload['primers']]
        self.assertTrue(primers[0].sequence.startswith('TTAAGTACccccAAACAGTA'))
        self.assertTrue(reverse_complement(primers[-1].sequence).endswith('TTAAGTACccccAAACAGTA'))
        self.assertTrue(reverse_complement(primers[1].sequence).endswith('GATTCTATaaaaGTTTACAA'))
        self.assertTrue(primers[2].sequence.startswith('GATTCTATaaaaGTTTACAA'))
        self.assertTrue(reverse_complement(primers[3].sequence).endswith('AAATGGAAttttAAGGACAA'))
        self.assertTrue(primers[4].sequence.startswith('AAATGGAAttttAAGGACAA'))

        # Test that wrong number of spacers fails
        response = client.post(
            '/primer_design/gibson_assembly',
            json={'pcr_templates': queries, 'spacers': ['aaaa', 'tttt']},
            params=params,
        )
        self.assertEqual(response.status_code, 422)

        # Test that non-DNA spacers fails
        response = client.post(
            '/primer_design/gibson_assembly',
            json={'pcr_templates': queries, 'spacers': ['hello', 'TTTT']},
            params=params,
        )
        self.assertEqual(response.status_code, 422)

    def test_simple_pair(self):
        from Bio.Restriction import EcoRI, BamHI

        # Create a test sequence
        dseqr = Dseqrecord('aaaGGCTTCACCAAGTCCTTGGAACAGccc')
        dseqr.name = 'test_sequence'
        json_seq = format_sequence_genbank(dseqr)
        json_seq.id = 0

        query = {
            'sequence': json_seq.model_dump(),
            'location': PydanticSimpleLocation(start=3, end=27).model_dump(),
            'forward_orientation': True,
        }

        params = {
            'minimal_hybridization_length': 10,
            'target_tm': 30,
            'left_enzyme': 'EcoRI',
            'right_enzyme': 'BamHI',
            'filler_bases': 'GC',
        }

        response = client.post('/primer_design/simple_pair', json={'pcr_template': query}, params=params)
        self.assertEqual(response.status_code, 200)

        payload = response.json()
        self.assertEqual(len(payload['primers']), 2)  # 2 primers (forward and reverse)

        fwd_primer = PrimerModel.model_validate(payload['primers'][0])
        rvs_primer = PrimerModel.model_validate(payload['primers'][1])

        self.assertEqual(fwd_primer.name, 'test_sequence_EcoRI_fwd')
        self.assertEqual(rvs_primer.name, 'test_sequence_BamHI_rvs')

        self.assertTrue(fwd_primer.sequence.startswith('GC' + str(EcoRI.site)))
        self.assertTrue(rvs_primer.sequence.startswith('GC' + str(BamHI.site)))

        # Same without enzymes
        params_no_enzymes = copy.deepcopy(params)
        params_no_enzymes.pop('left_enzyme')
        params_no_enzymes.pop('right_enzyme')
        response = client.post('/primer_design/simple_pair', json={'pcr_template': query}, params=params_no_enzymes)
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['primers']), 2)  # 2 primers (forward and reverse)
        fwd_primer = PrimerModel.model_validate(payload['primers'][0])
        rvs_primer = PrimerModel.model_validate(payload['primers'][1])
        self.assertEqual(fwd_primer.name, 'test_sequence_fwd')
        self.assertEqual(rvs_primer.name, 'test_sequence_rvs')
        self.assertTrue(fwd_primer.sequence.startswith('GGCTT'))
        self.assertTrue(rvs_primer.sequence.startswith('CTGTT'))

        # Same, now inverted
        query2 = copy.deepcopy(query)
        query2['forward_orientation'] = False
        response = client.post('/primer_design/simple_pair', json={'pcr_template': query2}, params=params_no_enzymes)
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['primers']), 2)  # 2 primers (forward and reverse)
        fwd_primer = PrimerModel.model_validate(payload['primers'][0])
        rvs_primer = PrimerModel.model_validate(payload['primers'][1])
        self.assertEqual(fwd_primer.name, 'test_sequence_fwd')
        self.assertEqual(rvs_primer.name, 'test_sequence_rvs')
        self.assertTrue(fwd_primer.sequence.startswith('CTGTT'))
        self.assertTrue(rvs_primer.sequence.startswith('GGCTT'))

        # Test primer name when sequence name is unset
        query2 = copy.deepcopy(query)
        dseqr2 = copy.deepcopy(dseqr)
        dseqr2.name = 'name'
        json_seq2 = format_sequence_genbank(dseqr2)
        json_seq2.id = 0
        query2['sequence'] = json_seq2.model_dump()
        response = client.post('/primer_design/simple_pair', json={'pcr_template': query2}, params=params)
        payload = response.json()
        fwd_primer = PrimerModel.model_validate(payload['primers'][0])
        rvs_primer = PrimerModel.model_validate(payload['primers'][1])

        self.assertEqual(fwd_primer.name, 'seq_0_EcoRI_fwd')
        self.assertEqual(rvs_primer.name, 'seq_0_BamHI_rvs')

        # Test with spacers
        response = client.post(
            '/primer_design/simple_pair',
            json={'pcr_template': query, 'spacers': ['ATTT', 'CCAG']},
            params=params,
        )
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['primers']), 2)

        fwd_primer = PrimerModel.model_validate(payload['primers'][0])
        rvs_primer = PrimerModel.model_validate(payload['primers'][1])

        self.assertEqual(fwd_primer.name, 'test_sequence_EcoRI_fwd')
        self.assertEqual(rvs_primer.name, 'test_sequence_BamHI_rvs')

        self.assertTrue(fwd_primer.sequence.startswith('GC' + str(EcoRI.site) + 'ATTT'))
        self.assertTrue(rvs_primer.sequence.startswith('GC' + str(BamHI.site) + 'CTGG'))

        # Test error case with invalid parameters
        params['minimal_hybridization_length'] = 100  # Too long
        response = client.post('/primer_design/simple_pair', json={'pcr_template': query}, params=params)
        self.assertEqual(response.status_code, 400)
        self.assertIn('Primers could not be designed', response.json()['detail'])

        # Test error case with invalid enzyme
        params['minimal_hybridization_length'] = 10  # Reset to valid value
        params['left_enzyme'] = 'InvalidEnzyme'
        response = client.post('/primer_design/simple_pair', json={'pcr_template': query}, params=params)
        self.assertEqual(response.status_code, 404)
        self.assertIn('These enzymes do not exist', response.json()['detail'])

        # Test error case with wrong filler bases
        params['left_enzyme'] = 'EcoRI'
        params['filler_bases'] = 'zAA'
        response = client.post('/primer_design/simple_pair', json={'pcr_template': query}, params=params)
        self.assertEqual(response.status_code, 400)
        self.assertIn('Filler bases can only contain ACGT bases.', response.json()['detail'])

        # Test error case with wrong spacer number
        params['filler_bases'] = 'GC'
        response = client.post(
            '/primer_design/simple_pair', json={'pcr_template': query, 'spacers': ['ATGC']}, params=params
        )
        self.assertEqual(response.status_code, 422)
        self.assertIn('The number of spacers must be', response.json()['detail'])


class ValidateEndPointTest(unittest.TestCase):

    def test_validate(self):
        with open('test_files/homologous_recombination.json') as ins:
            # Read it to json
            data = json.load(ins)
        response = client.post('/validate', json=data)
        self.assertEqual(response.status_code, 200)

        data['dummy'] = 'dummy'
        response = client.post('/validate', json=data)
        self.assertEqual(response.status_code, 422)


class RenameSequenceTest(unittest.TestCase):

    def test_rename(self):
        dseqr = Dseqrecord('ACGT')
        dseqr.name = 'original'
        json_seq = format_sequence_genbank(dseqr)
        json_seq.id = 0
        response = client.post('/rename_sequence?name=hello', json=json_seq.model_dump())
        self.assertEqual(response.status_code, 200)

        payload = response.json()
        dseq_resp = read_dsrecord_from_json(TextFileSequence.model_validate(payload))
        self.assertEqual(dseq_resp.name, 'hello')

    def test_error(self):
        """Does not allow spaces"""
        dseqr = Dseqrecord('ACGT')
        dseqr.name = 'original'
        json_seq = format_sequence_genbank(dseqr)
        json_seq.id = 0
        response = client.post('/rename_sequence?name=hello world', json=json_seq.model_dump())
        self.assertEqual(response.status_code, 422)


class RestrictionEnzymeListTest(unittest.TestCase):

    def test_restriction_enzyme_list(self):
        response = client.get('/restriction_enzyme_list')
        assert response.status_code == 200
        assert 'EcoRI' in response.json()['enzyme_names']


class BenchlingUrlSourceTest(unittest.TestCase):

    def test_valid_url(self):
        url = 'https://benchling.com/siverson/f/lib_B94YxDHhQh-cidar-moclo-library/seq_dh1FrJTc-b0015_dh.gb'
        source = BenchlingUrlSource(id=0, repository_id=url, repository_name='benchling')
        response = client.post('/repository_id/benchling', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['sequences']), 1)
        self.assertEqual(len(payload['sources']), 1)
        self.assertEqual(payload['sources'][0], source.model_dump())

    def test_invalid_url(self):
        # We have to initialize the object with a valid url
        url = 'https://benchling.com/siverson/f/lib_B94YxDHhQh-cidar-moclo-library/seq_dh1FrJTc-b0015_dh.gb'
        source_object = BenchlingUrlSource(id=0, repository_id=url, repository_name='benchling')

        # In the dict, we can then edit
        source_dict = source_object.model_dump()

        # url missing the .gb
        source_dict['repository_id'] = (
            'https://benchling.com/siverson/f/lib_B94YxDHhQh-cidar-moclo-library/seq_dh1FrJTc-b0015_dh'
        )
        response = client.post('/repository_id/benchling', json=source_dict)
        self.assertEqual(response.status_code, 422)

        # The /edit url
        source_dict['repository_id'] = (
            'https://benchling.com/siverson/f/lib_B94YxDHhQh-cidar-moclo-library/seq_dh1FrJTc-b0015_dh/edit'
        )
        response = client.post('/repository_id/benchling', json=source_dict)
        self.assertEqual(response.status_code, 422)

        # One that matches the pattern but does not exist
        url = 'https://benchling.com/bluh/blah.gb'
        source = BenchlingUrlSource(id=0, repository_id=url, repository_name='benchling')
        response = client.post('/repository_id/benchling', json=source.model_dump())
        self.assertEqual(response.status_code, 404)
        self.assertIn('file requested from url not found', response.json()['detail'])


class SnapGenePlasmidSourceTest(unittest.TestCase):

    def test_valid_url(self):

        source = SnapGenePlasmidSource(
            id=0, repository_id='basic_cloning_vectors/pEASY-T1_(linearized)', repository_name='snapgene'
        )
        response = client.post('/repository_id/snapgene', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['sequences']), 1)
        self.assertEqual(len(payload['sources']), 1)
        out_source = payload['sources'][0]
        # It should have been renamed
        self.assertEqual(out_source['output_name'], 'pEASY-T1_(linearized)')
        # Remove that and compare with original source
        out_source['output_name'] = None
        self.assertEqual(out_source, source.model_dump())
        seq = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertEqual(seq.name, 'pEASY-T1_(linearized)')

        # We can also provide a name
        source2 = SnapGenePlasmidSource(
            id=0,
            repository_id='basic_cloning_vectors/pEASY-T1_(linearized)',
            repository_name='snapgene',
            output_name='my_name',
        )
        response = client.post('/repository_id/snapgene', json=source2.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(payload['sources'][0]['output_name'], 'my_name')
        seq = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertEqual(seq.name, 'my_name')

    def test_invalid_url(self):
        # Invalid url
        source = SnapGenePlasmidSource(id=0, repository_id='hello/world', repository_name='snapgene')
        response = client.post('/repository_id/snapgene', json=source.model_dump())
        self.assertEqual(response.status_code, 404)

        # Wrongly formatted url
        source_dict = source.model_dump()
        source_dict['repository_id'] = 'hello'
        response = client.post('/repository_id/snapgene', json=source_dict)
        self.assertEqual(response.status_code, 422)


class EuroscarfSourceTest(unittest.TestCase):

    def test_valid_url(self):
        source = EuroscarfSource(id=0, repository_id='P30174', repository_name='euroscarf')
        response = client.post('/repository_id/euroscarf', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['sequences']), 1)
        self.assertEqual(len(payload['sources']), 1)
        out_source = payload['sources'][0]
        self.assertEqual(out_source, source.model_dump())
        sequence = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertEqual(sequence.name, 'pKT128')
        self.assertEqual(len(sequence), 4738)
        self.assertTrue(any('yEGFP' in f.qualifiers['gene'] for f in sequence.features))

        # Ensure that linear files are circularised
        source = EuroscarfSource(id=0, repository_id='P30555', repository_name='euroscarf')
        response = client.post('/repository_id/euroscarf', json=source.model_dump())
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['sequences']), 1)
        sequence = read_dsrecord_from_json(TextFileSequence.model_validate(payload['sequences'][0]))
        self.assertTrue(sequence.circular)

    def test_invalid_url(self):
        # Compatible with regex, but does not exist
        source = EuroscarfSource(id=0, repository_id='P99999999999999', repository_name='euroscarf')
        response = client.post('/repository_id/euroscarf', json=source.model_dump())
        self.assertEqual(response.status_code, 404)

        # Not compatible with regex
        source_dict = source.model_dump()
        source_dict['repository_id'] = 'hello'
        response = client.post('/repository_id/euroscarf', json=source_dict)
        self.assertEqual(response.status_code, 422)


if __name__ == '__main__':
    unittest.main()
