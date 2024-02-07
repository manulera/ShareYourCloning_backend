
from dna_functions import format_sequence_genbank, read_dsrecord_from_json
from main import app
from fastapi.testclient import TestClient
from pydna.parsers import parse as pydna_parse
from Bio.Restriction.Restriction import CommOnly
from pydantic_models import RepositoryIdSource, PCRSource, PrimerModel, \
    RestrictionEnzymeDigestionSource, SequenceEntity, LigationSource, \
    UploadedFileSource, HomologousRecombinationSource, GibsonAssemblySource, \
    RestrictionAndLigationSource
from pydna.dseqrecord import Dseqrecord
import unittest
from pydna.dseq import Dseq

client = TestClient(app)

# TODO further tests are needed (combinations)


class ReadFileTest(unittest.TestCase):

    def test_read_files(self):
        """Test that uploading files with single and multiple sequences works."""

        examples = [
            {
                'file': './examples/sequences/pFA6a-hphMX6.gb',
                'format': 'genbank',
                'nb_sequences': 1
            },
            {
                'file': './examples/sequences/dummy_EcoRI.fasta',
                'format': 'fasta',
                'nb_sequences': 1
            },
            {
                'file': './examples/sequences/dummy_multi_fasta.fasta',
                'format': 'fasta',
                'nb_sequences': 2
            },
            {
                'file': './examples/sequences/addgene-plasmid-39296-sequence-49545.dna',
                'format': 'snapgene',
                'nb_sequences': 1
            }

        ]

        for example in examples:
            with open(example['file'], 'rb') as f:
                response = client.post("/read_from_file", files={"file": f})

            self.assertEqual(response.status_code, 200)
            payload = response.json()

            resulting_sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
            sources = [UploadedFileSource.model_validate(s) for s in payload['sources']]

            self.assertEqual(len(sources), example['nb_sequences'])
            self.assertEqual(len(resulting_sequences), example['nb_sequences'])
            for seq in resulting_sequences:
                self.assertGreater(len(seq), 3)
            for source in sources:
                self.assertEqual(source.file_format, example['format'])

    def test_errors_read_files(self):

        examples = [
            {
                'file': './test_endpoints.py',
                'format': None,
                'error_message': 'could not guess',
            },
            {
                'file': './test_endpoints.py',
                'format': 'snapgene',
                'error_message': 'snapgene reader cannot',
            },
            {
                'file': './test_endpoints.py',
                'format': 'genbank',
                'error_message': 'Pydna parser',
            }
        ]
        for example in examples:
            with open(example['file'], 'rb') as f:
                if example["format"] is not None:
                    response = client.post(f'/read_from_file?file_format={example["format"]}', files={'file': f})
                else:
                    response = client.post('/read_from_file', files={'file': f})

            self.assertEqual(response.status_code, 422)
            self.assertTrue(example['error_message'] in response.json()['detail'])


class GenBankTest(unittest.TestCase):

    # TODO these tests will not work off-line, so the case where connection cannot be established should be handled in some way
    def test_request_gene(self):
        """Test whether the gene is requested from GenBank"""
        source = RepositoryIdSource(
            repository='genbank',
            repository_id='NM_001018957.2',
        )
        response = client.post("/repository_id", json=source.model_dump())
        self.assertEqual(response.status_code, 200)

    def test_request_wrong_id(self):
        """Test a wrong Genbank id"""
        source = RepositoryIdSource(
            repository='genbank',
            repository_id='wrong_id',
        )
        response = client.post("/repository_id", json=source.model_dump())
        self.assertEqual(response.status_code, 404)


class AddGeneTest(unittest.TestCase):
    def test_request_plasmid(self):
        """Test whether the gene is requested from AddGene and returns the right info"""
        examples = [
            {
                'id': '39282',
                'url': 'https://media.addgene.org/snapgene-media/v1.7.9-0-g88a3305/sequences/240599/4936a6ae-6b4d-4d24-b7ac-2339fad5755d/addgene-plasmid-39282-sequence-240599.gbk',
                'type': 'addgene-full'
            },
            {
                'id': '39289',
                'url': 'https://media.addgene.org/snapgene-media/v1.7.9-0-g88a3305/sequences/49640/4fa9f18b-d5ca-4ac6-a50c-76a6cd15cbab/addgene-plasmid-39289-sequence-49640.gbk',
                'type': 'depositor-full'
            },
        ]
        for example in examples:
            source = RepositoryIdSource(
                repository='addgene',
                repository_id=example['id'],
            )

            response = client.post("/repository_id", json=source.model_dump())
            self.assertEqual(response.status_code, 200)
            payload = response.json()
            resulting_sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
            sources = [RepositoryIdSource.model_validate(s) for s in payload['sources']]

            self.assertEqual(len(resulting_sequences), 1)
            self.assertEqual(len(sources), 1)
            self.assertEqual(sources[0].info['type'], example['type'])
            self.assertEqual(sources[0].info['url'], example['url'])

            # We get the same response when making the response with the url
            response2 = client.post("/repository_id", json=payload['sources'][0])
            self.assertEqual(response.json(), response2.json())

    def test_missing_sequences(self):
        # Non-existing id
        source = RepositoryIdSource(
            repository='addgene',
            repository_id='DUMMYTEST',
        )

        response = client.post("/repository_id", json=source.model_dump())
        self.assertEqual(response.status_code, 404)
        self.assertIn('wrong addgene id', response.json()['detail'])

        # Id that has no full-sequences
        source = RepositoryIdSource(
            repository='addgene',
            repository_id='39291',
        )
        response = client.post("/repository_id", json=source.model_dump())
        self.assertEqual(response.status_code, 404)
        self.assertIn('The requested plasmid does not exist, or does not have', response.json()['detail'])

        # url does not exist
        source = RepositoryIdSource(
            repository='addgene',
            repository_id='39282',
            info={'url': 'https://media.addgene.org/snapgene-media/wrongggggggg.gbk'}
        )
        response = client.post("/repository_id", json=source.model_dump())
        self.assertEqual(response.status_code, 404)

        # TODO url exists but does not match id, or does not match


class StickyLigationTest(unittest.TestCase):

    def test_ligation(self):
        """Test whether assembly is executed when the order is provided"""
        # Load dummy sequence

        initial_sequence = pydna_parse('examples/sequences/dummy_EcoRI.fasta')[0]

        # Restriction cut
        enzyme = CommOnly.format('EcoRI')
        output_list: list[Dseqrecord] = initial_sequence.cut([enzyme])

        # Convert to json to use as input
        json_seqs = [format_sequence_genbank(seq) for seq in output_list]
        json_seqs[0].id = 1
        json_seqs[1].id = 2
        json_seqs = [seq.model_dump() for seq in json_seqs]

        # Assign ids to define deterministic assembly
        source = LigationSource(
            input=[1, 2],
            assembly=[(1, 2, '8..11', '1..4')],
            circular=False
        )
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post("/ligation", json=data)
        payload = response.json()

        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
        sources = [LigationSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)

        # Check that the assembly is correct
        self.assertEqual(resulting_sequences[0].seq, initial_sequence.seq)
        self.assertEqual(sources[0], source)

        # Check that the inverse assembly will not pass
        source = LigationSource(
            input=[2, 1],
            assembly=[(1, 2, '8..11', '1..4')],
            circular=False
        )
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post("/ligation", json=data)
        self.assertEqual(response.status_code, 400)
        data = response.json()
        self.assertEqual(data['detail'], 'The provided assembly is not valid.')

        # Check that the circular assembly does not pass either
        source = LigationSource(
            input=[1, 2],
            assembly=[(1, 2, '8..11', '1..4')],
            circular=True
        )
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post("/ligation", json=data)
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
        source = LigationSource(
            input=[1, 2],
        )
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post("/ligation", json=data)
        payload = response.json()

        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
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

        source = LigationSource(
            input=[1, 2],
        )
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post("/ligation", json=data)
        payload = response.json()
        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
        sources = [LigationSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)

        # Check that the assembly is correct, the sequences cannot be compared with the equal operator since their origin is different
        # TODO there should be something in pydna to test circular molecules being equal but having different frames

        self.assertEqual(len(resulting_sequences[0]), len(initial_sequence))

    def test_circularisation(self):
        enzyme = CommOnly.format('EcoRI')
        fragment = Dseqrecord('AGAATTC', circular=True).cut(enzyme)[0]
        json_seqs = [format_sequence_genbank(fragment)]
        json_seqs[0].id = 1
        json_seqs = [seq.model_dump() for seq in json_seqs]
        source = LigationSource(
            input=[1],
        )
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post("/ligation", json=data)
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
        self.assertEqual(resulting_sequences[0].seq, fragment.looped().seq)


class BluntLigationTest(unittest.TestCase):

    def test_blunt_ligation(self):
        seqs =[
            Dseqrecord('ATCC', circular=False),
            Dseqrecord('TAAT', circular=False)
        ]
        json_seqs = [format_sequence_genbank(seq) for seq in seqs]
        json_seqs[0].id = 1
        json_seqs[1].id = 2
        json_seqs = [seq.model_dump() for seq in json_seqs]
        source = LigationSource(
            input=[1, 2],
        )
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post("/ligation", json=data, params={'blunt': True})
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
        self.assertEqual(len(resulting_sequences), 2)

        # We submit one of the resulting sources, to check that it does the
        # blunt ligation without adding the request parameter
        source = payload['sources'][0]
        response = client.post("/ligation", json={'source': source, 'sequences': json_seqs})
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        resulting_sequences2 = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
        self.assertEqual(len(resulting_sequences2), 1)
        self.assertEqual(resulting_sequences2[0].seq, resulting_sequences[0].seq)

    def test_circularisation(self):
        seqs =[
            Dseqrecord('ATCC', circular=False)
        ]
        json_seqs = [format_sequence_genbank(seq) for seq in seqs]
        json_seqs[0].id = 1
        json_seqs = [seq.model_dump() for seq in json_seqs]
        source = LigationSource(
            input=[1],
        )
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post("/ligation", json=data, params={'blunt': True})
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(resulting_sequences[0].seq, seqs[0].looped().seq)


class RestrictionTest(unittest.TestCase):

    def test_enzyme_doesnt_exist(self):
        dseq = Dseqrecord('AAAAAAGAATTCTTTTTT', circular=False)
        json_seq = format_sequence_genbank(dseq)
        json_seq.id = 1

        # One enzyme
        source = RestrictionEnzymeDigestionSource(
            input=[1],
            restriction_enzymes=['helloworld'],
        )
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data)
        self.assertEqual(response.status_code, 404)
        self.assertTrue(response.json()['detail'] == 'These enzymes do not exist: helloworld')

        # More than one
        source = RestrictionEnzymeDigestionSource(
            input=[1],
            restriction_enzymes=['helloworld', 'byebye'],
        )
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data)
        self.assertEqual(response.status_code, 404)
        self.assertIn('byebye', response.json()['detail'])
        self.assertIn('helloworld', response.json()['detail'])

    def test_enzymes_dont_cut(self):
        dseq = Dseqrecord('AAAAAAGAATTCTTTTTT', circular=False)
        json_seq = format_sequence_genbank(dseq)
        json_seq.id = 1

        # See if we get the right responses
        source = RestrictionEnzymeDigestionSource(
            input=[1],
            restriction_enzymes=['FbaI'],
        )
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data)
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json()['detail'], 'These enzymes do not cut: FbaI')

        # Returns when one does not cut
        source = RestrictionEnzymeDigestionSource(
            input=[1],
            restriction_enzymes=['EcoRI', 'BamHI'],
        )
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data)
        self.assertEqual(response.status_code, 400)
        self.assertTrue(response.json()['detail'] == 'These enzymes do not cut: BamHI')

        # Returns all that don't cut
        source = RestrictionEnzymeDigestionSource(
            input=[1],
            restriction_enzymes=['EcoRI', 'BamHI', 'FbaI'],
        )
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data)
        self.assertEqual(response.status_code, 400)
        self.assertIn('BamHI', response.json()['detail'])
        self.assertIn('FbaI', response.json()['detail'])

    def test_linear_single_restriction(self):

        dseq = Dseqrecord('AAAAAAGAATTCTTTTTT', circular=False)
        json_seq = format_sequence_genbank(dseq)
        json_seq.id = 1

        # See if we get the right responses
        source = RestrictionEnzymeDigestionSource(
            input=[1],
            restriction_enzymes=['EcoRI'],
        )
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data)
        payload = response.json()

        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
        sources = [RestrictionEnzymeDigestionSource.model_validate(s) for s in payload['sources']]

        # The right sequences are returned in the right order
        self.assertEqual(resulting_sequences[0].seq.watson.upper(), 'AAAAAAG')
        self.assertEqual(resulting_sequences[1].seq.watson.upper(), 'AATTCTTTTTT')

        # The edges of the fragments are correct:
        self.assertEqual(sources[0].left_edge, None)
        self.assertEqual(sources[0].right_edge, (7, -4))

        self.assertEqual(sources[1].left_edge, (7, -4))
        self.assertEqual(sources[1].right_edge, None)

        # The enzyme names are correctly returned:
        self.assertEqual(sources[0].restriction_enzymes, [None, 'EcoRI'])
        self.assertEqual(sources[1].restriction_enzymes, ['EcoRI', None])

        # Now we specify the output
        source = RestrictionEnzymeDigestionSource(
            input=[1],
            restriction_enzymes=[None, 'EcoRI'],
            left_edge=None,
            right_edge=(7, -4)
        )
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data)
        payload = response.json()

        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
        sources = [RestrictionEnzymeDigestionSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)
        self.assertEqual(sources[0].left_edge, None)
        self.assertEqual(sources[0].right_edge, (7, -4))
        self.assertEqual(sources[0].restriction_enzymes, [None, 'EcoRI'])

    def test_circular_single_restriction(self):

        dseq = Dseqrecord('AAAAAAGAATTCTTTTTT', circular=True)
        json_seq = format_sequence_genbank(dseq)
        json_seq.id = 1

        # See if we get the right responses
        source = RestrictionEnzymeDigestionSource(
            input=[1],
            restriction_enzymes=['EcoRI'],
        )
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data)
        payload = response.json()

        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
        sources = [RestrictionEnzymeDigestionSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)

        # The right sequences are returned in the right order
        self.assertEqual(resulting_sequences[0].seq.watson.upper(), 'AATTCTTTTTTAAAAAAG')

        # The edges of the fragments are correct:
        self.assertEqual(sources[0].left_edge, (7, -4))
        self.assertEqual(sources[0].right_edge, (7, -4))

        # The enzyme names are correctly returned:
        self.assertEqual(sources[0].restriction_enzymes, ['EcoRI', 'EcoRI'])

        # When the cutting site spans the origin
        sequences = ['AATTCTTTTTTG', 'ATTCTTTTTTGA']
        cut_positions = [(0, -4), (11, -4)]

        for s, pos in zip(sequences, cut_positions):
            dseq = Dseqrecord(s, circular=True)
            json_seq = format_sequence_genbank(dseq)
            json_seq.id = 1

            # See if we get the right responses
            source = RestrictionEnzymeDigestionSource(
                input=[1],
                restriction_enzymes=['EcoRI'],
            )
            data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
            response = client.post('/restriction', json=data)
            payload = response.json()

            resulting_sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
            sources = [RestrictionEnzymeDigestionSource.model_validate(s) for s in payload['sources']]

            self.assertEqual(len(resulting_sequences), 1)
            self.assertEqual(len(sources), 1)

            # The right sequences are returned in the right order
            self.assertEqual(resulting_sequences[0].seq.watson.upper(), 'AATTCTTTTTTG')

            # The edges of the fragments are correct:
            self.assertEqual(sources[0].left_edge, pos)
            self.assertEqual(sources[0].right_edge, pos)

            # The enzyme names are correctly returned:
            self.assertEqual(sources[0].restriction_enzymes, ['EcoRI', 'EcoRI'])

    def test_linear_multiple_restriction(self):

        dseq = Dseqrecord('AAAGGATCCAAAAGATATCAAAAA', circular=False)
        json_seq = format_sequence_genbank(dseq)
        json_seq.id = 1

        # We try EcoRV, which gives blunt ends
        source = RestrictionEnzymeDigestionSource(
            input=[1],
            restriction_enzymes=['BamHI', 'EcoRV'],
        )
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data)
        payload = response.json()

        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
        sources = [RestrictionEnzymeDigestionSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 3)
        self.assertEqual(len(sources), 3)

        # The right sequences are returned in the right order
        fragment_sequences = 'AAAG^GATCCAAAAGAT^ATCAAAAA'.split('^')
        for i, s in enumerate(fragment_sequences):
            self.assertEqual(resulting_sequences[i].seq.watson.upper(), s)

        # The edges of the fragments are correct:
        # AAAG^GATC_CAAAAGAT^ATCAAAAA
        edges = [(None, (4, -4)), ((4, -4), (16, 0)), ((16, 0), None)]
        for i, e in enumerate(edges):
            self.assertEqual(sources[i].left_edge, e[0])
            self.assertEqual(sources[i].right_edge, e[1])

        # The enzyme names are correct
        restriction_enzymes = [[None, 'BamHI'], ['BamHI', 'EcoRV'], ['EcoRV', None]]
        for i, e in enumerate(restriction_enzymes):
            self.assertEqual(sources[i].restriction_enzymes, e)

        # Submitting the known fragments

        for i in range(len(edges)):
            source = RestrictionEnzymeDigestionSource(
                input=[1],
                restriction_enzymes=restriction_enzymes[i],
                left_edge=edges[i][0],
                right_edge=edges[i][1]
            )
            data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
            response = client.post('/restriction', json=data)
            payload = response.json()

            resulting_sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
            sources = [RestrictionEnzymeDigestionSource.model_validate(s) for s in payload['sources']]

            self.assertEqual(len(resulting_sequences), 1)
            self.assertEqual(len(sources), 1)
            self.assertEqual(resulting_sequences[0].seq.watson, fragment_sequences[i])

    def test_circular_multiple_restriction(self):

        dseq = Dseqrecord('AAAGGATCCAAAAGATATCAAAAA', circular=True)
        json_seq = format_sequence_genbank(dseq)
        json_seq.id = 1

        # We try EcoRV, which gives blunt ends
        source = RestrictionEnzymeDigestionSource(
            input=[1],
            restriction_enzymes=['BamHI', 'EcoRV'],
        )
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data)
        payload = response.json()

        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
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
        edges = [((4, -4), (16, 0)), ((16, 0), (4, -4))]
        for i, e in enumerate(edges):
            self.assertEqual(sources[i].left_edge, e[0])
            self.assertEqual(sources[i].right_edge, e[1])

        # The enzyme names are correct
        restriction_enzymes = [['BamHI', 'EcoRV'], ['EcoRV', 'BamHI']]
        for i, e in enumerate(restriction_enzymes):
            self.assertEqual(sources[i].restriction_enzymes, e)

        # Submitting the known fragments
        for i in range(len(edges)):
            source = RestrictionEnzymeDigestionSource(
                input=[1],
                restriction_enzymes=restriction_enzymes[i],
                left_edge=edges[i][0],
                right_edge=edges[i][1]
            )
            data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
            response = client.post('/restriction', json=data)
            payload = response.json()

            resulting_sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
            sources = [RestrictionEnzymeDigestionSource.model_validate(s) for s in payload['sources']]

            self.assertEqual(len(resulting_sequences), 1)
            self.assertEqual(len(sources), 1)
            self.assertEqual(resulting_sequences[0].seq.watson, fragment_sequences[i])


class PCRTest(unittest.TestCase):

    def test_pcr(self):

        template = Dseqrecord(Dseq('TTTTACGTACGTAAAAAAGCGCGCGCTTTTT'))

        json_seq = format_sequence_genbank(template)
        json_seq.id = 1

        submitted_source = PCRSource(
            input=[1],
            forward_primer=2,
            reverse_primer=3,
        )

        primer_fwd = PrimerModel(
            sequence='ACGTACGT',
            id=2,
            name='forward'
        )

        primer_rvs = PrimerModel(
            sequence='GCGCGCGC',
            id=3,
            name='reverse'
        )

        data = {'source': submitted_source.model_dump(), 'sequences': [json_seq.model_dump()], 'primers': [primer_fwd.model_dump(), primer_rvs.model_dump()]}
        response = client.post("/pcr", json=data, params={'minimal_annealing': 8})
        payload = response.json()

        sources = [PCRSource.model_validate(s) for s in payload['sources']]
        sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]

        # Only one possible output
        self.assertEqual(len(sources), 1)
        self.assertEqual(len(sequences), 1)
        dseq1 = sequences[0]
        source1 = sources[0]

        # The sequence matches what we expect
        self.assertEqual(str(dseq1.seq), 'ACGTACGTAAAAAAGCGCGCGC')
        self.assertEqual(source1.forward_primer, primer_fwd.id)
        self.assertEqual(source1.reverse_primer, primer_rvs.id)

        # Now we submit the deterministic PCR (we already know which fragment we want)

        data = {'source': source1.model_dump(), 'sequences': [json_seq.model_dump()], 'primers': [primer_fwd.model_dump(), primer_rvs.model_dump()]}
        response = client.post("/pcr", json=data)
        payload = response.json()

        sources = [PCRSource.model_validate(s) for s in payload['sources']]
        sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]

        # Only one possible output
        self.assertEqual(len(sources), 1)
        self.assertEqual(len(sequences), 1)
        dseq2 = sequences[0]
        source2 = sources[0]
        self.assertEqual(source1, source2)
        self.assertEqual(str(dseq2.seq), 'ACGTACGTAAAAAAGCGCGCGC')

    def test_wrong_primers(self):

        template = Dseqrecord(Dseq('TTTTACGTACGTAAAAAAGCGCGCGCTTTTT'))
        json_seq = format_sequence_genbank(template)
        json_seq.id = 1

        submitted_source = PCRSource(input=[1], forward_primer=2, reverse_primer=3)

        primer_fwd = PrimerModel(
            sequence='CCCCCCCC',
            id=2,
            name='forward'

        )

        primer_rvs = PrimerModel(
            sequence='GCGCGCGC',
            id=3,
            name='reverse'
        )

        # Without specifying the pair
        data = {'source': submitted_source.model_dump(), 'sequences': [json_seq.model_dump()], 'primers': [primer_fwd.model_dump(), primer_rvs.model_dump()]}
        response = client.post("/pcr", json=data, params={'minimal_annealing': 8})
        payload = response.json()
        self.assertEqual(response.status_code, 400)
        self.assertEqual(payload['detail'], 'No pair of annealing primers was found. Try changing the annealing settings.')

        # We submit the right pair of primers, but the wrong annealing information

        primer_fwd = PrimerModel(
            sequence='ACGTACGT',
            id=2,
            name='forward'
        )

        # This would be the correct annealing info
        submitted_source.forward_primer = 2
        submitted_source.reverse_primer = 3
        submitted_source.assembly = [(1, 2, '1..8', '5..12'), (2, -3, '19..26', '1..8')]

        data = {'source': submitted_source.model_dump(), 'sequences': [json_seq.model_dump()], 'primers': [primer_fwd.model_dump(), primer_rvs.model_dump()]}
        response = client.post("/pcr", json=data)
        payload = response.json()
        self.assertEqual(response.status_code, 200)

        # This is the wrong annealing info
        submitted_source.forward_primer = 2
        submitted_source.reverse_primer = 3
        submitted_source.assembly = [(2, -3, '19..26', '1..8'), (1, 2, '1..8', '5..12')]

        data = {'source': submitted_source.model_dump(), 'sequences': [json_seq.model_dump()], 'primers': [primer_fwd.model_dump(), primer_rvs.model_dump()]}
        response = client.post("/pcr", json=data)
        payload = response.json()
        self.assertEqual(response.status_code, 400)
        self.assertEqual(payload['detail'], 'The provided assembly is not valid.')

        # Test clashing primers
        template = Dseqrecord(Dseq('ACGTACGTGCGCGCGC'))

        json_seq = format_sequence_genbank(template)
        json_seq.id = 1

        submitted_source = PCRSource(
            input=[1],
            forward_primer=2,
            reverse_primer=3,
        )

        primer_fwd = PrimerModel(
            sequence='ACGTACGTG',
            id=2,
            name='forward'
        )

        primer_rvs = PrimerModel(
            sequence='GCGCGCGCA',
            id=3,
            name='reverse'
        )

        data = {'source': submitted_source.model_dump(), 'sequences': [json_seq.model_dump()], 'primers': [primer_fwd.model_dump(), primer_rvs.model_dump()]}
        response = client.post("/pcr", json=data, params={'minimal_annealing': 8})
        payload = response.json()
        self.assertEqual(response.status_code, 400)
        self.assertIn('Clashing primers', payload['detail'])


class HomologousRecombinationTest(unittest.TestCase):

    def test_homologous_recombination(self):
        template = Dseqrecord('TTTTacgatAAtgctccCCCC', circular=False)
        json_template = format_sequence_genbank(template)
        json_template.id = 1

        insert = Dseqrecord('acgatCCCtgctcc', circular=False)
        json_insert = format_sequence_genbank(insert)
        json_insert.id = 2
        # One enzyme
        source = HomologousRecombinationSource(
            input=[1, 2],
        )
        data = {'source': source.model_dump(), 'sequences': [json_template.model_dump(), json_insert.model_dump()]}
        response = client.post('/homologous_recombination', params={'minimal_homology': 5}, json=data)
        self.assertEqual(response.status_code, 200)

        payload = response.json()
        sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
        self.assertEqual(len(sequences), 1)
        self.assertEqual(str(sequences[0].seq), 'TTTTacgatCCCtgctccCCCC'.upper())

        response = client.post('/homologous_recombination', json={'source': payload['sources'][0], 'sequences': [json_template.model_dump(), json_insert.model_dump()]})
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
        self.assertEqual(len(sequences), 1)
        self.assertEqual(str(sequences[0].seq), 'TTTTacgatCCCtgctccCCCC'.upper())


class GibsonAssemblyTest(unittest.TestCase):

    def test_gibson_assembly(self):
        # A circular one
        fragments = [
            Dseqrecord('TTTTacgatAAtgctccCCCC', circular=False),
            Dseqrecord('CCCCtcatGGGG', circular=False),
            Dseqrecord('GGGGatataTTTT', circular=False)
        ]

        json_fragments = [format_sequence_genbank(f) for f in fragments]
        for i, f in enumerate(json_fragments):
            f.id = i + 1

        source = GibsonAssemblySource(
            input=[1, 2, 3],
        )

        data = {'source': source.model_dump(), 'sequences': [f.model_dump() for f in json_fragments]}
        response = client.post('/gibson_assembly', json=data, params={'minimal_homology': 4})
        self.assertEqual(response.status_code, 200)
        payload = response.json()
        sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]

        self.assertEqual(len(sequences), 2)
        self.assertEqual(str(sequences[0].seq), 'TTTTacgatAAtgctccCCCCtcatGGGGatata'.upper())
        self.assertEqual(str(sequences[1].seq), 'TTTTacgatAAtgctccCCCCatgaGGGGatata'.upper())

        # Circularisation works
        f1 = Dseqrecord('AGAGACCaaaAGAGACC')
        json_fragments = [format_sequence_genbank(f1)]
        for i, f in enumerate(json_fragments):
            f.id = i + 1

        source = GibsonAssemblySource(
            input=[1],
        )
        data = {'source': source.model_dump(), 'sequences': [f.model_dump() for f in json_fragments]}
        response = client.post('/gibson_assembly', json=data, params={'minimal_homology': 7})

        self.assertEqual(response.status_code, 200)
        payload = response.json()
        sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]

        self.assertEqual(len(sequences), 1)
        self.assertEqual(str(sequences[0].seq), 'AGAGACCaaa'.upper())


class RestrictionAndLigationTest(unittest.TestCase):

    def test_restriction_and_ligation(self):
        fragments = [Dseqrecord('AAAGAATTCAAA'), Dseqrecord('CCCCGAATTCCCC')]
        [format_sequence_genbank(f) for f in fragments]
        json_fragments = [format_sequence_genbank(f) for f in fragments]
        for i, f in enumerate(json_fragments):
            f.id = i + 1

        source = RestrictionAndLigationSource(
            input=[1, 2, 3],
            restriction_enzymes=['EcoRI'],
        )

        data = {'source': source.model_dump(), 'sequences': [f.model_dump() for f in json_fragments]}
        response = client.post('/restriction_and_ligation', json=data, params={'minimal_homology': 4})

        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['sequences']), 4)
        self.assertEqual(len(payload['sources']), 4)

    def test_golden_gate(self):
        fragments =[ Dseqrecord('GGTCTCAattaAAAAAttaaAGAGACC'),
            Dseqrecord('GGTCTCAttaaCCCCCatatAGAGACC'),
            Dseqrecord('GGTCTCAatatGGGGGccggAGAGACC'),
            Dseqrecord('TTTTattaAGAGACCTTTTTGGTCTCAccggTTTT', circular=True),
        ]

        json_fragments = [format_sequence_genbank(f) for f in fragments]
        for i, f in enumerate(json_fragments):
            f.id = i + 1

        source = RestrictionAndLigationSource(
            input=[1, 2, 3, 4],
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
            input=[1],
            restriction_enzymes=['EcoRI'],
        )

        data = {'source': source.model_dump(), 'sequences': [f.model_dump() for f in json_fragments]}
        response = client.post('/restriction_and_ligation', json=data)

        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(len(payload['sequences']), 2)
        self.assertEqual(len(payload['sources']), 2)

        sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]

        self.assertEqual(str(sequences[0].seq), 'AATTCAAAG')
        self.assertEqual(str(sequences[1].seq), 'AAAGAATTCAAAA')



if __name__ == "__main__":
    unittest.main()
