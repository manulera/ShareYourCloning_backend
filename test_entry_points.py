
from dna_functions import format_sequence_genbank, read_dsrecord_from_json
from main import app
from fastapi.testclient import TestClient
from pydna.parsers import parse as pydna_parse
from Bio.Restriction.Restriction import CommOnly
from pydantic_models import RepositoryIdSource, PCRSource, PrimerModel, \
    RestrictionEnzymeDigestionSource, SequenceEntity, StickyLigationSource, UploadedFileSource, HomologousRecombinationSource
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
                'file': './test_entry_points.py',
                'format': None,
                'error_message': 'could not guess',
            },
            {
                'file': './test_entry_points.py',
                'format': 'snapgene',
                'error_message': 'snapgene reader cannot',
            },
            {
                'file': './test_entry_points.py',
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

    def test_sticky_ligation(self):
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
        source = StickyLigationSource(
            input=[1, 2],
            assembly=[(1, 2, '8..11', '1..4')],
            circular=False
        )
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post("/sticky_ligation", json=data)
        payload = response.json()

        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
        sources = [StickyLigationSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)

        # Check that the assembly is correct
        self.assertEqual(resulting_sequences[0].seq, initial_sequence.seq)
        self.assertEqual(sources[0], source)

        # Check that the inverse assembly will not pass
        source = StickyLigationSource(
            input=[2, 1],
            assembly=[(1, 2, '8..11', '1..4')],
            circular=False
        )
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post("/sticky_ligation", json=data)
        self.assertEqual(response.status_code, 400)
        data = response.json()
        self.assertEqual(data['detail'], 'The provided assembly is not valid.')

        # Check that the circular assembly does not pass either
        source = StickyLigationSource(
            input=[1, 2],
            assembly=[(1, 2, '8..11', '1..4')],
            circular=True
        )
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post("/sticky_ligation", json=data)
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
        source = StickyLigationSource(
            input=[1, 2],
        )
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post("/sticky_ligation", json=data)
        payload = response.json()

        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
        sources = [StickyLigationSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)

        # Check that the assembly is correct
        self.assertEqual(resulting_sequences[0].seq, initial_sequence.seq)

    def test_circular_assembly_no_order(self):
        initial_sequence = pydna_parse('examples/sequences/pFA6a-hphMX6.gb')[0]

        # Restriction cut
        output_list: list[Dseqrecord] = initial_sequence.cut([CommOnly.format('AscI'), CommOnly.format('SacI')])
        for seq in output_list:
            print(seq.seq.__repr__())

        # Convert to json to use as input
        json_seqs = [format_sequence_genbank(seq) for seq in output_list]
        json_seqs[0].id = 1
        json_seqs[1].id = 2
        json_seqs = [seq.model_dump() for seq in json_seqs]

        # We don't set the fragments_inverted, so we will get all possibilities (in this case only one)
        source = StickyLigationSource(
            input=[1, 2],
        )
        data = {'source': source.model_dump(), 'sequences': json_seqs}
        response = client.post("/sticky_ligation", json=data)
        payload = response.json()

        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
        sources = [StickyLigationSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)

        # Check that the assembly is correct, the sequences cannot be compared with the equal operator since their origin is different
        # TODO there should be something in pydna to test circular molecules being equal but having different frames

        self.assertEqual(len(resulting_sequences[0]), len(initial_sequence))


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
        self.assertEqual(sources[0].right_edge, (7, 11))

        self.assertEqual(sources[1].left_edge, (7, 11))
        self.assertEqual(sources[1].right_edge, None)

        # The enzyme names are correctly returned:
        self.assertEqual(sources[0].restriction_enzymes, [None, 'EcoRI'])
        self.assertEqual(sources[1].restriction_enzymes, ['EcoRI', None])

        # Now we specify the output
        source = RestrictionEnzymeDigestionSource(
            input=[1],
            restriction_enzymes=[None, 'EcoRI'],
            left_edge=None,
            right_edge=(7, 11)
        )
        data = {'source': source.model_dump(), 'sequences': [json_seq.model_dump()]}
        response = client.post('/restriction', json=data)
        payload = response.json()

        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]
        sources = [RestrictionEnzymeDigestionSource.model_validate(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)
        self.assertEqual(sources[0].left_edge, None)
        self.assertEqual(sources[0].right_edge, (7, 11))
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
        self.assertEqual(sources[0].left_edge, (7, 11))
        self.assertEqual(sources[0].right_edge, (7, 11))

        # The enzyme names are correctly returned:
        self.assertEqual(sources[0].restriction_enzymes, ['EcoRI', 'EcoRI'])

        # When the cutting site spans the origin
        sequences = ['AATTCTTTTTTG', 'ATTCTTTTTTGA']
        cut_positions = [(0, 4), (11, 3)]

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
        edges = [(None, (4, 8)), ((4, 8), (16, 16)), ((16, 16), None)]
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
        edges = [((4, 8), (16, 16)), ((16, 16), (4, 8))]
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
        template = pydna_parse('examples/sequences/pFA6a-hphMX6.gb')[0]
        json_seq = format_sequence_genbank(template)
        json_seq.id = 1

        submitted_source = PCRSource(
            input=[1]
        )

        primer_fwd = PrimerModel(
            sequence='AGTTTTCATATCTTCCTTTATATTCTATTAATTGAATTTCAAACATCGTTTTATTGAGCTCATTTACATCAACCGGTTCACGGATCCCCGGGTTAATTAA',
            id=2,
            name='ase1_forward'
        )

        primer_rvs = PrimerModel(
            sequence='CTTTTATGAATTATCTATATGCTGTATTCATATGCAAAAATATGTATATTTAAATTTGATCGATTAGGTAAATAAGAAGCGAATTCGAGCTCGTTTAAAC',
            id=3,
            name='ase1_reverse'
        )

        data = {'source': submitted_source.model_dump(), 'sequences': [json_seq.model_dump()], 'primers': [primer_fwd.model_dump(), primer_rvs.model_dump()]}
        response = client.post("/pcr", json=data, params={'minimal_annealing': 13})
        payload = response.json()

        sources = [PCRSource.model_validate(s) for s in payload['sources']]
        sequences = [read_dsrecord_from_json(SequenceEntity.model_validate(s)) for s in payload['sequences']]

        # Only one possible output
        self.assertEqual(len(sources), 1)
        self.assertEqual(len(sequences), 1)
        dseq1 = sequences[0]
        source1 = sources[0]

        # The sequence matches what we expect
        predicted_seq = primer_fwd.sequence + template[source1.fragment_boundaries[0]:source1.fragment_boundaries[1]] + Dseq(primer_rvs.sequence).reverse_complement()
        self.assertEqual(dseq1.seq, predicted_seq.seq)
        self.assertEqual(source1.primers[0], primer_fwd.id)
        self.assertEqual(source1.primers[1], primer_rvs.id)

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
        self.assertEqual(dseq2.seq, predicted_seq.seq)

    def test_wrong_primers(self):

        template = pydna_parse('examples/sequences/pFA6a-hphMX6.gb')[0]
        json_seq = format_sequence_genbank(template)
        json_seq.id = 1

        submitted_source = PCRSource(input=[1])

        primer_fwd = PrimerModel(
            sequence='AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
            id=2,
            name='ase1_forward'

        )

        primer_rvs = PrimerModel(
            sequence='CTTTTATGAATTATCTATATGCTGTATTCATATGCAAAAATATGTATATTTAAATTTGATCGATTAGGTAAATAAGAAGCGAATTCGAGCTCGTTTAAAC',
            id=3,
            name='ase1_reverse'
        )

        # Without specifying the pair
        data = {'source': submitted_source.model_dump(), 'sequences': [json_seq.model_dump()], 'primers': [primer_fwd.model_dump(), primer_rvs.model_dump()]}
        response = client.post("/pcr", json=data, params={'minimal_annealing': 13})
        payload = response.json()
        self.assertEqual(response.status_code, 400)
        self.assertEqual(payload['detail'], 'No pair of annealing primers was found. Try changing the annealing settings.')

        # We submit the right pair of primers, but the wrong annealing information

        primer_fwd = PrimerModel(
            sequence='AGTTTTCATATCTTCCTTTATATTCTATTAATTGAATTTCAAACATCGTTTTATTGAGCTCATTTACATCAACCGGTTCACGGATCCCCGGGTTAATTAA',
            id=2,
            name='ase1_forward'
        )

        # This would be the correct annealing info
        submitted_source.primers = [2, 3]
        submitted_source.fragment_boundaries = [59, 1718]
        submitted_source.primer_footprints = [21, 20]

        data = {'source': submitted_source.model_dump(), 'sequences': [json_seq.model_dump()], 'primers': [primer_fwd.model_dump(), primer_rvs.model_dump()]}
        response = client.post("/pcr", json=data, params={'minimal_annealing': 13})
        payload = response.json()
        self.assertEqual(response.status_code, 200)

        # This is the wrong annealing info
        submitted_source.primers = [2, 3]
        submitted_source.fragment_boundaries = [59, 1200]  # 1200 instead of 1718
        submitted_source.primer_footprints = [21, 20]

        data = {'source': submitted_source.model_dump(), 'sequences': [json_seq.model_dump()], 'primers': [primer_fwd.model_dump(), primer_rvs.model_dump()]}
        response = client.post("/pcr", json=data, params={'minimal_annealing': 13})
        payload = response.json()
        self.assertEqual(response.status_code, 400)
        self.assertEqual(payload['detail'], 'The annealing positions of the primers seem to be wrong.')

        # This is the wrong annealing info
        submitted_source.primers = [2, 3]
        submitted_source.fragment_boundaries = [59, 1718]

        # Wrong footprint (To return this error, the 'wrong one' cannot mess with the annealing settings
        # e.g. it could not be [21,40] since it would not align with footprint 20)
        submitted_source.primer_footprints = [21, 12]
        data = {'source': submitted_source.model_dump(), 'sequences': [json_seq.model_dump()], 'primers': [primer_fwd.model_dump(), primer_rvs.model_dump()]}
        response = client.post("/pcr", json=data, params={'minimal_annealing': 13})
        payload = response.json()
        self.assertEqual(response.status_code, 400)
        self.assertEqual(payload['detail'], 'The annealing positions of the primers seem to be wrong.')


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


if __name__ == "__main__":
    unittest.main()
