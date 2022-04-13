
from dna_functions import format_sequence_genbank, read_dsrecord_from_json
from main import app
from fastapi.testclient import TestClient
from pydna.parsers import parse as pydna_parse
from Bio.Restriction.Restriction import CommOnly
from pydantic_models import GenbankIdSource, PCRSource, PrimerAnnealingSettings, PrimerModel,\
    RestrictionEnzymeDigestionSource, SequenceEntity, StickyLigationSource, UploadedFileSource
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

            resulting_sequences = [read_dsrecord_from_json(SequenceEntity.parse_obj(s)) for s in payload['sequences']]
            sources = [UploadedFileSource.parse_obj(s) for s in payload['sources']]

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


class GenbankTest(unittest.TestCase):

    def test_request_gene(self):
        """Test whether the gene is requested from GenBank"""
        source = GenbankIdSource(
            genbank_id='NM_001018957.2',
        )
        response = client.post("/genbank_id", json=source.dict())
        self.assertEqual(response.status_code, 200)

    def test_request_wrong_id(self):
        """Test a wrong Genbank id"""
        source = GenbankIdSource(
            genbank_id='wrong_id',
        )
        response = client.post("/genbank_id", json=source.dict())
        self.assertEqual(response.status_code, 404)


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
        json_seqs = [seq.dict() for seq in json_seqs]

        # Assign ids to define deterministic assembly
        source = StickyLigationSource(
            input=[1, 2],
            fragments_inverted=[False, False],
            circularised=False
        )
        data = {'source': source.dict(), 'sequences': json_seqs}
        response = client.post("/sticky_ligation", json=data)
        payload = response.json()

        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.parse_obj(s)) for s in payload['sequences']]
        sources = [StickyLigationSource.parse_obj(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)

        # Check that the assembly is correct
        self.assertEqual(resulting_sequences[0].seq, initial_sequence.seq)
        self.assertEqual(sources[0], source)

        # Check that the inverse assembly will not pass
        source = StickyLigationSource(
            input=[2, 1],
            fragments_inverted=[False, False],
            circularised=False
        )
        data = {'source': source.dict(), 'sequences': json_seqs}
        response = client.post("/sticky_ligation", json=data)
        self.assertEqual(response.status_code, 400)
        data = response.json()
        self.assertEqual(data['detail'], 'Fragments are not compatible for sticky ligation')

        # Check that the circular assembly does not pass either
        source = StickyLigationSource(
            input=[1, 2],
            fragments_inverted=[False, False],
            circularised=True
        )
        data = {'source': source.dict(), 'sequences': json_seqs}
        response = client.post("/sticky_ligation", json=data)
        self.assertEqual(response.status_code, 400)
        data = response.json()
        self.assertEqual(data['detail'], 'Fragments are not compatible for sticky ligation')


class RestrictionTest(unittest.TestCase):

    def test_linear_single_restriction(self):

        dseq = Dseqrecord('AAAAAAGAATTCTTTTTT', linear=True)
        json_seq = format_sequence_genbank(dseq)
        json_seq.id = 1

        # See if we get the right responses
        source = RestrictionEnzymeDigestionSource(
            input=[1],
            restriction_enzymes=['EcoRI'],
        )
        data = {'source': source.dict(), 'sequences': [json_seq.dict()]}
        response = client.post('/restriction', json=data)
        payload = response.json()

        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.parse_obj(s)) for s in payload['sequences']]
        sources = [RestrictionEnzymeDigestionSource.parse_obj(s) for s in payload['sources']]

        # The right sequences are returned in the right order
        self.assertEqual(resulting_sequences[0].seq.watson.upper(), 'AAAAAAG')
        self.assertEqual(resulting_sequences[1].seq.watson.upper(), 'AATTCTTTTTT')

        # The edges of the fragments are correct:
        self.assertEqual(sources[0].fragment_boundaries, [0, 11])
        self.assertEqual(sources[1].fragment_boundaries, [7, 18])

        # The enzyme names are correctly returned:
        self.assertEqual(sources[0].restriction_enzymes, ['', 'EcoRI'])
        self.assertEqual(sources[1].restriction_enzymes, ['EcoRI', ''])

        # Now we specify the output
        source = RestrictionEnzymeDigestionSource(
            input=[1],
            restriction_enzymes=['', 'EcoRI'],
            fragment_boundaries=[0, 11]
        )
        data = {'source': source.dict(), 'sequences': [json_seq.dict()]}
        response = client.post('/restriction', json=data)
        payload = response.json()

        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.parse_obj(s)) for s in payload['sequences']]
        sources = [RestrictionEnzymeDigestionSource.parse_obj(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)
        self.assertEqual(sources[0].fragment_boundaries, [0, 11])
        self.assertEqual(sources[0].restriction_enzymes, ['', 'EcoRI'])

    def test_circular_single_restriction(self):

        dseq = Dseqrecord('AAAAAAGAATTCTTTTTT', linear=False)
        json_seq = format_sequence_genbank(dseq)
        json_seq.id = 1

        # See if we get the right responses
        source = RestrictionEnzymeDigestionSource(
            input=[1],
            restriction_enzymes=['EcoRI'],
        )
        data = {'source': source.dict(), 'sequences': [json_seq.dict()]}
        response = client.post('/restriction', json=data)
        payload = response.json()

        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.parse_obj(s)) for s in payload['sequences']]
        sources = [RestrictionEnzymeDigestionSource.parse_obj(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 1)
        self.assertEqual(len(sources), 1)

        # The right sequences are returned in the right order
        self.assertEqual(resulting_sequences[0].seq.watson.upper(), 'AATTCTTTTTTAAAAAAG')

        # The edges of the fragments are correct:
        self.assertEqual(sources[0].fragment_boundaries, [7, 7 + len(resulting_sequences[0].seq)])

        # The enzyme names are correctly returned:
        self.assertEqual(sources[0].restriction_enzymes, ['EcoRI', 'EcoRI'])

        # When the cutting site spans the origin
        sequences = ['AATTCTTTTTTG', 'ATTCTTTTTTGA']
        cut_positions = [0, 11]

        for s, pos in zip(sequences, cut_positions):
            dseq = Dseqrecord(s, linear=False)
            json_seq = format_sequence_genbank(dseq)
            json_seq.id = 1

            # See if we get the right responses
            source = RestrictionEnzymeDigestionSource(
                input=[1],
                restriction_enzymes=['EcoRI'],
            )
            data = {'source': source.dict(), 'sequences': [json_seq.dict()]}
            response = client.post('/restriction', json=data)
            payload = response.json()

            resulting_sequences = [read_dsrecord_from_json(SequenceEntity.parse_obj(s)) for s in payload['sequences']]
            sources = [RestrictionEnzymeDigestionSource.parse_obj(s) for s in payload['sources']]

            self.assertEqual(len(resulting_sequences), 1)
            self.assertEqual(len(sources), 1)

            # The right sequences are returned in the right order
            self.assertEqual(resulting_sequences[0].seq.watson.upper(), 'AATTCTTTTTTG')

            # The edges of the fragments are correct:
            self.assertEqual(sources[0].fragment_boundaries, [pos, pos + len(resulting_sequences[0].seq)])

            # The enzyme names are correctly returned:
            self.assertEqual(sources[0].restriction_enzymes, ['EcoRI', 'EcoRI'])

    def test_linear_multiple_restriction(self):

        dseq = Dseqrecord('AAAGGATCCAAAAGATATCAAAAA', linear=True)
        json_seq = format_sequence_genbank(dseq)
        json_seq.id = 1

        # We try EcoRV, which gives blunt ends
        source = RestrictionEnzymeDigestionSource(
            input=[1],
            restriction_enzymes=['BamHI', 'EcoRV'],
        )
        data = {'source': source.dict(), 'sequences': [json_seq.dict()]}
        response = client.post('/restriction', json=data)
        payload = response.json()

        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.parse_obj(s)) for s in payload['sequences']]
        sources = [RestrictionEnzymeDigestionSource.parse_obj(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 3)
        self.assertEqual(len(sources), 3)

        # The right sequences are returned in the right order
        fragment_sequences = 'AAAG^GATCCAAAAGAT^ATCAAAAA'.split('^')
        for i, s in enumerate(fragment_sequences):
            self.assertEqual(resulting_sequences[i].seq.watson.upper(), s)

        # The edges of the fragments are correct:
        # AAAG^GATC_CAAAAGAT^ATCAAAAA
        fragment_boundaries = [[0, 8], [4, 16], [16, len(dseq)]]
        for i, e in enumerate(fragment_boundaries):
            self.assertEqual(sources[i].fragment_boundaries, e)

        # The enzyme names are correct
        restriction_enzymes = [['', 'BamHI'], ['BamHI', 'EcoRV'], ['EcoRV', '']]
        for i, e in enumerate(restriction_enzymes):
            self.assertEqual(sources[i].restriction_enzymes, e)

        # Submitting the known fragments

        for i in range(len(fragment_boundaries)):
            source = RestrictionEnzymeDigestionSource(
                input=[1],
                restriction_enzymes=restriction_enzymes[i],
                fragment_boundaries=fragment_boundaries[i],
            )
            data = {'source': source.dict(), 'sequences': [json_seq.dict()]}
            response = client.post('/restriction', json=data)
            payload = response.json()

            resulting_sequences = [read_dsrecord_from_json(SequenceEntity.parse_obj(s)) for s in payload['sequences']]
            sources = [RestrictionEnzymeDigestionSource.parse_obj(s) for s in payload['sources']]

            self.assertEqual(len(resulting_sequences), 1)
            self.assertEqual(len(sources), 1)
            self.assertEqual(resulting_sequences[0].seq.watson, fragment_sequences[i])

    def test_circular_multiple_restriction(self):

        dseq = Dseqrecord('AAAGGATCCAAAAGATATCAAAAA', linear=False)
        json_seq = format_sequence_genbank(dseq)
        json_seq.id = 1

        # We try EcoRV, which gives blunt ends
        source = RestrictionEnzymeDigestionSource(
            input=[1],
            restriction_enzymes=['BamHI', 'EcoRV'],
        )
        data = {'source': source.dict(), 'sequences': [json_seq.dict()]}
        response = client.post('/restriction', json=data)
        payload = response.json()

        resulting_sequences = [read_dsrecord_from_json(SequenceEntity.parse_obj(s)) for s in payload['sequences']]
        sources = [RestrictionEnzymeDigestionSource.parse_obj(s) for s in payload['sources']]

        self.assertEqual(len(resulting_sequences), 2)
        self.assertEqual(len(sources), 2)

        fragment_sequences = 'GATCCAAAAGAT^ATCAAAAAAAAG'.split('^')
        # The right sequences are returned in the right order
        for i, s in enumerate(fragment_sequences):
            self.assertEqual(resulting_sequences[i].seq.watson.upper(), s)

        # The edges of the fragments are correct:
        # AAAG^GATC_CAAAAGAT^ATCAAAAA
        # We use 8 + len because it's an origin-spanning sequence
        fragment_boundaries = [[4, 16], [16, 8 + len(dseq)]]
        for i, e in enumerate(fragment_boundaries):
            self.assertEqual(sources[i].fragment_boundaries, e)

        # The enzyme names are correct
        restriction_enzymes = [['BamHI', 'EcoRV'], ['EcoRV', 'BamHI']]
        for i, e in enumerate(restriction_enzymes):
            self.assertEqual(sources[i].restriction_enzymes, e)

        # Submitting the known fragments
        for i in range(len(fragment_boundaries)):
            source = RestrictionEnzymeDigestionSource(
                input=[1],
                restriction_enzymes=restriction_enzymes[i],
                fragment_boundaries=fragment_boundaries[i],
            )
            data = {'source': source.dict(), 'sequences': [json_seq.dict()]}
            response = client.post('/restriction', json=data)
            payload = response.json()

            resulting_sequences = [read_dsrecord_from_json(SequenceEntity.parse_obj(s)) for s in payload['sequences']]
            sources = [RestrictionEnzymeDigestionSource.parse_obj(s) for s in payload['sources']]

            self.assertEqual(len(resulting_sequences), 1)
            self.assertEqual(len(sources), 1)
            self.assertEqual(resulting_sequences[0].seq.watson, fragment_sequences[i])


class PCRTest(unittest.TestCase):

    def test_pcr(self):
        template = pydna_parse('examples/sequences/pFA6a-hphMX6.gb')[0]
        json_seq = format_sequence_genbank(template)
        json_seq.id = 1

        submitted_source = PCRSource(
            input=[1],
            primer_annealing_settings=PrimerAnnealingSettings(minimum_annealing=13)
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

        data = {'source': submitted_source.dict(), 'sequences': [json_seq.dict()], 'primers': [primer_fwd.dict(), primer_rvs.dict()]}
        response = client.post("/pcr", json=data)
        payload = response.json()

        sources = [PCRSource.parse_obj(s) for s in payload['sources']]
        sequences = [read_dsrecord_from_json(SequenceEntity.parse_obj(s)) for s in payload['sequences']]

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

        data = {'source': source1.dict(), 'sequences': [json_seq.dict()], 'primers': [primer_fwd.dict(), primer_rvs.dict()]}
        response = client.post("/pcr", json=data)
        payload = response.json()

        sources = [PCRSource.parse_obj(s) for s in payload['sources']]
        sequences = [read_dsrecord_from_json(SequenceEntity.parse_obj(s)) for s in payload['sequences']]

        # Only one possible output
        self.assertEqual(len(sources), 1)
        self.assertEqual(len(sequences), 1)
        dseq2 = sequences[0]
        source2 = sources[0]

        # When we submit the information, the annealing_settings are changed
        self.assertEqual(source2.primer_annealing_settings.minimum_annealing, min(source1.primer_footprints))

        # To compare, we unset the annealing settings (see PCRSource)
        source1.primer_annealing_settings = None
        source2.primer_annealing_settings = None
        self.assertEqual(source1, source2)
        self.assertEqual(dseq2.seq, predicted_seq.seq)

    def test_wrong_primers(self):

        template = pydna_parse('examples/sequences/pFA6a-hphMX6.gb')[0]
        json_seq = format_sequence_genbank(template)
        json_seq.id = 1

        submitted_source = PCRSource(
            input=[1],
            primer_annealing_settings=PrimerAnnealingSettings(minimum_annealing=13)
        )

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
        data = {'source': submitted_source.dict(), 'sequences': [json_seq.dict()], 'primers': [primer_fwd.dict(), primer_rvs.dict()]}
        response = client.post("/pcr", json=data)
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

        data = {'source': submitted_source.dict(), 'sequences': [json_seq.dict()], 'primers': [primer_fwd.dict(), primer_rvs.dict()]}
        response = client.post("/pcr", json=data)
        payload = response.json()
        self.assertEqual(response.status_code, 200)

        # This is the wrong annealing info
        submitted_source.primers = [2, 3]
        submitted_source.fragment_boundaries = [59, 1200]  # 1200 instead of 1718
        submitted_source.primer_footprints = [21, 20]

        data = {'source': submitted_source.dict(), 'sequences': [json_seq.dict()], 'primers': [primer_fwd.dict(), primer_rvs.dict()]}
        response = client.post("/pcr", json=data)
        payload = response.json()
        self.assertEqual(response.status_code, 400)
        self.assertEqual(payload['detail'], 'The annealing positions of the primers seem to be wrong.')

        # This is the wrong annealing info
        submitted_source.primers = [2, 3]
        submitted_source.fragment_boundaries = [59, 1718]

        # Wrong footprint (To return this error, the 'wrong one' cannot mess with the annealing settings
        # e.g. it could not be [21,40] since it would not align with footprint 20)
        submitted_source.primer_footprints = [21, 12]
        data = {'source': submitted_source.dict(), 'sequences': [json_seq.dict()], 'primers': [primer_fwd.dict(), primer_rvs.dict()]}
        response = client.post("/pcr", json=data)
        payload = response.json()
        self.assertEqual(response.status_code, 400)
        self.assertEqual(payload['detail'], 'The annealing positions of the primers seem to be wrong.')


if __name__ == "__main__":
    unittest.main()
