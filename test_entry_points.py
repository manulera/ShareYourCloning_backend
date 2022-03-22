from typing import List
from dna_functions import format_sequence_genbank, read_dsrecord_from_json
from main import app
from fastapi.testclient import TestClient
from pydna.parsers import parse as pydna_parse
from Bio.Restriction.Restriction import CommOnly
from pydantic_models import PCRSource, PrimerAnnealing, PrimerAnnealingSettings, PrimerModel, PrimerPair, SequenceEntity, StickyLigationSource
from pydna.dseqrecord import Dseqrecord
import unittest
from pydna.dseq import Dseq

client = TestClient(app)

# TODO further tests are needed (combinations)


class StickyLigationTest(unittest.TestCase):

    def test_sticky_ligation(self):
        """Test whether assembly is executed when the order is provided"""
        # Load dummy sequence

        initial_sequence = pydna_parse('examples/sequences/dummy_EcoRI.fasta')[0]

        # Restriction cut
        enzyme = CommOnly.format('EcoRI')
        output_list: List[Dseqrecord] = initial_sequence.cut([enzyme])

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
        resulting_sequence = read_dsrecord_from_json(SequenceEntity.parse_obj(payload['sequence']))

        # Check that the assembly is correct
        self.assertEqual(resulting_sequence.seq, initial_sequence.seq)
        resulting_source = StickyLigationSource.parse_obj(payload['source'])

        # No extra field should have been set on the source
        self.assertEqual(resulting_source, source)

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

        source1 = PCRSource.parse_obj(payload['source'])

        # Only one possible output
        self.assertEqual(len(source1.output_list), 1)
        self.assertEqual(len(source1.possible_primer_pairs), 1)

        # The sequence matches what we expect
        dseq1 = read_dsrecord_from_json(source1.output_list[0])
        primer_pair = source1.possible_primer_pairs[0]
        predicted_seq = primer_fwd.sequence + template[primer_pair.forward.position:primer_pair.reverse.position] + Dseq(primer_rvs.sequence).reverse_complement()
        self.assertEqual(dseq1.seq, predicted_seq.seq)
        self.assertEqual(primer_pair.forward.primer_id, primer_fwd.id)
        self.assertEqual(primer_pair.reverse.primer_id, primer_rvs.id)

        # Now we submit the deterministic PCR (we already know which fragment we want)
        submitted_source.primer_pair = primer_pair
        data = {'source': submitted_source.dict(), 'sequences': [json_seq.dict()], 'primers': [primer_fwd.dict(), primer_rvs.dict()]}
        response = client.post("/pcr", json=data)
        payload = response.json()
        source2 = PCRSource.parse_obj(payload['source'])
        dseq2 = read_dsrecord_from_json(SequenceEntity.parse_obj(payload['sequence']))
        self.assertEqual(dseq2.seq, predicted_seq.seq)
        self.assertEqual(len(source2.output_list), 0)
        self.assertEqual(len(source2.possible_primer_pairs), 0)

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

        # This is the correct annealing info
        submitted_source.primer_pair = PrimerPair(
            forward=PrimerAnnealing(primer_id=2, position=59, footprint_length=21),
            reverse=PrimerAnnealing(primer_id=3, position=1718, footprint_length=20)
        )

        data = {'source': submitted_source.dict(), 'sequences': [json_seq.dict()], 'primers': [primer_fwd.dict(), primer_rvs.dict()]}
        response = client.post("/pcr", json=data)
        payload = response.json()
        self.assertEqual(response.status_code, 200)

        # This is the wrong annealing info
        submitted_source.primer_pair = PrimerPair(
            forward=PrimerAnnealing(primer_id=2, position=59, footprint_length=21),
            reverse=PrimerAnnealing(primer_id=3, position=1710, footprint_length=20)
        )

        data = {'source': submitted_source.dict(), 'sequences': [json_seq.dict()], 'primers': [primer_fwd.dict(), primer_rvs.dict()]}
        response = client.post("/pcr", json=data)
        payload = response.json()
        self.assertEqual(response.status_code, 400)
        self.assertEqual(payload['detail'], 'The annealing positions of the primers seem to be wrong.')
