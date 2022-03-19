from typing import List
from dna_functions import format_sequence_genbank, read_dsrecord_from_json
from main import app
from fastapi.testclient import TestClient
from pydna.parsers import parse as pydna_parse
from Bio.Restriction.Restriction import CommOnly
from pydantic_models import SequenceEntity, StickyLigationSource
from pydna.dseqrecord import Dseqrecord
import unittest

client = TestClient(app)


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
