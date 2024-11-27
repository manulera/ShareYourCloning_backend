import unittest
from dna_functions import custom_file_parser


class TestPermisiveParserWithApe(unittest.TestCase):
    def test_permisive_parser_with_ape_circular(self):
        with open('test_files/P2RP3.ape', 'r') as f:
            plasmid = custom_file_parser(f, 'gb')[0]
            # Since APE files are not correctly gb formatted (as of 2024-11-27)
            # the Bio.SeqIO.parse may not recognize the topology of the plasmid
            # Our custom permissive parser should be then used and the topology
            # parameter properly recognized
            self.assertEqual(plasmid.circular, True)

    def test_permisive_parser_with_ape_linear(self):
        with open('test_files/P2RP3_linear.ape', 'r') as f:
            # I manually changed the topology of the plasmid to linear
            plasmid = custom_file_parser(f, 'gb')[0]
            self.assertEqual(plasmid.circular, False)
