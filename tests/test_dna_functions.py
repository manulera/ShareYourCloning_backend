from shareyourcloning.dna_functions import find_sequence_regex, custom_file_parser
from shareyourcloning.dna_utils import sum_is_sticky
import unittest
from pydna.dseq import Dseq
import os

test_files = os.path.join(os.path.dirname(__file__), 'test_files')


# Tests for sum_is_sticky() in dna_functions.py
class TestPartialSticky(unittest.TestCase):
    # General test functions
    def expectTrue(self, seq_left, seq_right, partial):
        with self.subTest():
            self.assertTrue(sum_is_sticky(seq_left.three_prime_end(), seq_right.five_prime_end(), partial))

    def expectFalse(self, seq_left, seq_right, partial):
        with self.subTest():
            self.assertFalse(sum_is_sticky(seq_left.three_prime_end(), seq_right.five_prime_end(), partial))


class MultiTestPartialSticky(TestPartialSticky):
    # Specific cases
    def test_blunt_ends(self):
        for partial in [False, True]:
            self.expectFalse(Dseq('ACGT'), Dseq('ACGT'), partial)

    def test_sticky_ends_full_overlap_3(self):
        seq1 = Dseq('ACGTAAA', 'ACGT', ovhg=0)
        seq2 = Dseq('ACGT', 'ACGTTTT', ovhg=3)

        for partial in [False, True]:
            self.expectTrue(seq1, seq2, partial)

    def test_sticky_ends_full_overlap_5(self):
        seq1 = Dseq('ACGT', 'TTTACGT', ovhg=0)
        seq2 = Dseq('AAAACGT', 'ACGT', ovhg=-3)

        for partial in [False, True]:
            self.expectTrue(seq1, seq2, partial)

    def test_sticky_ends_partial_overlap_3(self):
        seq1 = Dseq('ACGTAA', 'ACGT', ovhg=0)
        seq2 = Dseq('ACGT', 'ACGTTTT', ovhg=3)

        self.expectTrue(seq1, seq2, True)
        self.expectFalse(seq1, seq2, False)

        seq3 = Dseq('ACGTAAA', 'ACGT', ovhg=0)
        seq4 = Dseq('ACGT', 'ACGTTT', ovhg=2)

        self.expectTrue(seq3, seq4, True)
        self.expectFalse(seq3, seq4, False)

    def test_sticky_ends_partial_overlap_5(self):
        seq1 = Dseq('ACGT', 'TTACGT', ovhg=0)
        seq2 = Dseq('AAAACGT', 'ACGT', ovhg=-3)

        self.expectTrue(seq1, seq2, True)
        self.expectFalse(seq1, seq2, False)

        seq3 = Dseq('ACGT', 'TTTACGT', ovhg=0)
        seq4 = Dseq('AAACGT', 'ACGT', ovhg=-2)

        self.expectTrue(seq3, seq4, True)
        self.expectFalse(seq3, seq4, False)

    def test_sticky_ends_max_len(self):
        # Ensures that all possible overlapping lengths are covered
        seq1 = Dseq('ACGT', 'GTACGT', ovhg=0)
        seq2 = Dseq('ACAACGT', 'ACGT', ovhg=-3)

        self.expectTrue(seq1, seq2, True)
        self.expectFalse(seq1, seq2, False)


class SequenceRegexTest(unittest.TestCase):
    def test_regex(self):

        # Features spanning the whole sequence
        regex_pattern = 'AA.*AA'
        template_seq = 'AATTAA'
        template_seq2 = 'TTAATT'
        for circular in [False, True]:
            features = find_sequence_regex(regex_pattern, template_seq, circular)
            self.assertEqual(len(features), 1)
            self.assertEqual(features[0].start, 0)
            self.assertEqual(features[0].end, 6)
            self.assertEqual(features[0].strand, 1)

            # Find in the reverse strand
            features = find_sequence_regex(regex_pattern, template_seq2, circular)

            self.assertEqual(len(features), 1)
            self.assertEqual(features[0].start, 0)
            self.assertEqual(features[0].end, 6)
            self.assertEqual(features[0].strand, -1)

        # Nested features are found and returned in the correct order
        regex_pattern = 'AA.*AA'
        template_seq = 'AATTAATTAA'
        for circular in [False, True]:
            features = find_sequence_regex(regex_pattern, template_seq, circular)
            self.assertEqual(len(features), 4)

            # First AATTAA
            self.assertEqual([features[0].start, features[0].end], [0, 6])
            self.assertEqual(features[0].extract(template_seq), 'AATTAA')
            # Entire sequence
            self.assertEqual([features[1].start, features[1].end], [0, 10])
            self.assertEqual(features[1].extract(template_seq), 'AATTAATTAA')
            # reverse match
            self.assertEqual([features[2].start, features[2].end], [2, 8])
            self.assertEqual(features[2].extract(template_seq), 'AATTAA')
            # Second AATTAA
            self.assertEqual([features[3].start, features[3].end], [4, 10])
            self.assertEqual(features[3].extract(template_seq), 'AATTAA')

        # Features that span the origin, the order in which they are returned is
        # a bit arbitrary, see the documentation of location_sorter
        regex_pattern = 'AA.*CC'
        template_seq = 'TTCCTTAAGG'
        features = find_sequence_regex(regex_pattern, template_seq, False)
        self.assertEqual(len(features), 0)
        features = find_sequence_regex(regex_pattern, template_seq, True)
        self.assertEqual(len(features), 3)

        # match: AAGGAACC
        # TTCCTTAAGG
        # <<<<<<--<<
        f1, f2 = features[0].parts
        self.assertEqual([f2.start, f2.end], [8, 10])
        self.assertEqual([f1.start, f1.end], [0, 6])
        features[0].strand = -1
        self.assertEqual(features[0].extract(template_seq), 'AAGGAACC')

        # match: AACC
        # TTCCTTAAGG
        # <<------<<
        f1, f2 = features[1].parts
        self.assertEqual([f1.start, f1.end], [0, 2])
        self.assertEqual([f2.start, f2.end], [8, 10])
        self.assertEqual(features[1].extract(template_seq), 'AACC')

        # match: AAGGTTCC
        # TTCCTTAAGG
        # >>>>-->>>>
        f1, f2 = features[2].parts
        self.assertEqual([f1.start, f1.end], [6, 10])
        self.assertEqual([f2.start, f2.end], [0, 4])
        self.assertEqual(features[2].extract(template_seq), 'AAGGTTCC')


class TestPermisiveParserWithApe(unittest.TestCase):
    def test_permisive_parser_with_ape_circular(self):
        with open(f'{test_files}/P2RP3.ape', 'r') as f:
            plasmid = custom_file_parser(f, 'genbank')[0]
            # Since APE files are not correctly gb formatted (as of 2024-11-27)
            # the Bio.SeqIO.parse may not recognize the topology of the plasmid
            # Our custom permissive parser should be then used and the topology
            # parameter properly recognized
            self.assertEqual(plasmid.circular, True)

    def test_permisive_parser_with_ape_linear(self):
        with open(f'{test_files}/P2RP3_linear.ape', 'r') as f:
            # I manually changed the topology of the plasmid to linear
            plasmid = custom_file_parser(f, 'genbank')[0]
            self.assertEqual(plasmid.circular, False)

    def test_permisive_parser_no_topology(self):
        with open(f'{test_files}/ase1_no_topology.gb', 'r') as f:
            plasmid = custom_file_parser(f, 'genbank')[0]
            self.assertEqual(plasmid.circular, False)
