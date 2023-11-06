from dna_functions import dseq_from_both_overhangs, both_overhangs_from_dseq, \
    format_sequence_genbank, read_dsrecord_from_json, sum_is_sticky, find_sequence_regex
import unittest
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from Bio.SeqFeature import SimpleLocation, SeqFeature, CompoundLocation
from typing import OrderedDict


class DseqFromBothOverhangsTest(unittest.TestCase):

    def test_conversion(self):
        # Here we try both with longer watson and longer crick
        for watson, crick in [('AAAAAA', 'TTTT'), ('TTTT', 'AAAAAA')]:
            for ovhg in [-2, 0, 3]:
                with self.subTest(ovhg=ovhg):
                    dseq_original = Dseq(watson, crick=crick, ovhg=ovhg)

                    crick_overhang_3p, watson_overhang_3p = both_overhangs_from_dseq(
                        dseq_original)
                    dseq_2 = dseq_from_both_overhangs(
                        str(dseq_original), crick_overhang_3p, watson_overhang_3p)

                    # We check that the elements of Dseq are transferred properly
                    self.assertEqual(dseq_original.watson,
                                     dseq_2.watson)
                    self.assertEqual(dseq_original.crick, dseq_2.crick)
                    self.assertEqual(dseq_original.ovhg, dseq_2.ovhg)

                    # Now we check for the features
                    dseq_original = Dseqrecord(dseq_original)
                    dseq_2 = Dseqrecord(dseq_2)

                    # We add some features:
                    # TODO document somewhere the fact that the strand must
                    # be specified. The file readers assume +1 strand for
                    # all features when reading from GenBank files
                    for a, start, end in [('a', 0, 2), ('b', 1, 2), ('c', 4, 7)]:
                        dseq_original.features.append(
                            SeqFeature(
                                location=SimpleLocation(start, end),
                                type="misc_feature",
                                qualifiers=OrderedDict({"label": [a]}),
                                strand=1)
                        )
                    dseq_2.features = dseq_original.features

                    # We check that the features are transferred normally
                    for i in range(len(dseq_2.features)):
                        feature_original: SeqFeature = dseq_original.features[i]
                        feature_2: SeqFeature = dseq_2.features[i]
                        self.assertEqual(
                            feature_original.extract(dseq_original),
                            feature_2.extract(dseq_2))

                    # Finally we test with pydantic models
                    seq_entity = format_sequence_genbank(dseq_original)
                    # Default value when creating DseqRecord
                    seq_entity.id = 'id'
                    dseq_3 = read_dsrecord_from_json(seq_entity)

                    self.assertEqual(dseq_original.seq.watson,
                                     dseq_3.seq.watson)
                    self.assertEqual(dseq_original.seq.crick, dseq_3.seq.crick)
                    self.assertEqual(dseq_original.seq.ovhg, dseq_3.seq.ovhg)

                    # We check that the features are transferred normally
                    for i in range(len(dseq_3.features)):
                        feature_original: SeqFeature = dseq_original.features[i]
                        feature_3: SeqFeature = dseq_3.features[i]
                        self.assertEqual(
                            feature_original.extract(dseq_original),
                            feature_3.extract(dseq_3))


# Tests for sum_is_sticky() in dna_functions.py
class TestPartialSticky(unittest.TestCase):
    # General test functions
    def expectTrue(self, seq_left, seq_right, partial):
        with self.subTest():
            self.assertTrue(sum_is_sticky(seq_left, seq_right, partial))

    def expectFalse(self, seq_left, seq_right, partial):
        with self.subTest():
            self.assertFalse(sum_is_sticky(seq_left, seq_right, partial))


class MultiTestPartialSticky(TestPartialSticky):
    # Specific cases
    def test_blunt_ends(self):
        for partial in [False, True]:
            self.expectFalse(Dseq("ACGT"), Dseq("ACGT"), partial)
    

    def test_sticky_ends_full_overlap_3(self):
        seq1 = Dseq("ACGTAAA", "ACGT", ovhg=0)
        seq2 = Dseq("ACGT", "ACGTTTT", ovhg=3)

        for partial in [False, True]:
            self.expectTrue(seq1, seq2, partial)


    def test_sticky_ends_full_overlap_5(self):
        seq1 = Dseq("ACGT", "TTTACGT", ovhg=0)
        seq2 = Dseq("AAAACGT", "ACGT", ovhg=-3)

        for partial in [False, True]:
            self.expectTrue(seq1, seq2, partial)


    def test_sticky_ends_partial_overlap_3(self):
        seq1 = Dseq("ACGTAA", "ACGT", ovhg=0)
        seq2 = Dseq("ACGT", "ACGTTTT", ovhg=3)

        self.expectTrue(seq1, seq2, True)
        self.expectFalse(seq1, seq2, False)

        seq3 = Dseq("ACGTAAA", "ACGT", ovhg=0)
        seq4 = Dseq("ACGT", "ACGTTT", ovhg=2)

        self.expectTrue(seq3, seq4, True)
        self.expectFalse(seq3, seq4, False)


    def test_sticky_ends_partial_overlap_5(self):
        seq1 = Dseq("ACGT", "TTACGT", ovhg=0)
        seq2 = Dseq("AAAACGT", "ACGT", ovhg=-3)

        self.expectTrue(seq1, seq2, True)
        self.expectFalse(seq1, seq2, False)

        seq3 = Dseq("ACGT", "TTTACGT", ovhg=0)
        seq4 = Dseq("AAACGT", "ACGT", ovhg=-2)

        self.expectTrue(seq3, seq4, True)
        self.expectFalse(seq3, seq4, False)

    def test_sticky_ends_max_len(self):
        # Ensures that all possible overlapping lengths are covered
        seq1 = Dseq("ACGT", "GTACGT", ovhg=0)
        seq2 = Dseq("ACAACGT", "ACGT", ovhg=-3)

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
        self.assertEqual([f1.start, f1.end], [8, 10])
        self.assertEqual([f2.start, f2.end], [0, 6])
        print(features[0])
        self.assertEqual(features[0].extract(template_seq), 'AAGGAACC')

        # match: AACC
        # TTCCTTAAGG
        # <<------<<
        f1, f2 = features[1].parts
        self.assertEqual([f1.start, f1.end], [8, 10])
        self.assertEqual([f2.start, f2.end], [0, 2])

        # match: AAGGTTCC
        # TTCCTTAAGG
        # >>>>-->>>>
        f1, f2 = features[2].parts
        self.assertEqual([f1.start, f1.end], [6, 10])
        self.assertEqual([f2.start, f2.end], [0, 4])
