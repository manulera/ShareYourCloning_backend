from primer_design import homologous_recombination_primers
from Bio.SeqFeature import SimpleLocation, SeqFeature
from unittest import TestCase
from pydna.dseqrecord import Dseqrecord


class TestHomologousRecombinationPrimers(TestCase):

    def test_normal_examples(self):
        """
        Test the homologous_recombination_primers function.
        """

        #                      >>>>>>>>>>> <<<<<<<<<<
        pcr_seq = Dseqrecord('GAAATGGAACAGTGCCAGAAATTTTT')
        pcr_loc = SimpleLocation(1, 23)
        hr_seq = Dseqrecord('AAACGTTT')
        homology_length = 3
        minimal_hybridization_length = 10
        insert_forward = True
        target_tm = 30

        # First we replace the CG
        hr_loc_replace = SimpleLocation(3, 5)
        primers = homologous_recombination_primers(
            pcr_seq,
            pcr_loc,
            hr_seq,
            hr_loc_replace,
            homology_length,
            minimal_hybridization_length,
            insert_forward,
            target_tm,
        )
        self.assertEqual(primers, ('aaaAAATGGAACAG', 'aaaAATTTCTGGC'))

        # An insertion
        hr_loc_insert = SimpleLocation(3, 3)
        primers = homologous_recombination_primers(
            pcr_seq,
            pcr_loc,
            hr_seq,
            hr_loc_insert,
            homology_length,
            minimal_hybridization_length,
            insert_forward,
            target_tm,
        )
        self.assertEqual(primers, ('aaaAAATGGAACAG', 'acgAATTTCTGGC'))

        # Same, circular and we loop
        pcr_seq = pcr_seq.looped()
        hr_seq = hr_seq.looped()
        pcr_seq.features.append(SeqFeature(pcr_loc))
        hr_seq.features.append(SeqFeature(hr_loc_replace))
        hr_seq.features.append(SeqFeature(hr_loc_insert))

        for shift_pcr in range(len(pcr_seq)):
            pcr_shifted = pcr_seq.shifted(shift_pcr)
            pcr_loc = pcr_shifted.features[0].location

            hr_shifted = hr_seq
            hr_loc_replace = hr_shifted.features[0].location
            hr_loc_insert = hr_shifted.features[1].location
            solutions = (
                ('aaaAAATGGAACAG', 'aaaAATTTCTGGC'),
                ('aaaAAATGGAACAG', 'acgAATTTCTGGC'),
            )
            for shift_hr in range(len(hr_seq)):
                hr_shifted = hr_seq.shifted(shift_hr)
                hr_loc_replace = hr_shifted.features[0].location
                hr_loc_insert = hr_shifted.features[1].location
                solutions = (
                    ('aaaAAATGGAACAG', 'aaaAATTTCTGGC'),
                    ('aaaAAATGGAACAG', 'acgAATTTCTGGC'),
                )
                for hr_loc, solution in zip([hr_loc_replace, hr_loc_insert], solutions):
                    for insert_forward in (True, False):
                        if not insert_forward:
                            solution = (
                                solution[0][:homology_length] + solution[1][homology_length:],
                                solution[1][:homology_length] + solution[0][homology_length:],
                            )

                        primers = homologous_recombination_primers(
                            pcr_shifted,
                            pcr_loc,
                            hr_shifted,
                            hr_loc,
                            homology_length,
                            minimal_hybridization_length,
                            insert_forward,
                            target_tm,
                        )
                        self.assertEqual(primers, solution)

    def test_clashing_homology(self):

        pcr_seq = Dseqrecord('GAAATGGAACAGTGCCAGAAATTTTT')
        pcr_loc = SimpleLocation(1, 23)
        hr_seq = Dseqrecord('AAACGTTT')
        homology_length = 5
        minimal_hybridization_length = 10
        insert_forward = True
        target_tm = 30
        hr_loc_insert = SimpleLocation(3, 3)

        try:
            homologous_recombination_primers(
                pcr_seq,
                pcr_loc,
                hr_seq,
                hr_loc_insert,
                homology_length,
                minimal_hybridization_length,
                insert_forward,
                target_tm,
            )
            self.fail('Expected ValueError.')

        except ValueError as e:
            self.assertEqual(str(e), 'Forward homology region is out of bounds.')

        hr_loc_insert = SimpleLocation(5, 5)

        try:
            homologous_recombination_primers(
                pcr_seq,
                pcr_loc,
                hr_seq,
                hr_loc_insert,
                homology_length,
                minimal_hybridization_length,
                insert_forward,
                target_tm,
            )
            self.fail('Expected ValueError.')
        except ValueError as e:
            self.assertEqual(str(e), 'Reverse homology region is out of bounds.')

        # circular case with clashing homology
        hr_seq = hr_seq.looped()
        homology_length = 4
        # exact match
        homologous_recombination_primers(
            pcr_seq,
            pcr_loc,
            hr_seq,
            hr_loc_insert,
            homology_length,
            minimal_hybridization_length,
            insert_forward,
            target_tm,
        )
        # clashing homology
        homology_length = 5
        try:
            homologous_recombination_primers(
                pcr_seq,
                pcr_loc,
                hr_seq,
                hr_loc_insert,
                homology_length,
                minimal_hybridization_length,
                insert_forward,
                target_tm,
            )
            self.fail('Expected ValueError.')
        except ValueError as e:
            self.assertEqual(str(e), 'Homology arms overlap.')
