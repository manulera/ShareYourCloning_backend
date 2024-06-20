from primer_design import homologous_recombination_primers
from Bio.SeqFeature import SimpleLocation
from Bio.Seq import reverse_complement
from unittest import TestCase
from pydna.dseqrecord import Dseqrecord
from pydna.tm import tm_default


class TestHomologousRecombinationPrimers(TestCase):

    def test_linear_sequence(self):
        """
        Test the homologous_recombination_primers function.
        """

        #                      >>>>>>>>>>  <<<<<<<<<<
        pcr_seq = Dseqrecord('GAAATGGAACAGTGCCAGAAATTTTT')
        pcr_loc = SimpleLocation(23, 43)
        hr_seq = Dseqrecord('AAACGTTT')
        print(tm_default('ATGGAACAGT'))
        # First we replace the CG
        hr_loc = SimpleLocation(3, 5)
        homology_length = 3
        minimal_hybridization_length = 10
        insert_forward = True
        target_tm = 30

        primers = homologous_recombination_primers(
            pcr_seq, pcr_loc, hr_seq, hr_loc, homology_length, minimal_hybridization_length, insert_forward, target_tm
        )
        print(primers)
        print(reverse_complement(primers[1]))
