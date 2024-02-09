from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
import assembly2 as assembly

template = Dseqrecord('AATTAGCAGCGATCGAGT')
# Mismatches

primer = Dseqrecord('GCGATCGAAAAA')
assembly.alignment_sub_strings(template, primer, True, 8, 0)

primer = Dseqrecord('GCGTTCGA')
assembly.alignment_sub_strings(template, primer, True, 8, 1)