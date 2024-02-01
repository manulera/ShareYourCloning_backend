from pydna.dseqrecord import Dseqrecord
import assembly2 as assembly
from Bio.Restriction import EcoRI

f1 = Dseqrecord('aaGAATTCtttGAATTCaa', circular=False)
algo = lambda x, y, l : assembly.restriction_ligation_overlap(x, y, [EcoRI], False)
f = assembly.SingleFragmentAssembly([f1], algorithm=algo)
for a in f.get_insertion_assemblies():
    print(assembly.assembly2str(a))


f1 = Dseqrecord('AGAGACCaaaAGAGACC')
f = assembly.SingleFragmentAssembly([f1], limit=7)
for a in f.get_insertion_assemblies():
    print(assembly.assembly2str(a))
    