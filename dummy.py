import assembly2 as assembly
from Bio.SeqFeature import SeqFeature, SimpleLocation
from pydna.dseqrecord import Dseqrecord

f1 = Dseqrecord('ATTTA', circular=True)
f1.features = [SeqFeature(SimpleLocation(1, 4))]
f1_shifted = f1.shifted(2)
print(f1_shifted.seq)
for f in f1_shifted.features:
    print(f)
dummy_cut = 
print(f1_shifted.apply_cut())