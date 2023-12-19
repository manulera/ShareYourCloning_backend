from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from Bio.Restriction.Restriction import RestrictionBatch
from assembly2 import restriction_ligation_overlap, Assembly, assemble, assembly2str
from itertools import product

a = Dseqrecord('AAAGAATTCAAA')
b = Dseqrecord('CCCCGAATTCCCC')
enzymes = RestrictionBatch(first=['EcoRI'])
algo = lambda x, y, l : restriction_ligation_overlap(x, y, enzymes)
f = Assembly([a, b], algorithm=algo, use_fragment_order=False)
for assembly in f.get_linear_assemblies():
    print(assembly2str(assembly))
    out = assemble([a, b], assembly, False)
    print(out.seq, out.reverse_complement().seq)
exit()
print(restriction_ligation_overlap(a, b, enzymes))
print(a.seq.get_cutsites(enzymes))
print(b.seq.get_cutsites(enzymes))


a = Dseqrecord('AAAGCGATCGCAAA')
b = Dseqrecord('CCCCGCGATCGCCCC')
enzymes = RestrictionBatch(first=['RgaI'])
algo = lambda x, y, l : restriction_ligation_overlap(x, y, enzymes)
f = Assembly([a, b], algorithm=algo)
print(f.assemble_linear()[0].seq)
print(f.assemble_linear()[1].seq)

print(restriction_ligation_overlap(a, b, enzymes))
print(a.seq.get_cutsites(enzymes))
print(b.seq.get_cutsites(enzymes))
