from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from Bio.Restriction.Restriction import RestrictionBatch
from assembly2 import restriction_ligation_overlap, Assembly, assemble, assembly2str
from itertools import product

a = Dseqrecord('AAAGAATTCTTT', circular=True)
b = Dseqrecord('TTTGAATTCAAAAGAATTCCCC')
enzymes = RestrictionBatch(first=['EcoRI'])

algo = lambda x, y, l : restriction_ligation_overlap(x, y, enzymes)
f = Assembly([b, a], algorithm=algo, use_fragment_order=False)
for assembly in f.get_circular_assemblies():
    print('>>', assembly2str(assembly))
    out = assemble([a, b], assembly, False)
    print(out.seq, out.reverse_complement().seq)
