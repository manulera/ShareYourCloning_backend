from pydna.dseqrecord import Dseqrecord
import assembly2 as assembly
from Bio.Restriction import BsaI

# Partial overlaps -> enzyme with negative overhang
fragments = [Dseqrecord('GGTCTCCCCAATT'), Dseqrecord('GGTCTCCAACCAA')]


# Allowing partial overlaps
def algo(x, y, _l):
    return assembly.restriction_ligation_overlap(x, y, [BsaI], partial=True)


f = assembly.Assembly(fragments, algorithm=algo, use_fragment_order=False)

a1, a2 = fragments[0].cut(BsaI)
b1, b2 = fragments[1].cut(BsaI)

print(repr(a1.seq))
print(repr(a2.seq))
print(repr(b1.seq))
print(repr(b2.seq))

print(*f.G.edges, sep='\n')

print(*[assembly.assembly2str(f) for f in f.get_linear_assemblies()], sep='\n')
