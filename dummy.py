from pydna.dseqrecord import Dseqrecord
import assembly2 as assembly
from Bio.Restriction import BsaI

# Partial overlaps -> enzyme with negative overhang
fragments = [Dseqrecord('GGTCTCCCCAATT'), Dseqrecord('GGTCTCCAACCAA'), Dseqrecord('ggTCTCCCCAATT')]


# Allowing partial overlaps
def algo(x, y, _l):
    return assembly.restriction_ligation_overlap(x, y, [BsaI], partial=True)


f = assembly.Assembly(fragments, algorithm=algo, use_fragment_order=False)

print(*f.G.edges, sep='\n')
