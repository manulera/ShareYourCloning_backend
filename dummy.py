from pydna.dseqrecord import Dseqrecord
import assembly2 as assembly
from Bio.Restriction import EcoRI

f1 = Dseqrecord('GAATTCaaaGAATTC')
f2 = Dseqrecord('ATTCCCCCCCCGA', circular=True)


fragments = [f1, f2]

algo = lambda x, y, l : assembly.restriction_ligation_overlap(x, y, [EcoRI])
f = assembly.Assembly(fragments, algorithm=algo, use_fragment_order=False)

a1, a2 = f.get_circular_assemblies()

print(*f.G.edges, sep='\n')

subfragment_representation = assembly.edge_representation2subfragment_representation(a1, True)
print([f'{a[0]}, {a[1]}, {a[2]}' for a in subfragment_representation])

# Length of the overlaps between consecutive assembly fragments
fragment_overlaps = [len(e[-1]) for e in a1]

for e in a1:
    print('>', e[-1], len(e[-1]))

subfragments = assembly.get_assembly_subfragments(fragments, subfragment_representation)

for sf in subfragments:
    print(sf.seq)