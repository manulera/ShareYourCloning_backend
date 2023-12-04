# %%
from assembly2 import Assembly
from pydna.assembly import Assembly as Assembly_pydna
from assembly2 import example_fragments, linear_results
from dna_functions import get_homologous_recombination_locations, perform_homologous_recombination
from pydna.dseqrecord import Dseqrecord
from pydna.common_sub_strings import common_sub_strings
import networkx as _nx
from Bio.SeqFeature import SeqFeature, SimpleLocation
from pydna.parsers import parse
from pydna.amplify import pcr
from pydna.readers import read


primer = parse("test_files/primers.fas", ds=False)

a = parse('dummy_a.gb')[0]
b = parse('dummy_b.gb')[0]

asm = Assembly([a, b], limit=30)

candidates = asm.assemble_circular()

new_assembly = set()
for c in candidates:
    print(c.seq.cseguid())
    new_assembly.add(c.seq.cseguid())

asm2 = Assembly_pydna([a, b], limit=30)

print('==')
candidates2 = asm2.assemble_circular()

for c in candidates2:
    if c.seq.cseguid() not in new_assembly:
        print(c.seq.cseguid())
        new_assembly.add(c.seq.cseguid())
        # c.write(f'{c.seq.cseguid()}.gb')
        print(c.figure())

