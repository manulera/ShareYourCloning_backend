from pydna.dseqrecord import Dseqrecord
from assembly2 import Assembly, assembly2str


homology = 'ATGCAAACAGTAATGATGGATGACATTCAAAGCACTGATT'
template = Dseqrecord(f'aaaaaa{homology}aattggaattttttt', circular=False)

# This will not give a circular assembly
insert = Dseqrecord(f'{homology}cccc{homology}', circular=False)

# This will give a circular assembly
insert = Dseqrecord(f'{homology}acaa{homology}', circular=False)

asm = Assembly((template, insert), limit=40, use_all_fragments=True)


# The condition is that the first and last fragments are the template
possible_assemblies = [a for a in asm.get_insertion_assemblies() if a[0][0] == 1]

for a in possible_assemblies:
    print(assembly2str(a))

for a in asm.assemble_insertion():
    print(len(a.seq))

for a in asm.assemble_circular():
    print(len(a.seq))
