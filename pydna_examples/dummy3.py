from pydna.assembly import Assembly as Assembly_pydna
from pydna.dseqrecord import Dseqrecord

a = Dseqrecord('ttttttCCtttttt')
b = Dseqrecord('ttttttGGtttttt')


asm = Assembly_pydna([a, b], limit=6)

# These circular assemblies are all circularisations of the single fragments,
# no sequence containing a and b is returned
print('>> circular')
for c in asm.assemble_circular():
    print(c.seq)
    print(c.figure())


print('>> linear')
for c in asm.assemble_linear():
    print(c.seq)

# Ignore code below
# from assembly2 import Assembly
# asm = Assembly([a, b], limit=6)

# print('>> circular')
# for c in asm.assemble_circular():
#     print(c.seq)

# print('>> linear')

# for c in asm.assemble_linear():
#     print(c.seq)

# print()



