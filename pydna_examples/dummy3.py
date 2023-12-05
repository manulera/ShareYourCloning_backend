from assembly2 import Assembly
from pydna.assembly import Assembly as Assembly_pydna
from pydna.dseqrecord import Dseqrecord

a = Dseqrecord('ttttttCCtttttt')
b = Dseqrecord('ttttttGGtttttt')

asm = Assembly([a, b], limit=6)

print('>> circular')
for c in asm.assemble_circular():
    print(c.seq)

print('>> linear')

for c in asm.assemble_linear():
    print(c.seq)

print()
print('> Old')

asm = Assembly_pydna([a, b], limit=6)

print('>> circular')
for c in asm.assemble_circular():
    print(c.seq)
    print(c.figure())


print('>> linear')
for c in asm.assemble_linear():
    print(c.seq)

# for c in asm.get_linear_assemblies():
#     print(c)
#     print(asm.execute_assembly(c).seq)

