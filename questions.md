* Same molecule returned multiple times
* Same fragment used several times (different edges?), does this make sense? In the end same fragment multiple times
would give infinite sequences for those that can assemble circularly.

```python
primer = parse("test_files/primers.fas", ds=False)

a = parse('test_files/dummy_a.gb')[0]
b = parse('tes_files/dummy_b.gb')[0]

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
```