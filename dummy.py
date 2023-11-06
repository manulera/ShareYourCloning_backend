import regex
from pydna.parsers import parse

seqrecord = parse('examples/sequences/addgene-plasmid-39296-sequence-49545.gbk')
print('a')
list(regex.finditer(r'AA.*AA', str(seqrecord[0].seq), overlapped=True))
print('b')
