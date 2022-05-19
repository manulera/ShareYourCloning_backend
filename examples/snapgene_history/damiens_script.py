#!/usr/bin/env python3
# Code by @gouttegd, see #12


from Bio.SeqIO.SnapGeneIO import _iterate
from lzma import decompress, FORMAT_XZ
from xml.etree.ElementTree import parse, indent
from io import BytesIO
import sys

XZ_SIGNATURE = b"\xFD\x37\x7A\x58\x5A\x00"
XML_SIGNATURE = b"\x3C\x3F\x78\x6D\x6C\x20"

snapfile = open(sys.argv[1], 'rb')
iterator = _iterate(snapfile)

for (ptype, length, data) in iterator:
    if ptype != 7: # Not a history packet
        continue

    if data[0:6] == XZ_SIGNATURE:
        data = decompress(data, format=FORMAT_XZ)

    if data[0:6] == XML_SIGNATURE:
        with BytesIO(data) as xml:
            history = parse(xml)
        indent(history)
        history.write(sys.stdout, encoding='unicode', xml_declaration=True)
    else:
        print("Unexpected payload in history packet", file=sys.stderr)

snapfile.close()
print(history)