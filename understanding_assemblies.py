from pydna.assembly import Assembly
from pydna.dseqrecord import Dseqrecord
from dna_functions import get_homologous_recombination_locations

template = Dseqrecord('TTTCGTCCTGAAAAATGTCGAAGT', circular=False)
insert = Dseqrecord('CGTCCTGAttttTGTCGAAGT', circular=False)

products = Assembly([template, insert, template], limit=6)

print(products.assemble_linear())
products = get_homologous_recombination_locations(template, insert, 6)

print(products)
