# %%
from dna_functions import get_homologous_recombination_locations, perform_homologous_recombination
from pydna.dseqrecord import Dseqrecord
from pydna.assembly import Assembly
from pydna.common_sub_strings import common_sub_strings
from pydna.contig import Contig
import networkx as _nx


# acgatgctatactgg 15
a = Dseqrecord("acgatgctatactggCCCCCtgtgctgtgctctaGG", name="one36")
# tgtgctgtgctcta 14
b = Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatct", name="two35")
# tattctggctgtatct 16
c = Dseqrecord("tattctggctgtatctGGGGGTacgatgctatactgg", name="three37").reverse_complement()

asm = Assembly((a, c, b), limit=14)

# contig : Contig = asm.assemble_linear()[0]

# print(contig.seq)
# print(contig.detailed_figure())
# print(asm.G.nodes)


# asm = Assembly((a, b, c), limit=14)
# print(asm.assemble_linear())
# contig : Contig = asm.assemble_linear()[0]


# print(contig.seq)
# print(contig.detailed_figure())
# print(asm.G.nodes)

