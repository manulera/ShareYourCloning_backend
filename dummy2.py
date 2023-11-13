# %%
from assembly2 import Assembly
from dna_functions import get_homologous_recombination_locations, perform_homologous_recombination
from pydna.dseqrecord import Dseqrecord
from pydna.common_sub_strings import common_sub_strings
import networkx as _nx

# acgatgctatactgg 15
a = Dseqrecord("acgatgctatactggCCCCCtgtgctgtgctctaGG", name="one36")
# tgtgctgtgctcta 14
b = Dseqrecord("tgtgctgtgctctaTTTTTtattctggctgtatct", name="two35")
# tattctggctgtatct 16
c = Dseqrecord("tattctggctgtatctGGGGGTacgatgctatactgg", name="three37").reverse_complement()

asm = Assembly((a, b, c), limit=14)

print(list(_nx.cycles.find_cycle(asm.G, orientation='ignore')))

# asm.G.add_node('begin', seq=Dseqrecord(''))
# asm.G.add_node('end', seq=Dseqrecord(''))

asm.G.add_edge('begin', 0, pair=(None, None))
asm.G.add_edge('begin', 1, pair=(None, None))
asm.G.add_edge('begin', 2, pair=(None, None))

asm.G.add_edge(0, 'end', pair=(None, None))
asm.G.add_edge(1, 'end', pair=(None, None))
asm.G.add_edge(2, 'end', pair=(None, None))
asm.G.add_edge(-0, 'end', pair=(None, None))
asm.G.add_edge(-1, 'end', pair=(None, None))
asm.G.add_edge(-2, 'end', pair=(None, None))



print(list(_nx.all_simple_paths(asm.G, 'begin', 'end')))


print(asm.G.edges())
