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

b2 = Dseqrecord("CCCCCtattctggctgtatctTTTTTtgtgctgtgctcta", name="two35")

b3 = Dseqrecord("tgtgctgtgctctaTTTTTtgtgctgtgctctaCCCCtattctggctgtatct", name="two35")
# tattctggctgtatct 16
c = Dseqrecord("tattctggctgtatctGGGGGTacgatgctatactgg", name="three37").reverse_complement()
c2 = Dseqrecord("tattctggctgtatctGGGGGTcccccccccccccccccc", name="three37").reverse_complement()

asm = Assembly((a, b3, c), limit=14, use_fragment_order=True)

for a in asm.get_linear_assemblies():
    print(a)

for a in asm.get_circular_assemblies():
    print(a)
    # print(asm.execute_assembly(a).seq)


# print(asm.get_circular_assemblies())


# Node to apply the constrain that circular seqs should start with 0
# # asm.G.add_edge('circ', 0, pair=(None, None))

# print(asm.G.nodes)
# print(asm.G.edges)

# print(list(filter(lambda x: len(x) == len(set(x)) and len(x) == 4, list(_nx.cycles.find_cycle(asm.G, orientation='original')))))

# # asm.G.add_node('begin', seq=Dseqrecord(''))
# # asm.G.add_node('end', seq=Dseqrecord(''))


# asm.G.add_edge('begin', 1, pair=(None, None))
# asm.G.add_edge('begin', 2, pair=(None, None))

# asm.G.add_edge(0, 'end', pair=(None, None))
# asm.G.add_edge(1, 'end', pair=(None, None))
# asm.G.add_edge(2, 'end', pair=(None, None))
# asm.G.add_edge(-0, 'end', pair=(None, None))
# asm.G.add_edge(-1, 'end', pair=(None, None))
# asm.G.add_edge(-2, 'end', pair=(None, None))



# print(list(filter(lambda x: len(x) == 5, list(_nx.all_simple_paths(asm.G, 'begin', 'end')))))


# print(asm.G.edges())
