# %%
import networkx as _nx
from pydna.common_sub_strings import common_sub_strings
from Bio.SeqFeature import SimpleLocation
from Bio.Seq import reverse_complement

G = _nx.MultiDiGraph()

G.add_node('A', blah='1')
G.add_node('B', blah='2')
G.add_node('C', blah='3')
G.add_node('D', blah='4')
G.add_edge('A', 'B', edge_info='AB')
G.add_edge('A', 'B', edge_info='AB2')
G.add_edge('B', 'C', edge_info='BC')
G.add_edge('A', 'D', edge_info='AD')

_nx.write_network_text(G)

list(_nx.all_simple_paths(G, "A", "C"))
G.get_edge_data('A', 'B')

# %%
G = _nx.MultiDiGraph()

seqA = 'AAAAAAAAAAAatgcaaacagtaatgatggaTTTTTTTTTTacaacggcaatgaatgcccaCCCCCCCCCC'
seqB = 'atgcaaacagtaatgatggaTTAAacaacggcaatgaatgccca'

G.add_node(0, seq=seqA)
G.add_node(1, seq=seqB)


matches = common_sub_strings(seqA, seqB, 20)
location_pairs = [(SimpleLocation(x_start, x_start + length, 1), SimpleLocation(y_start, y_start + length, 1)) for x_start, y_start, length in matches]

matches2 = common_sub_strings(seqA, reverse_complement(seqB), 20)
location_pairs += [(SimpleLocation(x_start, x_start + length, 1), SimpleLocation(y_start, y_start + length, -1)) for x_start, y_start, length in matches]


for pair in location_pairs:
    G.add_edge(0, 1, pair=pair)

_nx.write_network_text(G)

# %%

a = SimpleLocation(1,4,1)
b = a._flip(6)
c = b._flip(6)
print(b, c)