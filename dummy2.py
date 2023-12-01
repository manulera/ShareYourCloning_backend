# %%
from assembly2 import Assembly
from pydna.assembly import Assembly as Assembly_pydna
from assembly2 import example_fragments, linear_results
from dna_functions import get_homologous_recombination_locations, perform_homologous_recombination
from pydna.dseqrecord import Dseqrecord
from pydna.common_sub_strings import common_sub_strings
import networkx as _nx
from Bio.SeqFeature import SeqFeature, SimpleLocation

# for f in example_fragments:
#     print(f.name, f.seq, '/', f.seq.reverse_complement())

# print('++++++')

# for i, rs in enumerate(linear_results):
#     print(i, rs.seq, '/', rs.seq.reverse_complement())

# print('======')
# asm_pydna = Assembly_pydna(example_fragments, limit=5)
# for f in asm_pydna.assemble_circular():
#     print(f.seq, '/', f.seq.reverse_complement())

# print('======')

# acgatgctatactgg 15
a = Dseqrecord("GGacgatgctatactggCCCCCtgtgctgtgctctaGG", name="one36")
# tgtgctgtgctcta 14                      xxxxxxxxx
b = Dseqrecord("AAtgtgctgtgctctaTTTTTacgatgctatactggCC", name="two35")

fa = SeqFeature(SimpleLocation(2, 6), type='misc_feature')
fb = SeqFeature(SimpleLocation(25, 34), type='misc_feature')
a.features.append(fa)
b.features.append(fb)




a = Dseqrecord("GGacgatgctatactggCCCCCtgtgctgtgctctaGG", name="one36")
# tgtgctgtgctcta 14
b = Dseqrecord("GGtgtgctgtgctctaTTTTTtattctggctgtatctGG", name="two35")

# tattctggctgtatct 16
c = Dseqrecord("tattctggctgtatctGGGGGTacgatgctatactgg", name="three37")


asm = Assembly((a,b,c), limit=14, use_fragment_order=True, use_all_fragments=True)

for f in asm.assemble_linear():
    print(f.seq)
# asm_exc = list(asm.get_circular_assemblies())[1]

# asm.execute_assembly(asm_exc)




# print(aa[0].features[0].extract(aa[0].seq))
# print(aa[0].features[1].extract(aa[0].seq))


# for i, f in enumerate(asm.get_circular_assemblies()):
#     print(i, f)
# for i, f in enumerate(asm.assemble_circular()):
#     print(i, f.seq, '/', f.seq.reverse_complement())

# for e in _nx.cycles.find_cycle(asm.G, source=1, orientation='original'):
#     print(e)



# a AacgatCAtgctcc / ggagcaTGatcgtT
# b TtgctccTAAattctgc / gcagaatTTAggagcaA
# c CattctgcGAGGacgatG / CatcgtCCTCgcagaatG