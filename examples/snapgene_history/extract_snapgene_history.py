# %%
from Bio.SeqIO.SnapGeneIO import _iterate

snapgene_string = ''

with open("final_plasmid.dna", "rb") as in_handle:
    for i in _iterate(in_handle):
        try:
            snapgene_string += i[2].decode("utf-8")
        except UnicodeDecodeError:
            pass

with open("final_plasmid_history.xml", "w") as out_handle:
    out_handle.write(snapgene_string)
