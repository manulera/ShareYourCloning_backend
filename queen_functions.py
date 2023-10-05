from QUEEN.qobj import QUEEN, DNAfeature
from QUEEN.qfunction import cutdna, cropdna, modifyends, joindna
from Bio.Seq import reverse_complement
from pydna.dseqrecord import Dseqrecord
from pydna.parsers import parse


def queen_from_dseqrecord(dseqerec: Dseqrecord) -> QUEEN:
    queen_seq = QUEEN(seq=str(dseqerec.seq), topology='circular' if dseqerec.circular else 'linear')
    queen_seq._dnafeatures = [DNAfeature(f, subject=queen_seq) for f in dseqerec.features]
    return queen_seq


dseqrec = parse('examples/sequences/pFA6a-hphMX6.gb')[0]
queen_seq = queen_from_dseqrecord(dseqrec)
queen_seq.printsequence(display=True, hide_middle=10)
queen_seq.printfeature()
queen_seq.searchsequence
