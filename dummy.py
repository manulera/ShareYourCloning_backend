from pydantic_models import SeqFeatureModel
from dna_functions import find_sequence_regex
from Bio.SeqIO import read, write
from pydna.parsers import parse


sf = SeqFeatureModel.model_validate({
    'type': 'CDS',
    'location': '1..2',
    'qualifiers': {
        'gene': ['ase1'],
        'product': ['Ase1']
    }
})

print(sf)
sff = sf.convert_to_seq_feature()
sff.strand

print(SeqFeatureModel.read_from_seq_feature(sff))

# print(sf.model_dump())

plasmid = read('examples/sequences/circular_dummy.gb', 'genbank')
print(plasmid.features)
write(plasmid, 'dummy.gb', 'genbank')
plasmid_pydna = parse('dummy.gb')[0]
print(plasmid_pydna.features)

# print(find_sequence_regex('taa', 'AAAATGAAT', False))
