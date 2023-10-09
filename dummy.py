from pydantic_models import SeqFeatureModel

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