import assembly2 as assembly
from Bio.SeqFeature import SeqFeature, SimpleLocation
from pydna.dseqrecord import Dseqrecord

f1 = Dseqrecord('aaaTTTctaGGGccc', circular=True)
f2 = Dseqrecord('ccccTTTatgGGGaaa')

f1_feat1 = SeqFeature(SimpleLocation(3, 6))
f1_feat2 = SeqFeature(SimpleLocation(9, 12))

f2_feat1 = SeqFeature(SimpleLocation(4, 7))
f2_feat2 = SeqFeature(SimpleLocation(10, 13))

f1.features = [f1_feat1, f1_feat2]
f2.features = [f2_feat1, f2_feat2]

def modify_feat(feat):
    if len(feat.location.parts) > 1:
        feat.location.parts = feat.location.parts[::-1]
    return feat

for shift in range(len(f1)):
    f1_shifted = f1.shifted(shift)

    list(map(modify_feat, f1_shifted.features))

    # Re-order the features so that TTT is first
    if str(f1_shifted.features[0].location.extract(f1_shifted.seq)) != 'TTT':
        f1_shifted.features = f1_shifted.features[::-1]

    # Linear assembly 2 - 1 - 2 (ccccTTTctaGGGaaa)
    assembly_plan = [
        (2, 1, f2.features[0].location, f1_shifted.features[0].location),
        (1, 2, f1_shifted.features[1].location, f2.features[1].location),
    ]

    result = assembly.assemble([f1_shifted, f2], assembly_plan, False).seq
    print(shift, str(result) == 'ccccTTTctaGGGaaa')

    # Circular assembly 1 - 2 (ccccTTTctaGGGaaa)
    assembly_plan = [
        (1, 2, f1_shifted.features[0].location, f2.features[0].location),
        (2, 1, f2.features[1].location, f1_shifted.features[1].location),
    ]

    result = assembly.assemble([f1_shifted, f2], assembly_plan, True).seq
    print(shift, str(result) == 'GGGcccaaaTTTatg')


