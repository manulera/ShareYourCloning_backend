from dna_functions import find_sequence_regex, location_edges, get_homologous_recombination_locations, perform_homologous_recombination
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import SimpleLocation, Location, CompoundLocation, ExactPosition
from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
import regex

# compiled_pattern = regex.compile('AA.*AA', regex.IGNORECASE)
# print(list(regex.finditer(compiled_pattern, 'GGGAAAACCC', overlapped=True)))
# exit()
# # print(find_sequence_regex('taa', 'AAAATGAAT', True))

locations = find_sequence_regex('att', 'AAAATGAAT', True)

for loc in locations:
    feat = SeqFeature(location=loc, type='CDS')
    # print(feat.extract('AAAATGAAT'), loc)
    # print(loc.start, loc.end)
    # If location is a SimpleLocation
    print('edges', location_edges(loc))


# loc = SimpleLocation(-3,-1, strand=1)
# feat = SeqFeature(location=locations[0], type='CDS')
# print(feat.extract('AAAATGAAT'))

dseq = Dseq('GGGAAAACCC', circular=True)
locations = find_sequence_regex('AAA', str(dseq), True)
print(locations)
edges = location_edges(locations[0])

# If not circular
fragment_1 = dseq[0:edges[0]]
fragment_2 = dseq[edges[1]:]

print(fragment_1)
print(fragment_2)

# If circular

fragment = dseq[edges[1]:edges[0]]
print(fragment)
print(fragment.circular)

template = Dseqrecord(dseq)
insert = Dseqrecord('AATTCCAA')


locs = get_homologous_recombination_locations(template, insert, 2)

print(perform_homologous_recombination(template, insert, locs[0]))

