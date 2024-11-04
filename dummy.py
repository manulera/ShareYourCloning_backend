from pydna.dseqrecord import Dseqrecord as _Dseqrecord
from Bio.Seq import reverse_complement
from Bio.Data.IUPACData import ambiguous_dna_values as _ambiguous_dna_values
import re
import itertools as _itertools
import assembly2 as assembly
from dna_functions import custom_file_parser
import glob

ambiguous_only_dna_values = {**_ambiguous_dna_values}
for normal_base in 'ACGT':
    del ambiguous_only_dna_values[normal_base]


def compute_regex_site(site: str) -> str:
    upper_site = site.upper()
    for k, v in ambiguous_only_dna_values.items():
        if len(v) > 1:
            upper_site = upper_site.replace(k, f"[{''.join(v)}]")

    # Make case insensitive
    upper_site = f'(?i){upper_site}'
    return upper_site


def dseqrecord_finditer(pattern: str, seq: _Dseqrecord) -> list[re.Match]:
    query = str(seq.seq) if not seq.circular else str(seq.seq) * 2
    matches = re.finditer(pattern, query)
    return (m for m in matches if m.start() <= len(seq))


# Taken from gateway_features.txt, shipped with ApE
raw_gateway_sites = {
    'attB4': 'CMASTwTGTATAGAAAAGYWG',
    'attB3': 'CMASTwTGTATAATAAAGYWG',
    'attB2': 'CMASTwTGTACAAGAAAGYWG',
    'attB1': 'CMASTwTGTACAAAAAAGYWG',
    'attP4': 'AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTwTGTATAGAAAAGYWGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATSCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTA',
    'attP3': 'AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTwTGTATAATAAAGYWGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATSCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTA',
    'attP2': 'AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTwTGTACAAGAAAGYWGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATSCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTA',
    'attP1': 'AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTwTGTACAAAAAAGYWGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATSCAGTCACTATGAATCAACTACTTAGATGGTATTAGTGACCTGTA',
    'attL4': 'AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTwTGTATAGAAAAGYWG',
    'attL3': 'AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTwTGTATAATAAAGYWG',
    'attL2': 'AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTwTGTACAAGAAAGYWG',
    'attL1': 'AAATAATGATTTTATTTTGACTGATAGTGACCTGTTCGTTGCAACAMATTGATRAGCAATGCTTTYTTATAATGCCMASTwTGTACAAAAAAGYWG',
    # In the attR ones, I have removed a trailing AATCAACTACTTAGATGGTATTAGTGACCTGTA from the original file
    'attR4': 'CMASTwTGTATAGAAAAGYWGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATSCAGTCACTATG',
    'attR3': 'CMASTwTGTATAATAAAGYWGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATSCAGTCACTATG',
    'attR2': 'CMASTwTGTACAAGAAAGYWGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATSCAGTCACTATG',
    'attR1': 'CMASTwTGTACAAAAAAGYWGAACGAGAAACGTAAAATGATATAAATATCAATATATTAAATTAGATTTTGCATAAAAAACAGACTACATAATACTGTAAAACACAACATATSCAGTCACTATG',
    'ccdB': 'ATGGTGATCCCCCTGGCCAGTGCACGTCTGCTGTCAGATAAAGTCNCCCGTGAACTTTACCCGGTGGTGCATATCGGGGATGAAAGCTGGCGCATGATGACCACCGATATGGCCAGTGTGCCGGTCTCCGTTATCGGGGAAGAAGTGGCTGATCTCAGCCACCGCGAAAATGACATCAAAAACGCCATTAACCTGATGTTCTGGGGAATA',
    'core_4': 'cmastwtGTATAGAaaagywg',
    'core_3': 'cmastwtGTATAATaaagywg',
    'core_2': 'cmastwtGTACAAGaaagywg',
    'core_1': 'cmastwtGTACAAAaaagywg',
    'overlap_4': 'twtGTATAGAaaag',
    'overlap_3': 'twtGTATAATaaag',
    'overlap_2': 'twtGTACAAGaaag',
    'overlap_1': 'twtGTACAAAaaag',
}

gateway_sites = {
    k: {
        'forward_regex': compute_regex_site(v),
        'reverse_regex': compute_regex_site(reverse_complement(v)),
        'consensus_sequence': v,
    }
    for k, v in raw_gateway_sites.items()
}


def gateway_overlap(seqx: _Dseqrecord, seqy: _Dseqrecord, type=None):
    """Find gateway overlaps"""
    if type not in ['BP', 'LR']:
        raise ValueError(f'Invalid overlap type: {type}')

    out = list()
    # Iterate over the four possible att sites
    for num in range(1, 5):
        # Iterate over the two possible orientations
        # The sites have to be in the same orientation (fwd + fwd or rev + rev)
        for pattern in ['forward_regex', 'reverse_regex']:
            # The overlap regex is the same for all types
            overlap_regex = gateway_sites[f'overlap_{num}'][pattern]

            # Iterate over pairs B, P and P, B for BP and L, R and R, L for LR
            for site_x, site_y in zip(type, type[::-1]):
                site_x_regex = gateway_sites[f'att{site_x}{num}'][pattern]
                matches_x = list(dseqrecord_finditer(site_x_regex, seqx))
                if len(matches_x) == 0:
                    continue

                site_y_regex = gateway_sites[f'att{site_y}{num}'][pattern]
                matches_y = list(dseqrecord_finditer(site_y_regex, seqy))
                if len(matches_y) == 0:
                    continue

                for match_x, match_y in _itertools.product(matches_x, matches_y):
                    # Find the overlap sequence within each match
                    overlap_x = re.search(overlap_regex, match_x.group())
                    overlap_y = re.search(overlap_regex, match_y.group())

                    # Sanity check
                    assert (
                        overlap_x is not None and overlap_y is not None
                    ), 'Something went wrong, no overlap found within the matches'

                    out.append(
                        (
                            match_x.start() + overlap_x.start(),
                            match_y.start() + overlap_y.start(),
                            len(overlap_x.group()),
                        )
                    )

    return out


with open('test_files/gateway_manual_cloning/pairing.tsv') as f:
    for line in f:
        line = line.strip().split('\t')
        if len(line) < 5:
            backbone, pcr_product, entry_vector = line
            backbone_expression = None
            expression = None
        else:
            backbone, pcr_product, entry_vector, backbone_expression, expression = line

        print(backbone, pcr_product)
        with open('test_files/gateway_manual_cloning/' + backbone, 'rb') as f:
            backbone = custom_file_parser(f, 'snapgene')[0]
        with open('test_files/gateway_manual_cloning/' + pcr_product, 'rb') as f:
            pcr_product = custom_file_parser(f, 'snapgene')[0]
        with open('test_files/gateway_manual_cloning/' + entry_vector, 'rb') as f:
            entry_vector = custom_file_parser(f, 'snapgene')[0]

        def algo(x, y, _):
            return gateway_overlap(x, y, 'BP')

        asm = assembly.Assembly([backbone, pcr_product], algorithm=algo)
        # print(*asm.G.edges, sep='\n')
        out = asm.assemble_circular()
        if len(out):
            seq = out[0]
            asem = asm.get_circular_assemblies()[0]
            print(assembly.assembly2str(asem))
            if seq.seguid() == entry_vector.seguid():
                print('OK')
            else:
                print('NOT OK')
        else:
            print('NO OVERLAP')

        if backbone_expression is not None and expression is not None and len(out):
            with open('test_files/gateway_manual_cloning/' + backbone_expression, 'rb') as f:
                backbone_expression = custom_file_parser(f, 'snapgene')[0]
            with open('test_files/gateway_manual_cloning/' + expression, 'rb') as f:
                expression = custom_file_parser(f, 'snapgene')[0]

            def algo(x, y, _):
                return gateway_overlap(x, y, 'LR')

            asm = assembly.Assembly([backbone_expression, entry_vector], algorithm=algo)

            out = asm.assemble_circular()
            if len(out):
                seq = out[0]
                asem = asm.get_circular_assemblies()[0]
                print(assembly.assembly2str(asem))
                if seq.seguid() == expression.seguid():
                    print('OK')
                else:
                    print('NOT OK')
            else:
                print('NO OVERLAP')


example_valerie = glob.glob('test_files/gateway_manual_cloning/example_valerie/*.dna')
inputs = list()
for file in example_valerie:
    with open(file, 'rb') as f:
        print(file)
        seq = custom_file_parser(f, 'snapgene')[0]
        seq.name = file.split('/')[-1]
        inputs.append(seq)


def algo(x, y, _):
    return gateway_overlap(x, y, 'LR')


asm = assembly.Assembly(inputs, algorithm=algo, use_all_fragments=True, use_fragment_order=False)

out = asm.assemble_circular()
print(len(out))

out[0].write('test_files/gateway_manual_cloning/example_valerie/output.gb', 'gb')
out[1].write('test_files/gateway_manual_cloning/example_valerie/output2.gb', 'gb')
