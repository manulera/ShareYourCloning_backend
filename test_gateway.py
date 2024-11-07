from gateway import gateway_overlap
import assembly2 as assembly
import glob
from dna_functions import custom_file_parser


def algoBP(x, y, _):
    return gateway_overlap(x, y, 'BP')


def algoLR(x, y, _):
    return gateway_overlap(x, y, 'LR')


def test_gateway_manual_cloning():

    with open('test_files/gateway_manual_cloning/pairing.tsv') as f:
        for line in f:
            line = line.strip().split('\t')
            if len(line) < 5:
                backbone, pcr_product, entry_vector = line
                backbone_expression = None
                expression = None
            else:
                backbone, pcr_product, entry_vector, backbone_expression, expression = line

            with open('test_files/gateway_manual_cloning/' + backbone, 'rb') as f:
                backbone = custom_file_parser(f, 'snapgene')[0]
            with open('test_files/gateway_manual_cloning/' + pcr_product, 'rb') as f:
                pcr_product = custom_file_parser(f, 'snapgene')[0]
            with open('test_files/gateway_manual_cloning/' + entry_vector, 'rb') as f:
                entry_vector = custom_file_parser(f, 'snapgene')[0]

            # Works with the right reaction
            asm = assembly.Assembly([backbone, pcr_product], algorithm=algoBP)
            out = asm.assemble_circular()
            seguids = [seq.seguid() for seq in out]
            assert entry_vector.seguid() in seguids

            # Does not work with the wrong reaction
            asm = assembly.Assembly([backbone, pcr_product], algorithm=algoLR)
            out = asm.assemble_circular()
            assert len(out) == 0

            if backbone_expression is not None and expression is not None and len(out):
                with open('test_files/gateway_manual_cloning/' + backbone_expression, 'rb') as f:
                    backbone_expression = custom_file_parser(f, 'snapgene')[0]
                with open('test_files/gateway_manual_cloning/' + expression, 'rb') as f:
                    expression = custom_file_parser(f, 'snapgene')[0]

                # Works with the right reaction
                asm = assembly.Assembly([backbone_expression, entry_vector], algorithm=algoLR)

                out = asm.assemble_circular()
                seguids = [seq.seguid() for seq in out]
                assert expression.seguid() in seguids

                # Does not work with the wrong reaction
                asm = assembly.Assembly([backbone_expression, entry_vector], algorithm=algoBP)

                out = asm.assemble_circular()
                assert len(out) == 0

    example_valerie = glob.glob('test_files/gateway_manual_cloning/example_valerie/*.dna')
    inputs = list()
    for file in example_valerie:
        with open(file, 'rb') as f:
            seq = custom_file_parser(f, 'snapgene')[0]
            seq.name = file.split('/')[-1]
            inputs.append(seq)

    def algo(x, y, _):
        return gateway_overlap(x, y, 'LR')

    asm = assembly.Assembly(inputs, algorithm=algo, use_all_fragments=True, use_fragment_order=False)

    out = asm.assemble_circular()
