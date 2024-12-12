import os
import json
from Bio.Seq import reverse_complement
from fastapi.responses import FileResponse
from typing import Annotated
from fastapi import Body, Query, HTTPException

from ...get_router import get_router
from ...pydantic_models import (
    BaseCloningStrategy,
    PrimerModel,
    TextFileSequence,
    PCRSource,
    RestrictionAndLigationSource,
    GatewaySource,
)
from ...endpoints.assembly import pcr, restriction_and_ligation, gateway


router = get_router()


def validate_protospacers(protospacers):
    """
    Ensure that the protospacers are compatible for the golden gate assembly, by checking that the central 4
    bps are different in all of them, and also different from the other assembly joints.
    """
    forbidden_joints = ['AACA', 'GCTT']

    for i, ps in enumerate(protospacers):
        # Must be only ACGT
        if not set(ps).issubset({'A', 'C', 'G', 'T'}):
            raise ValueError(f"Protospacer {i} {ps} contains invalid bases")
        if len(ps) != 20:
            raise ValueError(f"Protospacer {i} {ps} is not 20 bp long")
        if ps[8:12] in forbidden_joints:
            raise ValueError(
                f"Protospacer {i} has a forbidden joint {ps[8:12]}, which is used for the constant parts of the assembly"
            )
        # Find any other protospacers with the same joint sequence
        same_joint_indexes = [
            j + 1 for j, other_ps in enumerate(protospacers) if i != j and ps[8:12] == other_ps[8:12]
        ]
        if same_joint_indexes:
            raise ValueError(
                f"Protospacer {i + 1} has the same joint as protospacers {', '.join(map(str, same_joint_indexes))}"
            )


def design_primers(protospacers):

    fwd_prefix = 'taggtctcc'
    rvs_prefix = 'atggtctca'
    fwd_suffix = 'gttttagagctagaa'
    rvs_suffix = 'tgcaccagccgggaa'

    primers = list()
    for protospacer in protospacers:
        primers.append(rvs_prefix + reverse_complement(protospacer[:12]) + rvs_suffix)
        primers.append(fwd_prefix + protospacer[8:] + fwd_suffix)
    return primers


@router.get('/batch_cloning/ziqiang_et_al2024')
def ziqiang_et_al2024():
    return FileResponse(os.path.join(os.path.dirname(__file__), 'index.html'))


@router.post('/batch_cloning/ziqiang_et_al2024', response_model=BaseCloningStrategy)
async def ziqiang_et_al2024_post(
    protospacers: Annotated[list[str], Body(..., min_length=1)], until_bp: bool = Query(False)
):
    try:
        validate_protospacers(protospacers)
    except ValueError as e:
        raise HTTPException(400, str(e))
    primers = design_primers(protospacers)

    with open(os.path.join(os.path.dirname(__file__), 'ziqiang_et_al2024.json'), 'r') as f:
        template = BaseCloningStrategy.model_validate(json.load(f))

    max_primer_id = max([primer.id for primer in template.primers], default=0)

    for i, primer in enumerate(primers):
        max_primer_id += 1
        orientation = 'rvs' if i % 2 == 0 else 'fwd'
        template.primers.append(
            PrimerModel(id=max_primer_id, name=f"protospacer_{i // 2 + 1}_{orientation}", sequence=primer)
        )

    primer_ids_for_pcrs = [3, *[p.id for p in template.primers[-len(primers) :]], 12]
    next_node_id = max([s.id for s in template.sequences] + [s.id for s in template.sources]) + 1

    template_sequence = next(s for s in template.sequences if s.id == 18)
    for i, (fwd_primer_id, rvs_primer_id) in enumerate(zip(primer_ids_for_pcrs[::2], primer_ids_for_pcrs[1::2])):
        if i == 0:
            name = 'start_ps1'
        elif i == (len(primer_ids_for_pcrs) // 2) - 1:
            name = f'end_ps{i}'
        else:
            name = f'end_ps{i}_start_ps{i + 1}'

        pcr_source = PCRSource(id=next_node_id, output_name=name)
        fwd_primer = next(p for p in template.primers if p.id == fwd_primer_id)
        rvs_primer = next(p for p in template.primers if p.id == rvs_primer_id)

        next_node_id += 1
        resp = await pcr(pcr_source, [template_sequence], [fwd_primer, rvs_primer], 14, 0)
        pcr_product: TextFileSequence = TextFileSequence.model_validate(resp['sequences'][0])
        pcr_product.id = next_node_id
        pcr_source: PCRSource = PCRSource.model_validate(resp['sources'][0])
        pcr_source.output = next_node_id

        template.sequences.append(pcr_product)
        template.sources.append(pcr_source)

        next_node_id += 1

    # Find all PCR products
    # (we use type instead of isinstance because the BaseCloningStrategy does not
    #  have the newer source models with extra methods)
    pcr_product_ids = [s.output for s in template.sources if s.type == 'PCRSource']

    # Make all input of a Golden gate assembly
    golden_gate_source = RestrictionAndLigationSource(
        id=next_node_id, output_name='golden_gate_assembly', restriction_enzymes=['BsaI'], input=pcr_product_ids
    )

    next_node_id += 1
    # Make them
    input_sequences = [next(s for s in template.sequences if s.id == p) for p in pcr_product_ids]
    resp = await restriction_and_ligation(golden_gate_source, input_sequences, False, False)
    golden_gate_product: TextFileSequence = TextFileSequence.model_validate(resp['sequences'][0])
    golden_gate_product.id = next_node_id
    golden_gate_source: RestrictionAndLigationSource = RestrictionAndLigationSource.model_validate(resp['sources'][0])
    golden_gate_source.output = next_node_id
    next_node_id += 1

    template.sequences.append(golden_gate_product)
    template.sources.append(golden_gate_source)

    bp_target = next(s for s in template.sequences if s.id == 12)
    gateway_source = GatewaySource(id=next_node_id, output_name='entry_clone', reaction_type='BP', greedy=False)
    next_node_id += 1
    resp = await gateway(gateway_source, [golden_gate_product, bp_target], circular_only=True, only_multi_site=True)
    gateway_product: TextFileSequence = TextFileSequence.model_validate(resp['sequences'][0])
    gateway_product.id = next_node_id
    gateway_source: GatewaySource = GatewaySource.model_validate(resp['sources'][0])
    gateway_source.output = next_node_id
    next_node_id += 1

    template.sequences.append(gateway_product)
    template.sources.append(gateway_source)

    if until_bp:
        # Delete sources and sequences left
        ids2delete = list(range(5, 11))
        template.sources = [s for s in template.sources if s.id not in ids2delete]
        template.sequences = [s for s in template.sequences if s.id not in ids2delete]
        return template

    # Now we want to do a Gateway with everything, so we need to find all sequences that are not input of anything
    all_input_ids = sum([s.input for s in template.sources], [])
    sequences_to_clone = [s for s in template.sequences if s.id not in all_input_ids]

    gateway_source = GatewaySource(id=next_node_id, output_name='expression_clone', reaction_type='LR', greedy=False)
    next_node_id += 1
    resp = await gateway(gateway_source, sequences_to_clone, circular_only=True, only_multi_site=True)
    index_of_product = next(i for i, s in enumerate(resp['sequences']) if '/label="Cas9"' in s.file_content)
    expression_clone: TextFileSequence = TextFileSequence.model_validate(resp['sequences'][index_of_product])
    expression_clone.id = next_node_id
    gateway_source: GatewaySource = GatewaySource.model_validate(resp['sources'][index_of_product])
    gateway_source.output = next_node_id
    next_node_id += 1

    template.sequences.append(expression_clone)
    template.sources.append(gateway_source)

    return template
