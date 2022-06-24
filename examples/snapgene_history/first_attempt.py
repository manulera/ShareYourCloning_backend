# %%
from typing import Optional
import xmltodict
from pydantic import BaseModel, Field
from pydantic_models import RestrictionEnzymeDigestionSource
from dna_functions import get_restriction_enzyme_products_list
# %%

"""
<User Id="2138">
    <FirstName />
    <LoggedIn Value="true" />
</User>
"""

xmlstring = """
<Node name="Cloned.dna" type="DNA" seqLen="6112" strandedness="double" ID="3" circular="1" operation="insertFragment">
    <InputSummary manipulation="replace" name1="AscI" name2="SalI" val1="63" val2="37" siteCount1="1" siteCount2="1" />
    <InputSummary manipulation="insert" name1="SalI" name2="AscI" val1="7" val2="2207" siteCount1="1" siteCount2="1" />
    <Node name="addgene-plasmid-39296-sequence-49545.dna" type="DNA" seqLen="3938" strandedness="double" ID="0" circular="1" useCustomMapLabel="1" customMapLabel="pFA6a-kanMX6" resurrectable="1" operation="invalid"></Node>
    <Node name="pcr_product.dna" type="DNA" seqLen="2219" strandedness="double" ID="2" circular="0" upstreamModification="Unmodified" downstreamModification="Unmodified" resurrectable="1" operation="amplifyFragment"></Node>
</Node>
"""


class InputSummary(BaseModel):
    manipulation: str
    name1: str
    name2: str
    val1: str
    val2: str
    siteCount1: str
    siteCount2: str


class Node(BaseModel):
    name: str
    type: str
    seqLen: str
    strandedness: str
    ID: str
    circular: str
    operation: str
    # Self-referencing fields have to be declared as strings!
    # https://pydantic-docs.helpmanual.io/usage/postponed_annotations/#self-referencing-models
    node: Optional[list['Node']] = Field([])
    # We cannot have camelCase in field names
    input_summary: Optional[list[InputSummary]] = Field([])


# user = Node.from_orm(fromstring(xmlstring))
# force_list makes all children be
# It seems that they either focus on the minimal start max end for the restriction
the_dict = xmltodict.parse(xmlstring, force_list=True)


def replace_at_symbols_in_dict(input_dict):
    dict_keys = list(input_dict.keys())
    for key in dict_keys:
        if key[0] == "@":
            input_dict[key[1:]] = input_dict[key]
            del input_dict[key]
    for key in dict_keys:
        # We cannot have camelCase in field names
        if key == 'InputSummary':
            input_dict['input_summary'] = input_dict[key]
            del input_dict[key]
        if key == 'Node':
            input_dict['node'] = input_dict[key]
            del input_dict[key]
    for key in input_dict:
        if type(input_dict[key]) == list:
            for d in input_dict[key]:
                replace_at_symbols_in_dict(d)


the_dict = the_dict['Node'][0]
replace_at_symbols_in_dict(the_dict)

node = Node.parse_obj(the_dict)

# It seems that each input summary corresponds to a Node (in the same order as the list, I guess?)
# For restriction and ligation (it seems):

# For now, let's assume 'replace' can only be a restriction-ligation

# It seems the edges of fragments are always the 5', but I wonder what would happen in a partial ligation
# Example: SalI gtcgac
# Cut: g^tcgac
# What happens if we would hybridate with a shorter one (only overlaping with gac). Then I guess they would just indicate the one furthest
# after cgc()
if node.input_summary[0].manipulation == 'replace':
    restriction_enzymes = [node.input_summary[0].name1, node.input_summary[0].name2]
    fragment_boundaries = list()

    # We look for the closest ones to the ones indicated in val1. This could be problematic if sites are overlaping, but that would be
    # an unlikely scenario. In any case, we will always check that the cloning works fine.
    snapgene_boundaries = [int(node.input_summary[0].val1), int(node.input_summary[0].val2)]

    source = RestrictionEnzymeDigestionSource(
        restriction_enzymes=[node.input_summary[0].name1, node.input_summary[0].name2],
        fragment_boundaries=[int(node.input_summary[0].val1), int(node.input_summary[0].val2)],
        input=[1]
    )

    # %%
