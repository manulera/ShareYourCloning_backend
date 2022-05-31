# %%

import xmltodict
from xml.etree.ElementTree import fromstring
from pydantic.utils import GetterDict
from typing import Any, List, Optional
import lxml.objectify
import sys
from pydantic import BaseModel, Field
from typing import List
# We can use models for type-hinting
# class Node(BaseModel):


# Read file to object:
with open('snapgene_history_minimal.xml') as ins:
    whole_file = ins.read()

xml = lxml.objectify.fromstring(str(whole_file))

# We start from the newest molecule, and we have to build up the history from that
# This means that the length of xml.Node should always be 1

if len(xml.Node) != 1:
    sys.exit('there is more than one initial node')


# def appendHistory(parent_node):
#     # print(type(parent_node))

#     # for child_node in parent_node.Node:
#     #     print(type(child_node))
#     kind = parent_node.


xml.Node.values()

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
    Node: Optional[list[Node]]
    InputSummary: Optional[list[InputSummary]]


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

    for key in input_dict:
        if type(input_dict[key]) == list:
            for d in input_dict[key]:
                replace_at_symbols_in_dict(d)


the_dict = the_dict['Node'][0]
replace_at_symbols_in_dict(the_dict)

first_node = Node.parse_obj(the_dict)

first_node.InputSummary
