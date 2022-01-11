from pydantic import BaseModel, Field
from enum import Enum
from typing import List, Optional

from pydantic.types import conlist


# Enumerations:

class SourceType(str, Enum):
    genbank_id = 'genbank_id',
    file = 'file',
    restriction = 'restriction'
    sticky_ligation = 'sticky_ligation'


class SequenceFileFormat(str, Enum):
    fasta = 'fasta'
    genbank = 'genbank'
    snapgene = 'snapgene'


# Sequence: =========================================


class GenbankSequence(BaseModel):
    """A class to store sequences and features in genbank model
    """
    type: str = 'file'
    file_extension: str = 'gb'
    file_content: str = ''
    overhang_crick_3prime: int = Field(0, description='Taken from pydna\'s `dseq::ovhg`\
        An integer describing the length of the\
        crick strand overhang in the 5\' of the molecule, or 3\' of the crick strand')
    overhang_watson_3prime: int = Field(0, description='The equivalent of `overhang_crick_3prime`\
        but for the watson strand')


class SequenceEntity(BaseModel):
    id: Optional[int]
    kind: str = 'entity'
    sequence: GenbankSequence = None

# Sources =========================================


class Source(BaseModel):
    """A class to represent sources of DNA
    """
    # Fields required to execute a source step
    id: int = None
    input: List[int] = []
    output: int = None
    type: SourceType = None
    output_index: int = None

    # Fields used to choose between multiple outputs
    # and other client-side functionality
    output_list: List[SequenceEntity] = []


class UploadedFileSource(Source):
    """Describes a sequence from a file uploaded by the user
    """
    file_name: str
    file_format: SequenceFileFormat
    type: SourceType = SourceType('file')


class GenbankIdSource(Source):
    """Documents a request to GenBank
    """
    genbank_id: str
    type: SourceType = SourceType('genbank_id')


class RestrictionEnzymeDigestionSource(Source):
    """Documents a restriction enzyme digestion, and the selection of
    one of the fragments
    """

    type: SourceType = SourceType('restriction')
    # This can only take one input
    input: conlist(int, min_items=1, max_items=1)

    # Field required to execute the source step
    restriction_enzymes: conlist(str, min_items=1)

    # Field for client-side functionality
    fragment_boundaries: List[int] = []


class StickyLigationSource(Source):
    """Documents a ligation with sticky ends. This might consist of \
    a single fragment's circularisation"""

    # TODO: this should support at some point specifying the order of the fragments
    # of the assembly + whether there is circularization.
    input: conlist(int, min_items=1)
    type: SourceType = SourceType('sticky_ligation')
    fragments_inverted: List[bool] = []
    circularised: bool = None

    # TODO include this
    # @validator('fragments_inverted')
    # def lists_have_equal_length(cls, v, values):
    #     assert len(v) == len(values['input']) or len(v) == 0, '`fragments_inverted` must\
    #         be either empty, or have the same length as `input`'
