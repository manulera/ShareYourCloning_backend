from pydantic import BaseModel, Field
from enum import Enum
from typing import List, Optional

from pydantic.types import conlist


# Enumerations:

class SourceType(str, Enum):
    genbank_id = 'genbank_id',
    file = 'file',
    restriction = 'restriction'


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
    id: int = None
    input: List[int] = []
    output_list: List[SequenceEntity] = []
    output_index: int = None
    output: int = None
    type: SourceType = None


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
    restriction_enzymes: conlist(str, min_items=1)
    type: SourceType = SourceType('restriction')
    fragment_boundaries: List[int] = []
