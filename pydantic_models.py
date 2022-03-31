from pydantic import BaseModel, Field
from enum import Enum
from typing import Optional
from Bio.SeqFeature import SeqFeature, FeatureLocation

from pydantic.types import conlist


# Enumerations:

class SourceType(str, Enum):
    genbank_id = 'genbank_id',
    file = 'file',
    restriction = 'restriction'
    sticky_ligation = 'sticky_ligation'
    PCR = 'PCR'


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


class PrimerModel(BaseModel):
    """Called PrimerModel not to be confused with the class from pydna."""

    id: int
    name: str
    sequence: str


class PrimerAnnealing(BaseModel):
    """A class to represent the annealing of a primer. Minimal version of pydna's Primer."""

    primer_id: int = Field(..., description='The id of the primer that anneals.')
    position: int = Field(..., description='The position of the 3\' end of the\
        primer in the template sequence.')
    footprint_length: int = Field(..., description='The number of basepairs that are anealed\
        in the primer. Missmatch support should be added in the future.')


# The next two models are unused for now


class SequenceFeature(BaseModel):
    id: str
    type: str
    start: int
    end: int
    strand: int = None


def seq_feature2pydantic(sf: SeqFeature) -> SequenceFeature:
    if not isinstance(sf.location, FeatureLocation):
        raise TypeError(
            'Compound locations are not yet supported.'
        )
    return SequenceFeature(
        id=sf.id,
        type=sf.type,
        strand=sf.location.strand,
        start=sf.location.start,
        end=sf.location.end
    )

# Sources =========================================


class Source(BaseModel):
    """A class to represent sources of DNA
    """
    # Fields required to execute a source step
    id: int = None
    input: list[int] = []
    output: int = None
    type: SourceType = None
    output_index: int = None

    # Fields used to choose between multiple outputs
    # and other client-side functionality
    output_list: list[SequenceEntity] = []


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


# TODO There is some abstract common thing between restriction and PCR, since
# they select a subset of the molecule, perhaps they can be merged in some way.

class SequenceSubsetSource(Source):
    """An abstract class for sources that select a subset of a sequence, such as PCR and digestion."""

    # This can only take one input
    input: conlist(int, min_items=1, max_items=1)

    # Boundaries of a fragment (length should be either empty, or length = 2)
    fragment_boundaries: list[int] = []


class RestrictionEnzymeDigestionSource(SequenceSubsetSource):
    """Documents a restriction enzyme digestion, and the selection of one of the fragments."""

    type: SourceType = SourceType('restriction')

    # The order of the enzymes in the list corresponds to the fragment_boundaries.
    # For instance, if a fragment 5' is cut with EcoRI and the 3' with BamHI,
    # restriction_enzymes = ['EcoRI', 'BamHI']
    restriction_enzymes: conlist(str, min_items=1)


class PrimerAnnealingSettings(BaseModel):
    """Settings to find annealing sites for the primer"""
    minimum_annealing: int = Field(..., description='The minimum number of \
    overlaping basepairs for an annealing to be considered.')


class PrimerPair(BaseModel):
    forward: PrimerAnnealing
    reverse: PrimerAnnealing


class PCRSource(Source):
    """Documents a PCR, and the selection of one of the products."""

    type: SourceType = SourceType('PCR')

    # This can only take one input
    input: conlist(int, min_items=1, max_items=1)

    primer_pair: PrimerPair = None

    # TODO test this
    primer_annealing_settings: PrimerAnnealingSettings = Field(None, description='This does not have\
        to be specified if the primer annealing fields are provided.')

    # Field for client-side functionality
    possible_primer_pairs: list[PrimerPair] = []


class StickyLigationSource(Source):
    """Documents a ligation with sticky ends. This might consist of \
    a single fragment's circularisation"""

    # TODO: this should support at some point specifying the order of the fragments
    # of the assembly + whether there is circularization.
    input: conlist(int, min_items=1)
    type: SourceType = SourceType('sticky_ligation')
    fragments_inverted: list[bool] = []
    circularised: bool = None

    # TODO include this
    # @validator('fragments_inverted')
    # def lists_have_equal_length(cls, v, values):
    #     assert len(v) == len(values['input']) or len(v) == 0, '`fragments_inverted` must\
    #         be either empty, or have the same length as `input`'
