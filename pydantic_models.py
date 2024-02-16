from pydantic import BaseModel, Field, ConfigDict
from pydantic.types import constr, conlist
from enum import Enum
from typing import Optional
from Bio.SeqFeature import SeqFeature, Location
from Bio.SeqIO.InsdcIO import _insdc_location_string as format_feature_location
from Bio.Restriction.Restriction import RestrictionType
# Enumerations:

class SourceType(str, Enum):
    repository_id = 'repository_id',
    file = 'file',
    restriction = 'restriction'
    ligation = 'ligation'
    PCR = 'PCR'
    homologous_recombination = 'homologous_recombination'
    gibson_assembly = 'gibson_assembly'
    restriction_and_ligation = 'restriction_and_ligation'


class SequenceFileFormat(str, Enum):
    fasta = 'fasta'
    genbank = 'genbank'
    snapgene = 'snapgene'


class RepositoryName(str, Enum):
    genbank = 'genbank'
    addgene = 'addgene'

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
    id: Optional[int] = Field(None, description='Unique identifier of the sequence')
    kind: str = Field('entity', description='The kind entity (always equal to "entity"). Should probably be removed.')
    sequence: Optional[GenbankSequence] = Field(..., description='The sequence in genbank format + some extra info that is not captured by the genbank format')


class PrimerModel(BaseModel):
    """Called PrimerModel not to be confused with the class from pydna."""

    id: int
    name: str
    # TODO: add this to the flake8 exceptions
    # TODO: implement constrains when there is an answer for https://github.com/pydantic/pydantic/issues/7745
    sequence: constr(pattern='^[acgtACGT]+$')
    # sequence: str


class SeqFeatureModel(BaseModel):
    type: str
    qualifiers: dict[str, list[str]] = {}
    location: str

    def convert_to_seq_feature(self) -> SeqFeature:
        return SeqFeature(
            location=Location.fromstring(self.location),
            type=self.type,
            qualifiers=self.qualifiers
        )

    def read_from_seq_feature(sf: SeqFeature) -> 'SeqFeatureModel':
        return SeqFeatureModel(
            type=sf.type,
            qualifiers=sf.qualifiers,
            location=format_feature_location(sf.location, None)
        )

# Sources =========================================


class Source(BaseModel):
    """A class to represent sources of DNA
    """
    # Fields required to execute a source step
    id: Optional[int] = Field(None, description='Unique identifier of the source')
    kind: str = Field('source', description='The kind entity (always equal to "source"). Should probably be removed.')
    input: list[int] = Field([], description='Identifiers of the sequences that are an input to this source. \
                             If the source represents external import of a sequence, it\'s empty.')
    output: Optional[int] = Field(None, description='Identifier of the sequence that is an output of this source.')
    type: Optional[SourceType] = Field(..., description='The type source (PCR, restriction, etc.)')
    info: dict = Field({}, description='Additional information about the source (not used much yet, and probably should be removed)')
    model_config = ConfigDict(extra='forbid')


class UploadedFileSource(Source):
    """Describes a sequence from a file uploaded by the user
    """
    file_name: str
    file_format: SequenceFileFormat
    type: SourceType = SourceType('file')
    index_in_file: Optional[int] = None


class RepositoryIdSource(Source):
    """Documents a request to a repository
    """
    repository: RepositoryName
    repository_id: str
    type: SourceType = SourceType('repository_id')


class GenomeCoordinatesSource(Source):
    """Documents a request to NCBI for genome sequence
    """
    genome_accession: str
    locus_tag: Optional[str] = None
    start: int
    end: int
    strand: int
    type: SourceType = SourceType('genome_coordinates')


class SequenceCut(Source):
    """A class to represent a cut in a sequence"""

    left_edge : Optional[tuple[int, int]] = Field(None, description='The left edge of the cut, in the format (cut_watson, ovhg)')
    right_edge : Optional[tuple[int, int]] = Field(None, description='The right edge of the cut, in the format (cut_watson, ovhg)')

class RestrictionEnzymeDigestionSource(SequenceCut):
    """Documents a restriction enzyme digestion, and the selection of one of the fragments."""

    type: SourceType = SourceType('restriction')

    # The order of the enzymes in the list corresponds to the fragment_boundaries.
    # For instance, if a fragment 5' is cut with EcoRI and the 3' with BamHI,
    # restriction_enzymes = ['EcoRI', 'BamHI']
    restriction_enzymes: conlist(str|None, min_length=1) = Field(..., description='Enzymes associated with the left and right sides of the cut. It can contain None to represent the edge the sequence in linear sequences.')

    def from_cutsites(left: tuple[tuple[int,int], RestrictionType], right: tuple[tuple[int,int], RestrictionType], input: list[int], id: int) -> 'RestrictionEnzymeDigestionSource':
        return RestrictionEnzymeDigestionSource(
            restriction_enzymes=[None if left is None else str(left[1]), None if right is None else str(right[1])],
            left_edge=None if left is None else left[0],
            right_edge=None if right is None else right[0],
            input=input
        )

class AssemblySource(Source):
    assembly:  Optional[conlist(tuple[int, int, str, str], min_length=1)] = Field(None, description='The assembly plan as a list of tuples (part_1_id, part_2_id, loc1, loc2)')
    circular: Optional[bool] = Field(None, description='Whether the assembly is circular or not')

    def minimal_overlap(self):
        """Returns the minimal overlap between the fragments in the assembly"""
        all_overlaps = list()
        for f in self.assembly:
            if f[2] is not None:
                all_overlaps.append(len(Location.fromstring(f[2])))
            if f[3] is not None:
                all_overlaps.append(len(Location.fromstring(f[3])))
        return min(all_overlaps)

    def get_assembly_plan(self):
        """Returns the assembly plan"""
        out = list()
        for p in self.assembly:
            out.append((p[0], p[1], Location.fromstring(p[2]), Location.fromstring(p[3])))
        return tuple(out)

class PCRSource(AssemblySource):
    """Documents a PCR, and the selection of one of the products."""

    type: SourceType = SourceType('PCR')
    circular: bool = False
    forward_primer: int = Field(..., description='The forward primer')
    reverse_primer: int = Field(..., description='The reverse primer')

    # This can only take one input
    input: conlist(int, min_length=1, max_length=1)

    def from_assembly(assembly: list[tuple[int, int, Location, Location]], input: list[int], id: int, forward_primer: int, reverse_primer: int) -> 'PCRSource':
        """Creates a PCRSource from an assembly, input and id"""
        return PCRSource(
            id=id,
            assembly=[(part[0], part[1], format_feature_location(part[2], None), format_feature_location(part[3], None)) for part in assembly],
            input=input,
            forward_primer=forward_primer,
            reverse_primer=reverse_primer
        )

class LigationSource(AssemblySource):

    type: SourceType = SourceType('ligation')

    def from_assembly(assembly: list[tuple[int, int, Location, Location]], input: list[int], circular: bool, id: int) -> 'LigationSource':
        """Creates a StickyLigationSource from an assembly, input and circularity"""
        return LigationSource(
            id=id,
            assembly=[(part[0], part[1], format_feature_location(part[2], None), format_feature_location(part[3], None)) for part in assembly],
            input=input,
            circular=circular
        )


class HomologousRecombinationSource(AssemblySource):

    # This can only take two inputs, the first one is the template, the second one is the insert
    type: SourceType = SourceType('homologous_recombination')
    input: conlist(int, min_length=2, max_length=2)

    def from_assembly(assembly: list[tuple[int, int, Location, Location]], input: list[int], circular: bool, id: int) -> 'HomologousRecombinationSource':
        return HomologousRecombinationSource(
            id=id,
            assembly=[(part[0], part[1], format_feature_location(part[2], None), format_feature_location(part[3], None)) for part in assembly],
            input=input,
            circular=circular
        )

class GibsonAssemblySource(AssemblySource):

    type: SourceType = SourceType('gibson_assembly')
    input: conlist(int, min_length=1)

    def from_assembly(assembly: list[tuple[int, int, Location, Location]], input: list[int], circular: bool, id: int) -> 'GibsonAssemblySource':
        return GibsonAssemblySource(
            id=id,
            assembly=[(part[0], part[1], format_feature_location(part[2], None), format_feature_location(part[3], None)) for part in assembly],
            input=input,
            circular=circular
        )

class RestrictionAndLigationSource(AssemblySource):
    type: SourceType = SourceType('restriction_and_ligation')
    input: conlist(int, min_length=1)
    restriction_enzymes: conlist(str, min_length=1) = Field(..., description='The list of restriction enzymes used in the digestion')

    def from_assembly(assembly: list[tuple[int, int, Location, Location]], input: list[int], circular: bool, id: int, restriction_enzymes=list['str']) -> 'RestrictionAndLigationSource':
        return RestrictionAndLigationSource(
            id=id,
            assembly=[(part[0], part[1], format_feature_location(part[2], None), format_feature_location(part[3], None)) for part in assembly],
            input=input,
            circular=circular,
            restriction_enzymes=restriction_enzymes
        )
