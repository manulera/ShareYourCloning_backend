from pydantic import BaseModel, Field, ConfigDict, model_validator
from enum import Enum
from typing import Optional, List
from Bio.SeqFeature import SeqFeature, Location, SimpleLocation as BioSimpleLocation
from Bio.SeqIO.InsdcIO import _insdc_location_string as format_feature_location
from Bio.Restriction.Restriction import RestrictionType, RestrictionBatch
from shareyourcloning_linkml.datamodel import (
    OligoHybridizationSource as _OligoHybridizationSource,
    PolymeraseExtensionSource as _PolymeraseExtensionSource,
    GenomeCoordinatesSource as _GenomeCoordinatesSource,
    RepositoryIdSource as _RepositoryIdSource,
    ManuallyTypedSource as _ManuallyTypedSource,
    UploadedFileSource as _UploadedFileSource,
    SequenceFileFormat as _SequenceFileFormat,
    RestrictionEnzymeDigestionSource as _RestrictionEnzymeDigestionSource,
    RestrictionSequenceCut as _RestrictionSequenceCut,
    TextFileSequence as _TextFileSequence,
    AssemblySource as _AssemblySource,
    PCRSource as _PCRSource,
    HomologousRecombinationSource as _HomologousRecombinationSource,
    GibsonAssemblySource as _GibsonAssemblySource,
    RestrictionAndLigationSource as _RestrictionAndLigationSource,
    LigationSource as _LigationSource,
    CRISPRSource as _CRISPRSource,
    Primer as _Primer,
    AssemblyJoin as _AssemblyJoin,
    AssemblyJoinComponent as _AssemblyJoinComponent,
    SimpleSequenceLocation as _SimpleSequenceLocation,
    AddGeneIdSource as _AddGeneIdSource,
)


SequenceFileFormat = _SequenceFileFormat


class SourceType(str, Enum):
    ligation = 'ligation'
    PCR = 'PCR'
    homologous_recombination = 'homologous_recombination'
    crispr = 'crispr'
    gibson_assembly = 'gibson_assembly'
    restriction_and_ligation = 'restriction_and_ligation'


class TextFileSequence(_TextFileSequence):
    pass


class PrimerModel(_Primer):
    """Called PrimerModel not to be confused with the class from pydna."""


class SeqFeatureModel(BaseModel):
    type: str
    qualifiers: dict[str, list[str]] = {}
    location: str

    def convert_to_seq_feature(self) -> SeqFeature:
        return SeqFeature(location=Location.fromstring(self.location), type=self.type, qualifiers=self.qualifiers)

    def read_from_seq_feature(sf: SeqFeature) -> 'SeqFeatureModel':
        return SeqFeatureModel(
            type=sf.type, qualifiers=sf.qualifiers, location=format_feature_location(sf.location, None)
        )


# Sources =========================================


class Source(BaseModel):
    """A class to represent sources of DNA"""

    # Fields required to execute a source step
    id: Optional[int] = Field(None, description='Unique identifier of the source')
    kind: str = Field('source', description='The kind entity (always equal to "source"). Should probably be removed.')
    input: list[int] = Field(
        [],
        description='Identifiers of the sequences that are an input to this source. \
                             If the source represents external import of a sequence, it\'s empty.',
    )
    output: Optional[int] = Field(None, description='Identifier of the sequence that is an output of this source.')
    type: Optional[SourceType] = Field(..., description='The type source (PCR, restriction, etc.)')
    info: dict = Field(
        {}, description='Additional information about the source (not used much yet, and probably should be removed)'
    )
    model_config = ConfigDict(extra='forbid')


class ManuallyTypedSource(_ManuallyTypedSource):
    """Describes a sequence that is typed manually by the user"""

    @model_validator(mode='after')
    def validate_circularity(self):
        # Do the validation instead of printing
        if self.circular:
            assert self.overhang_crick_3prime == 0, 'Circular sequences cannot have overhangs.'
            assert self.overhang_watson_3prime == 0, 'Circular sequences cannot have overhangs.'
        return self


class UploadedFileSource(_UploadedFileSource):
    pass


class RepositoryIdSource(_RepositoryIdSource):
    pass


class AddGeneIdSource(_AddGeneIdSource):
    # TODO: add this to LinkML
    # repository_name: RepositoryName = RepositoryName('addgene')
    pass


class GenomeCoordinatesSource(_GenomeCoordinatesSource):
    pass


class RestrictionSequenceCut(_RestrictionSequenceCut):

    @classmethod
    def from_cutsite_tuple(cls, cutsite_tuple: tuple[tuple[int, int], RestrictionType]):
        cut_watson, ovhg = cutsite_tuple[0]
        enzyme = str(cutsite_tuple[1])

        return cls(
            cut_watson=cut_watson,
            overhang=ovhg,
            restriction_enzyme=enzyme,
        )

    def to_cutsite_tuple(self) -> tuple[tuple[int, int], RestrictionType]:
        restriction_enzyme = RestrictionBatch(first=[self.restriction_enzyme]).pop()
        return ((self.cut_watson, self.overhang), restriction_enzyme)


class RestrictionEnzymeDigestionSource(_RestrictionEnzymeDigestionSource):
    """Documents a restriction enzyme digestion, and the selection of one of the fragments."""

    # TODO: maybe a better way? They have to be redefined here because
    # we have overriden the original class

    left_edge: Optional[RestrictionSequenceCut] = Field(None)
    right_edge: Optional[RestrictionSequenceCut] = Field(None)

    @classmethod
    def from_cutsites(
        cls,
        left: tuple[tuple[int, int], RestrictionType],
        right: tuple[tuple[int, int], RestrictionType],
        input: list[int],
        id: int,
    ):
        return cls(
            id=id,
            left_edge=None if left is None else RestrictionSequenceCut.from_cutsite_tuple(left),
            right_edge=None if right is None else RestrictionSequenceCut.from_cutsite_tuple(right),
            input=input,
        )

    # TODO could be made into a computed field?
    def get_enzymes(self) -> list[str]:
        """Returns the enzymes used in the digestion"""
        out = list()
        if self.left_edge is not None:
            out.append(self.left_edge.restriction_enzyme)
        if self.right_edge is not None:
            out.append(self.right_edge.restriction_enzyme)
        # Unique values, sorted the same way
        return sorted(list(set(out)), key=out.index)


class SimpleSequenceLocation(_SimpleSequenceLocation):
    @classmethod
    def from_simple_location(cls, location: BioSimpleLocation):
        return cls(
            start=location.start,
            end=location.end,
            strand=location.strand,
        )

    def to_simple_location(self) -> BioSimpleLocation:
        return BioSimpleLocation(self.start, self.end, self.strand)


class AssemblyJoinComponent(_AssemblyJoinComponent):
    location: SimpleSequenceLocation = Field(
        ...,
        description="""Location of the overlap in the fragment. Might be an empty location (start == end) to indicate blunt join.""",
    )
    pass


class AssemblyJoin(_AssemblyJoin):

    left: AssemblyJoinComponent = Field(...)
    right: AssemblyJoinComponent = Field(...)

    @classmethod
    def from_join_tuple(cls, join_tuple: tuple[int, int, BioSimpleLocation, BioSimpleLocation]):
        return cls(
            left=AssemblyJoinComponent(
                sequence=abs(join_tuple[0]),
                reverse_complemented=join_tuple[0] < 0,
                location=SimpleSequenceLocation.from_simple_location(join_tuple[2]),
            ),
            right=AssemblyJoinComponent(
                sequence=abs(join_tuple[1]),
                reverse_complemented=join_tuple[1] < 0,
                location=SimpleSequenceLocation.from_simple_location(join_tuple[3]),
            ),
        )

    def to_join_tuple(self) -> tuple[int, int, BioSimpleLocation, BioSimpleLocation]:
        return (
            self.left.sequence * (-1 if self.left.reverse_complemented else 1),
            self.right.sequence * (-1 if self.right.reverse_complemented else 1),
            self.left.location.to_simple_location(),
            self.right.location.to_simple_location(),
        )

    def __str__(self) -> str:
        u, v, lu, lv = self.to_join_tuple()
        return f'{u}{lu}:{v}{lv}'

    def __repr__(self) -> str:
        return self.__str__()


class AssemblySourceCommonClass:
    assembly: List[AssemblyJoin] = Field(
        default_factory=list, description="""The joins between the fragments in the assembly"""
    )

    def minimal_overlap(self):
        """Returns the minimal overlap between the fragments in the assembly"""
        all_overlaps = list()
        for f in self.assembly:
            # TODO: these conditions are probably not needed
            # if f.left.location is not None:
            all_overlaps.append(f.left.location.end - f.left.location.start)
            # if f.right.location is not None:
            all_overlaps.append(f.right.location.end - f.right.location.start)
        return min(all_overlaps)

    def get_assembly_plan(self):
        """Returns the assembly plan"""
        return tuple(j.to_join_tuple() for j in self.assembly)

    @classmethod
    def from_assembly(
        cls, assembly: list[tuple[int, int, Location, Location]], input: list[int], id: int, circular: bool, **kwargs
    ):
        return cls(
            id=id,
            assembly=[AssemblyJoin.from_join_tuple(join) for join in assembly],
            input=input,
            circular=circular,
            **kwargs,
        )


class AssemblySource(_AssemblySource, AssemblySourceCommonClass):
    pass


class PCRSource(_PCRSource, AssemblySourceCommonClass):
    """Documents a PCR, and the selection of one of the products."""

    # TODO: add this to LinkML
    # input: conlist(int, min_length=1, max_length=1)
    # circular has to be false

    @classmethod
    def from_assembly(
        cls,
        assembly: list[tuple[int, int, Location, Location]],
        input: list[int],
        id: int,
        forward_primer: int,
        reverse_primer: int,
    ):
        """Creates a PCRSource from an assembly, input and id"""
        return super().from_assembly(
            assembly, input, id, False, forward_primer=forward_primer, reverse_primer=reverse_primer
        )


class LigationSource(_LigationSource, AssemblySourceCommonClass):
    pass


class HomologousRecombinationSource(_HomologousRecombinationSource, AssemblySourceCommonClass):

    # TODO: add this to LinkML
    # This can only take two inputs, the first one is the template, the second one is the insert
    # input: conlist(int, min_length=2, max_length=2)
    pass


class GibsonAssemblySource(_GibsonAssemblySource, AssemblySourceCommonClass):

    # TODO: add this to LinkML
    # input: conlist(int, min_length=1)
    pass


class CRISPRSource(_CRISPRSource, AssemblySourceCommonClass):

    # TODO
    # input: conlist(int, min_length=2, max_length=2)
    # circular: bool = False

    @classmethod
    def from_assembly(
        cls,
        assembly: list[tuple[int, int, Location, Location]],
        input: list[int],
        id: int,
        guides: list[int],
    ):
        return super().from_assembly(assembly, input, id, False, guides=guides)


class RestrictionAndLigationSource(_RestrictionAndLigationSource, AssemblySourceCommonClass):
    # TODO: add this to LinkML
    # input: conlist(int, min_length=1)

    @classmethod
    def from_assembly(
        cls,
        assembly: list[tuple[int, int, Location, Location]],
        input: list[int],
        circular: bool,
        id: int,
        restriction_enzymes=list['str'],
    ):
        return super().from_assembly(assembly, input, id, circular, restriction_enzymes=restriction_enzymes)


class OligoHybridizationSource(_OligoHybridizationSource):
    pass


class PolymeraseExtensionSource(_PolymeraseExtensionSource):
    pass
