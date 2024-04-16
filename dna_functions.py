from functools import cmp_to_key
from urllib.error import HTTPError
from Bio.Restriction.Restriction import RestrictionBatch
from Bio.Seq import reverse_complement
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from pydantic_models import RepositoryIdSource, GenbankSequence, SequenceEntity
from pydna.parsers import parse as pydna_parse
import requests
from bs4 import BeautifulSoup
import regex
from Bio.SeqFeature import SimpleLocation, Location, CompoundLocation
from pydna.utils import shift_location
from pydna.common_sub_strings import common_sub_strings


def sum_is_sticky(three_prime_end: tuple[str, str], five_prime_end: tuple[str, str], partial: bool = False) -> int:
    """Return true if the 3' end of seq1 and 5' end of seq2 ends are sticky and compatible for ligation."""
    type_seq1, sticky_seq1 = three_prime_end
    type_seq2, sticky_seq2 = five_prime_end

    if 'blunt' != type_seq2 and type_seq2 == type_seq1 and str(sticky_seq2) == str(reverse_complement(sticky_seq1)):
        return len(sticky_seq1)

    if not partial:
        return 0

    if type_seq1 != type_seq2 or type_seq2 == 'blunt':
        return 0
    elif type_seq2 == "5'":
        sticky_seq1 = str(reverse_complement(sticky_seq1))
    elif type_seq2 == "3'":
        sticky_seq2 = str(reverse_complement(sticky_seq2))

    ovhg_len = min(len(sticky_seq1), len(sticky_seq2))
    # [::-1] to try the longest overhangs first
    for i in range(1, ovhg_len + 1)[::-1]:
        if sticky_seq1[-i:] == sticky_seq2[:i]:
            return i
    else:
        return 0


def format_sequence_genbank(seq: Dseqrecord) -> SequenceEntity:

    if seq.name.lower() == 'exported':
        correct_name(seq)
    # In principle here we do not set the id from the id of the Dseqrecord
    gb_seq = GenbankSequence(
        file_content=seq.format('genbank'),
        overhang_crick_3prime=seq.seq.ovhg,
        overhang_watson_3prime=seq.seq.watson_ovhg(),
    )
    return SequenceEntity(sequence=gb_seq)


def read_dsrecord_from_json(seq: SequenceEntity) -> Dseqrecord:
    initial_dseqrecord: Dseqrecord = pydna_parse(seq.sequence.file_content)[0]
    if seq.sequence.overhang_watson_3prime == 0 and seq.sequence.overhang_crick_3prime == 0:
        out_dseq_record = initial_dseqrecord
    else:
        out_dseq_record = Dseqrecord(
            Dseq.from_full_sequence_and_overhangs(
                str(initial_dseqrecord.seq), seq.sequence.overhang_crick_3prime, seq.sequence.overhang_watson_3prime
            ),
            features=initial_dseqrecord.features,
        )
    # We set the id to the integer converted to integer (this is only
    # useful for assemblies)
    out_dseq_record.id = str(seq.id)
    return out_dseq_record


def get_invalid_enzyme_names(enzyme_names_list: list[str | None]) -> list[str]:
    rest_batch = RestrictionBatch()
    invalid_names = list()
    for name in enzyme_names_list:
        # Empty enzyme names are the natural edges of the molecule
        if name is not None:
            try:
                rest_batch.format(name)
            except ValueError:
                invalid_names.append(name)
    return invalid_names


def get_sequences_from_gb_file_url(url: str) -> Dseqrecord:
    # TODO once pydna parse is fixed it should handle urls that point to non-gb files
    resp = requests.get(url)
    if resp.status_code != 200:
        raise HTTPError(url, 404, 'file requested from url not found', 'file requested from url not found', None)
    return pydna_parse(resp.text)


def request_from_addgene(source: RepositoryIdSource) -> tuple[list[Dseqrecord], list[RepositoryIdSource]]:
    # TODO here maybe it would be good to check that the addgeneID still returns the url requested.
    if 'url' in source.info:
        return [get_sequences_from_gb_file_url(source.info['url'])[0]], [source]

    url = f'https://www.addgene.org/{source.repository_id}/sequences/'
    resp = requests.get(url)
    if resp.status_code == 404:
        raise HTTPError(url, 404, 'wrong addgene id', 'wrong addgene id', None)
    soup = BeautifulSoup(resp.content, 'html.parser')
    sequence_file_url_dict = dict()
    for _type in ['depositor-full', 'depositor-partial', 'addgene-full', 'addgene-partial']:
        sequence_file_url_dict[_type] = []
        if soup.find(id=_type) is not None:
            sequence_file_url_dict[_type] = [
                a.get('href') for a in soup.find(id=_type).findAll(class_='genbank-file-download')
            ]

    # TODO provide addgene sequencing data supporting the sequence
    # We prefer to return addgene full if both available
    products = list()
    sources = list()
    for _type in ['addgene-full', 'depositor-full']:
        if len(sequence_file_url_dict[_type]) > 0:
            for seq_url in sequence_file_url_dict[_type]:
                new_source = source.model_copy()
                new_source.info = {'url': seq_url, 'type': _type}
                sources.append(new_source)
                # There should be only one sequence
                products.append(get_sequences_from_gb_file_url(seq_url)[0])
            break
    return products, sources


def correct_name(dseq: Dseqrecord):
    # Can set the name from keyword if locus is set to Exported
    if dseq.name.lower() == 'exported' and dseq.locus.lower() == 'exported' and 'keywords' in dseq.annotations:
        dseq.name = dseq.annotations['keywords'][0]


def location_sorter(x, y) -> int:
    """
    Sort by start, then length, then strand.
    """
    if x.parts[0].start != y.parts[0].start:
        return x.parts[0].start - y.parts[0].start
    elif x.parts[-1].end != y.parts[-1].end:
        return x.parts[-1].end - y.parts[-1].end
    return x.strand - y.strand


def get_all_regex_feature_edges(pattern: str, seq: str, is_circular: bool) -> list[tuple[int, int]]:

    subject = 2 * seq if is_circular else seq

    compiled_pattern = regex.compile(pattern, regex.IGNORECASE)
    compiled_pattern_rev = regex.compile('(?r)' + pattern, regex.IGNORECASE)

    matches = list(regex.finditer(compiled_pattern, subject, overlapped=True))
    matches += list(regex.finditer(compiled_pattern_rev, subject, overlapped=True))

    # In circular objects we remove the matches that span the sequence more than once: m.end() - m.start() <= len(seq)
    return list(set([(m.start(), m.end()) for m in matches if (m.end() - m.start() <= len(seq))]))


def find_sequence_regex(pattern: str, seq: str, is_circular: bool) -> list[Location]:

    feature_locations = list()

    # Strand 1
    feature_edges = get_all_regex_feature_edges(pattern, seq, is_circular)
    # We use shift_location to format origin-spanning features in circular DNA
    feature_locations += [shift_location(SimpleLocation(start, end, 1), 0, len(seq)) for start, end in feature_edges]

    # Strand -1
    feature_edges = get_all_regex_feature_edges(pattern, reverse_complement(seq), is_circular)
    feature_locations += [
        shift_location(SimpleLocation(start, end, 1)._flip(len(seq)), 0, len(seq)) for start, end in feature_edges
    ]

    # We return a unique list, cannot use a set because Location is not hashable
    return sorted(
        [x for i, x in enumerate(feature_locations) if x not in feature_locations[:i]], key=cmp_to_key(location_sorter)
    )


def location_edges(location: Location):
    if isinstance(location, SimpleLocation):
        return location.start, location.end
    if isinstance(location, CompoundLocation):
        return location.parts[0].start, location.parts[-1].end


def get_homologous_recombination_locations(
    template: Dseqrecord, insert: Dseqrecord, minimal_homology
) -> list[Location]:
    """Return the locations of the possible homologous recombination sites."""
    template_seq = str(template.seq)
    insert_seq = str(insert.seq)
    regex_pattern = insert_seq[:minimal_homology] + '.*' + insert_seq[-minimal_homology:]
    locations = find_sequence_regex(regex_pattern, template_seq, template.circular)
    return locations


def perform_homologous_recombination(template: Dseqrecord, insert: Dseqrecord, location: Location):
    edges = location_edges(location)
    if template.circular:
        return (template[edges[1] : edges[0]] + insert).looped()
    return template[0 : edges[0]] + insert + template[edges[1] :]


def seq_overlap_length(dseq: Dseq) -> int:
    return len(dseq) - abs(dseq.ovhg) - abs(dseq.watson_ovhg())


def oligo_hybridization_overhangs(fwd_oligo_seq: str, rvs_oligo_seq: str, minimal_annealing: int) -> list[int]:
    matches = common_sub_strings(fwd_oligo_seq.lower(), reverse_complement(rvs_oligo_seq.lower()), minimal_annealing)
    # Return possible overhangs
    return [start_on_rvs - start_on_fwd for start_on_fwd, start_on_rvs, length in matches]
