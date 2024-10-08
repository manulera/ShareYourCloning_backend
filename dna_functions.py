from functools import cmp_to_key
from urllib.error import HTTPError
from Bio.Restriction.Restriction import RestrictionBatch
from Bio.Seq import reverse_complement
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from pydantic_models import TextFileSequence, AddGeneIdSource, SequenceFileFormat
from pydna.parsers import parse as pydna_parse
import requests
from bs4 import BeautifulSoup
import regex
from Bio.SeqFeature import SimpleLocation, Location
from pydna.utils import shift_location
from pydna.common_sub_strings import common_sub_strings


def format_sequence_genbank(seq: Dseqrecord, seq_name: str = None) -> TextFileSequence:

    if seq_name is not None:
        seq.name = seq_name
    elif seq.name.lower() == 'exported':
        correct_name(seq)

    return TextFileSequence(
        id=0,
        file_content=seq.format('genbank'),
        sequence_file_format=SequenceFileFormat('genbank'),
        overhang_crick_3prime=seq.seq.ovhg,
        overhang_watson_3prime=seq.seq.watson_ovhg(),
    )


def read_dsrecord_from_json(seq: TextFileSequence) -> Dseqrecord:
    initial_dseqrecord: Dseqrecord = pydna_parse(seq.file_content)[0]
    if seq.overhang_watson_3prime == 0 and seq.overhang_crick_3prime == 0:
        out_dseq_record = initial_dseqrecord
    else:
        out_dseq_record = Dseqrecord(
            Dseq.from_full_sequence_and_overhangs(
                str(initial_dseqrecord.seq), seq.overhang_crick_3prime, seq.overhang_watson_3prime
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


def get_sequences_from_gb_file_url(url: str) -> list[Dseqrecord]:
    # TODO once pydna parse is fixed it should handle urls that point to non-gb files
    resp = requests.get(url)
    if resp.status_code != 200:
        raise HTTPError(url, 404, 'file requested from url not found', 'file requested from url not found', None)
    return pydna_parse(resp.text)


def request_from_addgene(source: AddGeneIdSource) -> tuple[Dseqrecord, AddGeneIdSource]:

    url = f'https://www.addgene.org/{source.repository_id}/sequences/'
    resp = requests.get(url)
    if resp.status_code == 404:
        raise HTTPError(url, 404, 'wrong addgene id', 'wrong addgene id', None)
    soup = BeautifulSoup(resp.content, 'html.parser')

    # Get a span.material-name from the soup, see https://github.com/manulera/ShareYourCloning_backend/issues/182
    plasmid_name = soup.find('span', class_='material-name').text

    if source.sequence_file_url:
        dseqr = get_sequences_from_gb_file_url(source.sequence_file_url)[0]
        dseqr.name = plasmid_name
        return dseqr, source

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
                new_source.sequence_file_url = seq_url
                new_source.addgene_sequence_type = _type
                sources.append(new_source)
                # There should be only one sequence
                products.append(get_sequences_from_gb_file_url(seq_url)[0])

    if len(products) == 0:
        # They may have only partial sequences
        raise HTTPError(
            url,
            404,
            f'The requested plasmid does not have full sequences, see https://www.addgene.org/{source.repository_id}/sequences/',
            f'The requested plasmid does not have full sequences, see https://www.addgene.org/{source.repository_id}/sequences/',
            None,
        )

    # Rename the plasmid
    for p in products:
        p.name = plasmid_name
    return products[0], sources[0]


def correct_name(dseq: Dseqrecord):
    # Can set the name from keyword if locus is set to Exported
    if dseq.name.lower() == 'exported' and dseq.locus.lower() == 'exported' and 'keywords' in dseq.annotations:
        dseq.name = dseq.annotations['keywords'][0].replace(' ', '_')


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


def get_homologous_recombination_locations(
    template: Dseqrecord, insert: Dseqrecord, minimal_homology
) -> list[Location]:
    """Return the locations of the possible homologous recombination sites."""
    template_seq = str(template.seq)
    insert_seq = str(insert.seq)
    regex_pattern = insert_seq[:minimal_homology] + '.*' + insert_seq[-minimal_homology:]
    locations = find_sequence_regex(regex_pattern, template_seq, template.circular)
    return locations


# Could be useful at some point
# def seq_overlap_length(dseq: Dseq) -> int:
#     return len(dseq) - abs(dseq.ovhg) - abs(dseq.watson_ovhg())


def oligonucleotide_hybridization_overhangs(
    fwd_oligo_seq: str, rvs_oligo_seq: str, minimal_annealing: int
) -> list[int]:
    matches = common_sub_strings(fwd_oligo_seq.lower(), reverse_complement(rvs_oligo_seq.lower()), minimal_annealing)

    for m in matches:
        if not (
            (m[0] == 0 and m[1] + m[2] == len(fwd_oligo_seq)) or (m[1] == 0 and m[0] + m[2] == len(rvs_oligo_seq))
        ):
            raise ValueError('The oligonucleotides can anneal with mismatches')

    # Return possible overhangs
    return [start_on_rvs - start_on_fwd for start_on_fwd, start_on_rvs, length in matches]
