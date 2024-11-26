from pydna.dseqrecord import Dseqrecord
from pydantic_models import SequenceFileFormat
from Bio.SeqIO import parse as seqio_parse
import io
import warnings
from Bio.SeqIO.InsdcIO import GenBankIterator, GenBankScanner
import re


class MyGenBankScanner(GenBankScanner):
    def _feed_first_line(self, consumer, line):
        # A regex for LOCUS       pKM265       4536 bp    DNA   circular  SYN        21-JUN-2013
        m = re.match(
            r'(?i)LOCUS\s+(?P<name>\S+)\s+(?P<size>\d+ bp)\s+(?P<molecule_type>\S+)\s+(?P<topology>(?:circular|linear))\s+.+(?P<date>\d{2}-\w{3}-\d{4})',
            line,
        )
        if m is None:
            raise ValueError('LOCUS line cannot be parsed')
        name, size, molecule_type, topology, date = m.groups()
        # If topology not present, error
        if topology is None:
            raise ValueError('LOCUS line does not contain topology')
        consumer.locus(name)
        consumer.size(size[:-3])
        consumer.molecule_type(molecule_type)
        consumer.topology(topology.lower())
        consumer.date(date)


class MyGenBankIterator(GenBankIterator):

    def parse(self, handle):
        """Start parsing the file, and return a SeqRecord generator."""
        records = MyGenBankScanner(debug=0).parse_records(handle)
        return records


def custom_file_parser(
    file_streamer: io.BytesIO | io.StringIO, sequence_file_format: SequenceFileFormat, circularize: bool = False
) -> list[Dseqrecord]:
    """
    Parse a file with SeqIO.parse (specifying the format and using the topology annotation to determine circularity).

    If the format is genbank and the parsing of the LOCUS line fails, fallback to custom regex-based parsing.
    """

    out = list()

    with file_streamer as handle:
        try:
            for parsed_seq in seqio_parse(handle, sequence_file_format):
                circularize = circularize or (
                    'topology' in parsed_seq.annotations.keys() and parsed_seq.annotations['topology'] == 'circular'
                )

                if 'topology' not in parsed_seq.annotations.keys():
                    raise ValueError('LOCUS line does not contain topology')

                out.append(Dseqrecord(parsed_seq, circular=circularize))
        except ValueError as e:
            # If the error is about the LOCUS line, we try to parse with regex
            warnings.warn(
                'LOCUS line is wrongly formatted, we used a more permissive parser.',
                stacklevel=2,
            )
            # Reset the file handle position to the start since we consumed it in the first attempt
            if 'LOCUS line does not contain' in str(e):
                handle.seek(0)
                out = list()
                for parsed_seq in MyGenBankIterator(handle):
                    circularize = circularize or (
                        'topology' in parsed_seq.annotations.keys()
                        and parsed_seq.annotations['topology'] == 'circular'
                    )
                    print(f'Circularize: {circularize}')
                    out.append(Dseqrecord(parsed_seq, circular=circularize))
            else:
                raise e

    return out


if __name__ == '__main__':
    with open('/Users/dgruano/dev/learn/Genestorian/cnrs-training/Re Gateway feature/P2RP3.ape', 'r') as f:
        circularize = False
        # seq = custom_file_parser(f, 'gb')[0]
        out = list()
        for parsed_seq in MyGenBankIterator(f):
            # print(parsed_seq)
            # print(parsed_seq.annotations)
            # print(parsed_seq.annotations['topology'] == 'circular')
            pass

        f.seek(0)
        plasmid = custom_file_parser(f, 'gb')[0]
        print(plasmid.circular)
    pass
