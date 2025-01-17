"""
Utility functions moved here to avoid circular imports.
"""

from Bio.Seq import reverse_complement
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
import tempfile
import subprocess
import os
import shutil
from pydna.parsers import parse
from Bio.Align import PairwiseAligner

aligner = PairwiseAligner()
aligner.open_gap_score = -1.0
aligner.extend_gap_score = 0.0
aligner.mode = 'local'


def sum_is_sticky(three_prime_end: tuple[str, str], five_prime_end: tuple[str, str], partial: bool = False) -> int:
    """Return the overlap length if the 3' end of seq1 and 5' end of seq2 ends are sticky and compatible for ligation.
    Return 0 if they are not compatible."""
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


def get_alignment_shift(alignment: Dseq, shift: int) -> int:
    """Shift the alignment by the given number of positions, ignoring gap characters (-).

    Parameters
    ----------
    alignment : Dseq
        The alignment sequence that may contain gap characters (-)
    shift : int
        Number of positions to shift the sequence by

    """

    nucleotides_shifted = 0
    positions_shifted = 0
    corrected_shift = shift if shift >= 0 else len(alignment) + shift
    alignment_str = str(alignment)

    while nucleotides_shifted != corrected_shift:
        if alignment_str[positions_shifted] != '-':
            nucleotides_shifted += 1
        positions_shifted += 1

    return positions_shifted


def permutate_traces(dseqr: Dseqrecord, sanger_traces: list[str], tmpdir: str) -> list[str]:
    """Permutate the traces to account for circular DNA"""
    # As an input for MARS, we need the reference + all traces
    # We include traces in both directions, since MARS does not handle
    # reverse complements - see https://github.com/lorrainea/MARS/issues/17#issuecomment-2598314356
    all_path = os.path.join(tmpdir, 'all.fa')
    with open(all_path, 'w') as f:
        f.write(f">ref\n{dseqr.seq}\n")
        for i, sanger_trace in enumerate(sanger_traces):
            f.write(f">trace-{i+1}\n{sanger_trace}\n")
            f.write(f">trace-rc-{i+1}\n{reverse_complement(sanger_trace)}\n")

    all_permutated_path = os.path.join(tmpdir, 'all_permutated.fa')
    result = subprocess.run(['mars', '-a', 'DNA', '-m', '0', '-i', all_path, '-o', all_permutated_path, '-q', '5', '-l', '20', '-P', '1'], capture_output=True, text=True)  # fmt: skip
    if result.returncode != 0:
        raise RuntimeError(f'MARS failed:\n{result.stderr}')

    # read permutated traces
    return [str(s.seq) for s in parse(all_permutated_path, 'fasta')[1:]]


def align_sanger_traces(dseqr: Dseqrecord, sanger_traces: list[str]) -> list[str]:
    """Align a sanger track to a dseqr sequence"""
    # Check that required executables exist in PATH
    if not shutil.which('mars'):
        raise RuntimeError("'mars' executable not found in PATH")
    if not shutil.which('mafft'):
        raise RuntimeError("'mafft' executable not found in PATH")
    # Create temporary directory and file
    with tempfile.TemporaryDirectory() as tmpdir:

        # Write reference sequence to FASTA file
        reference_path = os.path.join(tmpdir, 'reference.fa')
        with open(reference_path, 'w') as f:
            f.write(f">ref\n{dseqr.seq}\n")

        # If the sequence is circular, use MARS to permutate the traces
        if dseqr.circular:
            traces = permutate_traces(dseqr, sanger_traces, tmpdir)
        else:
            traces = sum([[t, reverse_complement(t)] for t in sanger_traces], [])

        # alignments = []
        # Pairwise-align and keep the best alignment
        for fwd, rvs in zip(traces[::2], traces[1::2]):
            pass
            # TODO Do with mafft
            # if fwd_align.score > rvs_align.score:
            #     alignments.append(fwd_align)
            # else:
            #     alignments.append(rvs_align)

        # for a in alignments[0]:
        #     print(a)

        # Read the file and return the sequences
        return []
