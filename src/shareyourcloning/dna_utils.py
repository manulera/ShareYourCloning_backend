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


def align_sanger_track(dseqr: Dseqrecord, sanger: str) -> (str, str):
    """Align a sanger track to a dseqr sequence"""
    # Check that required executables exist in PATH
    if not shutil.which('mars'):
        raise RuntimeError("'mars' executable not found in PATH")
    if not shutil.which('mafft'):
        raise RuntimeError("'mafft' executable not found in PATH")
    # Create temporary directory and file
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = os.path.join(tmpdir, 'sequences.fa')
        permuted_fasta_path = os.path.join(tmpdir, 'sequences_permuted.fa')
        aln_path = os.path.join(tmpdir, 'aligned.fa')

        # Write sequences to FASTA file
        with open(fasta_path, 'w') as f:
            f.write(f">seq1\n{dseqr.seq}\n")
            f.write(f">seq2\n{sanger}\n")

        # If the sequence is circular, use MARS to permutate it
        if dseqr.circular:
            result = subprocess.run(['mars', '-a', 'DNA', '-m', '0', '-i', fasta_path, '-o', permuted_fasta_path, '-q', '5', '-l', '20', '-P', '1'], capture_output=True, text=True)  # fmt: skip
            if result.returncode != 0:
                raise RuntimeError(f'MARS failed:\n{result.stderr}')
            aln_input = permuted_fasta_path

        else:
            aln_input = fasta_path

        # Run alignment
        result = subprocess.run(['mafft', '--nuc', '--adjustdirection', aln_input], capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f'MAFFT alignment failed:\n{result.stderr}')
        with open(aln_path, 'w') as f:
            f.write(result.stdout)

        # Read the file and return the sequences
        records = parse(aln_path, 'fasta')

        return [str(records[0].seq), str(records[1].seq)]
