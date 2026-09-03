from pathlib import Path
from unittest.mock import MagicMock, patch

from engines.muscle import muscle_align


def test_muscle_align_uses_super5_for_more_than_500_sequences():
    fasta = MagicMock()
    fasta.n_seqs = 501
    fasta.path = Path("/tmp/proteins.fasta")
    fasta.alignment_path = Path("/tmp/proteins.alignment.fasta")

    with patch("engines.muscle.utils.run_cmd") as run_cmd:
        muscle_align(fasta)

    run_cmd.assert_called_once_with(
        f"muscle -super5 {fasta.path} -output {fasta.alignment_path}"
    )


def test_muscle_align_uses_align_for_500_sequences():
    fasta = MagicMock()
    fasta.n_seqs = 500
    fasta.path = Path("/tmp/proteins.fasta")
    fasta.alignment_path = Path("/tmp/proteins.alignment.fasta")

    with patch("engines.muscle.utils.run_cmd") as run_cmd:
        muscle_align(fasta)

    run_cmd.assert_called_once_with(
        f"muscle -align {fasta.path} -output {fasta.alignment_path}"
    )


def test_muscle_align_uses_align_for_less_than_500_sequences():
    fasta = MagicMock()
    fasta.n_seqs = 100
    fasta.path = Path("/tmp/proteins.fasta")
    fasta.alignment_path = Path("/tmp/proteins.alignment.fasta")

    with patch("engines.muscle.utils.run_cmd") as run_cmd:
        muscle_align(fasta)

    run_cmd.assert_called_once_with(
        f"muscle -align {fasta.path} -output {fasta.alignment_path}"
    )
