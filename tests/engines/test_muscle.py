import multiprocessing
from unittest.mock import MagicMock, patch

from engines.muscle import muscle_align


def test_muscle_align_uses_super5_for_more_than_500_sequences(tmp_path):
    fasta = MagicMock()
    fasta.n_seqs = 501
    fasta.path = tmp_path / "proteins.fasta"
    fasta.alignment_path = tmp_path / "proteins.alignment.fasta"
    num_threads = max(1, multiprocessing.cpu_count() - 1)

    with patch("engines.muscle.utils.run_cmd") as run_cmd:
        muscle_align(fasta)

    run_cmd.assert_called_once_with(
        f"muscle -super5 {fasta.path} -output {fasta.alignment_path} -threads {num_threads}"
    )


def test_muscle_align_uses_align_for_500_sequences(tmp_path):
    fasta = MagicMock()
    fasta.n_seqs = 500
    fasta.path = tmp_path / "proteins.fasta"
    fasta.alignment_path = tmp_path / "proteins.alignment.fasta"
    num_threads = max(1, multiprocessing.cpu_count() - 1)

    with patch("engines.muscle.utils.run_cmd") as run_cmd:
        muscle_align(fasta)

    run_cmd.assert_called_once_with(
        f"muscle -align {fasta.path} -output {fasta.alignment_path} -threads {num_threads}"
    )


def test_muscle_align_uses_align_for_less_than_500_sequences(tmp_path):
    fasta = MagicMock()
    fasta.n_seqs = 100
    fasta.path = tmp_path / "proteins.fasta"
    fasta.alignment_path = tmp_path / "proteins.alignment.fasta"
    num_threads = max(1, multiprocessing.cpu_count() - 1)

    with patch("engines.muscle.utils.run_cmd") as run_cmd:
        muscle_align(fasta)

    run_cmd.assert_called_once_with(
        f"muscle -align {fasta.path} -output {fasta.alignment_path} -threads {num_threads}"
    )
