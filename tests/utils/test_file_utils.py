import pytest
from Bio.Seq import Seq as BioSeq

from models.fasta import Fasta
from models.fasta_type import FastaType
from models.seq import Seq
from utils import (
    ensure_files_exist,
    ensure_files_have_content,
    ensure_paths_do_not_exist,
    get_combined_fasta,
    get_fastas,
    get_seqs_from_fasta_str,
)


def test_get_fastas_returns_files_with_matching_extension(tmp_path):
    fasta_1 = tmp_path / "a.fasta"
    fasta_2 = tmp_path / "b.fasta"
    other_file = tmp_path / "c.txt"

    fasta_1.touch()
    fasta_2.touch()
    other_file.touch()

    result = get_fastas(tmp_path, FastaType.GENERIC, ".fasta")

    assert {f.path for f in result} == {fasta_1, fasta_2}


def test_get_fastas_invalid_extension(tmp_path):
    with pytest.raises(ValueError):
        get_fastas(tmp_path, FastaType.GENERIC, "fasta")


def test_ensure_files_exist_existing_files(tmp_path):
    file_1 = tmp_path / "file_1.txt"
    file_2 = tmp_path / "file_2.txt"

    file_1.touch()
    file_2.touch()

    ensure_files_exist(file_1, file_2)


def test_ensure_files_exist_missing_file(tmp_path):
    missing_file = tmp_path / "missing.txt"

    with pytest.raises(FileNotFoundError):
        ensure_files_exist(missing_file)


def test_ensure_files_do_not_exist_missing_files(tmp_path):
    ensure_paths_do_not_exist(
        tmp_path / "missing_1.txt",
        tmp_path / "missing_2.txt",
    )


def test_ensure_files_do_not_exist_existing_file(tmp_path):
    existing_file = tmp_path / "existing.txt"
    existing_file.touch()

    with pytest.raises(FileExistsError):
        ensure_paths_do_not_exist(existing_file)


def test_ensure_files_have_content_non_empty_files(tmp_path):
    file_1 = tmp_path / "file_1.txt"
    file_2 = tmp_path / "file_2.txt"

    file_1.write_text("content")
    file_2.write_text("more content")

    ensure_files_have_content(file_1, file_2)


def test_ensure_files_have_content_empty_file(tmp_path):
    empty_file = tmp_path / "empty.txt"
    empty_file.touch()

    with pytest.raises(ValueError):
        ensure_files_have_content(empty_file)


def test_ensure_files_have_content_missing_file(tmp_path):
    missing_file = tmp_path / "missing.txt"

    with pytest.raises(FileNotFoundError):
        ensure_files_have_content(missing_file)


def test_get_seqs_from_fasta_str():
    fasta_str = """>seq1 description one
ATGC
>seq2 description two
GGTA
"""

    result = get_seqs_from_fasta_str(fasta_str)

    assert len(result) == 2

    assert result[0].id == "seq1"
    assert result[0].description == "seq1 description one"
    assert str(result[0].seq) == "ATGC"

    assert result[1].id == "seq2"
    assert result[1].description == "seq2 description two"
    assert str(result[1].seq) == "GGTA"


def test_get_seqs_from_fasta_str_empty_string():
    with pytest.raises(ValueError):
        get_seqs_from_fasta_str("")


def test_get_combined_fasta(tmp_path):
    fasta_1 = Fasta(tmp_path / "a.fasta", FastaType.GENERIC)
    fasta_1.add_seqs(
        Seq(
            seq=BioSeq("ATGC"),
            seq_id="seq1",
            seq_description="first",
        )
    )

    fasta_2 = Fasta(tmp_path / "b.fasta", FastaType.GENERIC)
    fasta_2.add_seqs(
        Seq(
            seq=BioSeq("GGTA"),
            seq_id="seq2",
            seq_description="second",
        )
    )

    output_path = tmp_path / "combined.fasta"

    result = get_combined_fasta(output_path, fasta_1, fasta_2)

    assert result.path == output_path
    assert result.fasta_type == FastaType.GENERIC
    assert result.n_seqs == 2

    result_seqs = list(result.seqs)

    assert result_seqs[0].id == "seq1"
    assert str(result_seqs[0].seq) == "ATGC"
    assert result_seqs[0].description == "seq1 first"

    assert result_seqs[1].id == "seq2"
    assert str(result_seqs[1].seq) == "GGTA"
    assert result_seqs[1].description == "seq2 second"


def test_get_combined_fasta_requires_at_least_one_fasta(tmp_path):
    with pytest.raises(ValueError):
        get_combined_fasta(tmp_path / "combined.fasta")


def test_get_combined_fasta_existing_output(tmp_path):
    output_path = tmp_path / "combined.fasta"
    output_path.touch()

    fasta = Fasta(tmp_path / "input.fasta", FastaType.GENERIC)

    with pytest.raises(ValueError):
        get_combined_fasta(output_path, fasta)
