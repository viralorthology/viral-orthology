from pathlib import Path

import pytest
from Bio.Seq import Seq as BioSeq

from models.fasta import Fasta
from models.seq import Seq


@pytest.fixture
def generic_fasta(tmp_path):
    path = tmp_path / "fasta.fasta"
    path.write_text(
        ">seq1 seqdescription\nATCG\n>seq2 seqdescription\nATGC\n>seq3 seqdescription\nCATG\n",
        encoding="utf-8",
    )
    return path


def test_fasta(generic_fasta):
    fasta = Fasta(generic_fasta)
    seqs = list(fasta.seqs)
    assert len(seqs) == 3
    assert seqs[0].id == "seq1"
    assert seqs[1].id == "seq2"
    assert seqs[2].id == "seq3"
    assert seqs[0].description == "seq1 seqdescription"
    assert seqs[1].description == "seq2 seqdescription"
    assert seqs[2].description == "seq3 seqdescription"
    assert seqs[0].seq == "ATCG"
    assert seqs[1].seq == "ATGC"
    assert seqs[2].seq == "CATG"


def test_fasta_duplicate_id(tmp_path):
    fasta_path = tmp_path / "fasta.fasta"
    fasta_path.write_text(">seq1\nATGC\n>seq1\nATAT\n")
    with pytest.raises(ValueError):
        _ = Fasta(fasta_path)


def test_read_nonexistent_fasta():
    fasta = Fasta(Path("test.fasta"))
    with pytest.raises(FileNotFoundError):
        _ = list(fasta.seqs)


def test_read_empty_fasta(tmp_path):
    fasta_path = tmp_path / "fasta.fasta"
    fasta_path.touch()
    fasta = Fasta(fasta_path)
    with pytest.raises(ValueError):
        _ = list(fasta.seqs)


def test_fasta_no_id(tmp_path):
    path = tmp_path / "fasta_no_id.fasta"
    path.write_text(">\nATCG\n", encoding="utf-8")
    with pytest.raises(ValueError):
        _ = Fasta(path)


def test_fasta_no_seq(tmp_path):
    path = tmp_path / "fasta_no_id.fasta"
    path.write_text(">seq1\n", encoding="utf-8")
    with pytest.raises(ValueError):
        _ = Fasta(path)


# Fasta.genome_ids


def test_genome_ids(tmp_path):
    path = tmp_path / "protein_fasta.fasta"
    path.write_text(
        ">seq1 genome1 [protein_id=123]\nATGC\n>seq2 genome2 [protein_id=234]\nATCG\n",
        encoding="utf-8",
    )
    fasta = Fasta(path)
    genome_ids = fasta.genome_ids
    assert genome_ids[0] == "genome1"
    assert genome_ids[1] == "genome2"


def test_genome_ids_in_generic_fasta(generic_fasta):
    fasta = Fasta(generic_fasta)
    with pytest.raises(AttributeError):
        _ = fasta.genome_ids


# Fasta.n_seqs


def test_n_seqs(generic_fasta):
    fasta = Fasta(generic_fasta)
    assert fasta.n_seqs == 3


# Fasta.get_seqs()


def test_get_seqs_right_order(generic_fasta):
    fasta = Fasta(generic_fasta)
    seqs = fasta.get_seqs("seq2", "seq1", "seq3")
    assert len(seqs) == 3
    assert seqs[0].id == "seq2"
    assert seqs[1].id == "seq1"
    assert seqs[2].id == "seq3"
    assert seqs[0].description == "seq2 seqdescription"
    assert seqs[1].description == "seq1 seqdescription"
    assert seqs[2].description == "seq3 seqdescription"
    assert seqs[0].seq == "ATGC"
    assert seqs[1].seq == "ATCG"
    assert seqs[2].seq == "CATG"


def test_get_seqs_no_ids(generic_fasta):
    fasta = Fasta(generic_fasta)
    with pytest.raises(ValueError):
        _ = fasta.get_seqs()


def test_get_seqs_duplicate_ids(generic_fasta):
    fasta = Fasta(generic_fasta)
    with pytest.raises(ValueError):
        _ = fasta.get_seqs("seq1", "seq1")


def test_get_seqs_sequence_not_found(generic_fasta):
    fasta = Fasta(generic_fasta)
    with pytest.raises(ValueError):
        _ = fasta.get_seqs("seq10")


# Fasta.add_seqs()


def test_add_seq(generic_fasta):
    fasta = Fasta(generic_fasta)
    seq = Seq(
        seq=BioSeq("AT"),
        seq_id="seq4",
        seq_description="",
    )
    fasta.add_seqs(seq)
    seqs = list(fasta.seqs)
    assert seqs[3].id == "seq4"
    assert seqs[3].seq == "AT"


def test_add_no_sequences(generic_fasta):
    fasta = Fasta(generic_fasta)
    with pytest.raises(ValueError):
        fasta.add_seqs()


# Fasta.remove_seqs()


def test_remove_seqs(generic_fasta):
    fasta = Fasta(generic_fasta)
    assert len(list(fasta.seqs)) == 3
    ids = fasta.ids
    assert "seq1" in ids
    fasta.remove_seqs("seq1")
    ids = fasta.ids
    assert "seq1" not in ids
    assert "seq2" in ids
    assert "seq3" in ids
    assert len(list(fasta.seqs)) == 2


def test_remove_all_seqs(generic_fasta):
    fasta = Fasta(generic_fasta)
    fasta.remove_seqs("seq1", "seq2", "seq3")
    assert not fasta.path.exists()


def test_remove_no_sequences(generic_fasta):
    fasta = Fasta(generic_fasta)
    with pytest.raises(ValueError):
        fasta.remove_seqs()


def test_remove_seqs_not_found(generic_fasta):
    fasta = Fasta(generic_fasta)
    with pytest.raises(ValueError):
        fasta.remove_seqs("seq10")


def test_remove_duplicate_seq_id(generic_fasta):
    fasta = Fasta(generic_fasta)
    with pytest.raises(ValueError):
        fasta.remove_seqs("seq1", "seq1")


# FILE OPERATIONS

# Fasta.move_fasta()


def test_move_fasta(generic_fasta, tmp_path):
    fasta = Fasta(generic_fasta)
    dir_to = tmp_path / "dir_to"
    dir_to.mkdir()
    original_path = fasta.path
    fasta.move_fasta(dir_to)
    assert fasta.path == dir_to / original_path.name
    assert fasta.path.exists()
    assert not original_path.exists()


def test_move_nonexistent_fasta(tmp_path):
    fasta_path = tmp_path / "fasta.fasta"
    fasta = Fasta(fasta_path)
    with pytest.raises(FileNotFoundError):
        fasta.move_fasta(tmp_path)


def test_move_empty_fasta(tmp_path):
    fasta_path = tmp_path / "fasta.fasta"
    fasta_path.touch()
    fasta = Fasta(fasta_path)
    with pytest.raises(ValueError):
        fasta.move_fasta(tmp_path)


def test_move_fasta_to_nonexistent_dir(generic_fasta, tmp_path):
    fasta = Fasta(generic_fasta)
    dir_to = tmp_path / "dir_to"
    with pytest.raises(FileNotFoundError):
        fasta.move_fasta(dir_to)


def test_move_fasta_to_existing_file_path(generic_fasta, tmp_path):
    fasta = Fasta(generic_fasta)
    with pytest.raises(FileExistsError):
        fasta.move_fasta(tmp_path)


# Fasta.rename_fasta()


def test_rename_fasta(generic_fasta):
    fasta = Fasta(generic_fasta)
    fasta.rename_fasta("testfile.txt")
    assert fasta.path.name == "testfile.txt"
    assert fasta.path.is_file()


def test_rename_nonexistent_fasta():
    fasta = Fasta(Path("fasta.fasta"))
    with pytest.raises(FileNotFoundError):
        fasta.rename_fasta("testfile.txt")


def test_rename_empty_filename(generic_fasta):
    fasta = Fasta(generic_fasta)
    with pytest.raises(ValueError):
        fasta.rename_fasta("")


def test_rename_filename_no_suffix(generic_fasta):
    fasta = Fasta(generic_fasta)
    with pytest.raises(ValueError):
        fasta.rename_fasta("test")


def test_rename_empty_fasta(tmp_path):
    fasta_path = tmp_path / "fasta.fasta"
    fasta_path.touch()
    fasta = Fasta(fasta_path)
    with pytest.raises(ValueError):
        fasta.rename_fasta("test")


def test_rename_fasta_to_existent_file_path(generic_fasta, tmp_path):
    other_file_path = tmp_path / "file.file"
    other_file_path.touch()
    fasta = Fasta(generic_fasta)
    other_file_name = other_file_path.name
    with pytest.raises(FileExistsError):
        fasta.rename_fasta(other_file_name)


# Fasta.delete_fasta()


def test_delete_fasta(generic_fasta):
    fasta = Fasta(generic_fasta)
    fasta.delete_fasta()
    assert fasta.path.exists() is False


def test_delete_nonexistent_fasta():
    fasta = Fasta(Path("fasta.fasta"))
    with pytest.raises(FileNotFoundError):
        fasta.delete_fasta()
