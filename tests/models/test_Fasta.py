import pytest
from Bio.Seq import Seq as BioSeq

from models.fasta import Fasta
from models.fasta_type import FastaType
from models.seq import Seq

# SEQUENCE OPERATIONS

# Fasta.seqs


def test_fasta(fasta):
    fasta = Fasta(fasta, FastaType.GENERIC)
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


def test_fasta_duplicate_id(fasta_duplicate_id):
    fasta = Fasta(fasta_duplicate_id, FastaType.GENERIC)
    with pytest.raises(ValueError):
        _ = list(fasta.seqs)


def test_read_nonexistent_fasta(nonexistent_fasta):
    fasta = Fasta(nonexistent_fasta, FastaType.GENERIC)
    with pytest.raises(FileNotFoundError):
        _ = list(fasta.seqs)


def test_read_empty_fasta(empty_fasta):
    fasta = Fasta(empty_fasta, FastaType.GENERIC)
    with pytest.raises(ValueError):
        _ = list(fasta.seqs)


def test_fasta_no_id(fasta_no_id):
    fasta = Fasta(fasta_no_id, FastaType.GENERIC)
    with pytest.raises(ValueError):
        _ = list(fasta.seqs)


def test_fasta_no_seq(fasta_no_seqs):
    fasta = Fasta(fasta_no_seqs, FastaType.GENERIC)
    with pytest.raises(ValueError):
        _ = list(fasta.seqs)


# Fasta.genome_ids


def test_genome_ids(protein_fasta):
    fasta = Fasta(protein_fasta, FastaType.PROTEIN)
    genome_ids = fasta.genome_ids
    assert genome_ids[0] == "genome1"
    assert genome_ids[1] == "genome2"


def test_genome_ids_in_generic_fasta(fasta):
    fasta = Fasta(fasta, FastaType.GENERIC)
    with pytest.raises(ValueError):
        _ = fasta.genome_ids


# Fasta.n_seqs


def test_n_seqs(fasta):
    fasta = Fasta(fasta, FastaType.GENERIC)
    assert fasta.n_seqs == 3


# Fasta.get_seqs()


def test_get_seqs_right_order(fasta):
    fasta = Fasta(fasta, FastaType.GENERIC)
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


def test_get_seqs_no_ids(fasta):
    fasta = Fasta(fasta, FastaType.GENERIC)
    with pytest.raises(ValueError):
        _ = fasta.get_seqs()


def test_get_seqs_duplicate_ids(fasta):
    fasta = Fasta(fasta, FastaType.GENERIC)
    with pytest.raises(ValueError):
        _ = fasta.get_seqs("seq1", "seq1")


def test_get_seqs_sequence_not_found(fasta):
    fasta = Fasta(fasta, FastaType.GENERIC)
    with pytest.raises(ValueError):
        _ = fasta.get_seqs("seq10")


# Fasta.add_seqs()


def test_add_seq(fasta):
    fasta = Fasta(fasta, FastaType.GENERIC)
    seq = Seq(
        seq=BioSeq("AT"),
        seq_id="seq4",
        seq_description="",
    )
    fasta.add_seqs(seq)
    seqs = list(fasta.seqs)
    assert seqs[3].id == "seq4"
    assert seqs[3].seq == "AT"


def test_add_no_sequences(fasta):
    fasta = Fasta(fasta, FastaType.GENERIC)
    with pytest.raises(ValueError):
        fasta.add_seqs()


# Fasta.remove_seqs()


def test_remove_seqs(fasta):
    fasta = Fasta(fasta, FastaType.GENERIC)
    assert len(list(fasta.seqs)) == 3
    ids = fasta.ids
    assert "seq1" in ids
    fasta.remove_seqs("seq1")
    ids = fasta.ids
    assert "seq1" not in ids
    assert "seq2" in ids
    assert "seq3" in ids
    assert len(list(fasta.seqs)) == 2


def test_remove_no_sequences(fasta):
    fasta = Fasta(fasta, FastaType.GENERIC)
    with pytest.raises(ValueError):
        fasta.remove_seqs()


def test_remove_seqs_not_found(fasta):
    fasta = Fasta(fasta, FastaType.GENERIC)
    with pytest.raises(ValueError):
        fasta.remove_seqs("seq10")


def test_remove_duplicate_seq_id(fasta):
    fasta = Fasta(fasta, FastaType.GENERIC)
    with pytest.raises(ValueError):
        fasta.remove_seqs("seq1", "seq1")


# FILE OPERATIONS

# Fasta.move_fasta()


def test_move_fasta(fasta, tmp_path):
    fasta = Fasta(fasta, FastaType.GENERIC)
    dir_to = tmp_path / "dir_to"
    dir_to.mkdir()
    original_path = fasta.path
    fasta.move_fasta(dir_to)
    assert fasta.path == dir_to / original_path.name
    assert fasta.path.exists()
    assert not original_path.exists()


def test_move_nonexistent_fasta(nonexistent_fasta, tmp_path):
    fasta = Fasta(nonexistent_fasta, FastaType.GENERIC)
    with pytest.raises(FileNotFoundError):
        fasta.move_fasta(tmp_path)


def test_move_empty_fasta(empty_fasta, tmp_path):
    fasta = Fasta(empty_fasta, FastaType.GENERIC)
    with pytest.raises(ValueError):
        fasta.move_fasta(tmp_path)


def test_move_fasta_to_nonexistent_dir(fasta, tmp_path):
    fasta = Fasta(fasta, FastaType.GENERIC)
    dir_to = tmp_path / "dir_to"
    with pytest.raises(FileNotFoundError):
        fasta.move_fasta(dir_to)


def test_move_fasta_to_existing_file_path(fasta, tmp_path):
    fasta = Fasta(fasta, FastaType.GENERIC)
    with pytest.raises(FileExistsError):
        fasta.move_fasta(tmp_path)


# Fasta.rename_fasta()


def test_rename_fasta(fasta):
    fasta = Fasta(fasta, FastaType.GENERIC)
    fasta.rename_fasta("testfile.txt")
    assert fasta.path.name == "testfile.txt"
    assert fasta.path.is_file()


def test_rename_nonexistent_fasta(nonexistent_fasta):
    fasta = Fasta(nonexistent_fasta, FastaType.GENERIC)
    with pytest.raises(FileNotFoundError):
        fasta.rename_fasta("testfile.txt")


def test_rename_empty_filename(fasta):
    fasta = Fasta(fasta, FastaType.GENERIC)
    with pytest.raises(ValueError):
        fasta.rename_fasta("")


def test_rename_filename_no_suffix(fasta):
    fasta = Fasta(fasta, FastaType.GENERIC)
    with pytest.raises(ValueError):
        fasta.rename_fasta("test")


def test_rename_empty_fasta(empty_fasta):
    fasta = Fasta(empty_fasta, FastaType.GENERIC)
    with pytest.raises(ValueError):
        fasta.rename_fasta("test")


def test_rename_fasta_to_existent_file_path(empty_fasta, fasta):
    fasta = Fasta(fasta, FastaType.GENERIC)
    other_file_name = empty_fasta.name
    with pytest.raises(FileExistsError):
        fasta.rename_fasta(other_file_name)


# Fasta.delete_fasta()


def test_delete_fasta(fasta):
    fasta = Fasta(fasta, FastaType.GENERIC)
    fasta.delete_fasta()
    assert fasta.path.is_file() is False


def test_delete_nonexistent_fasta(nonexistent_fasta):
    fasta = Fasta(nonexistent_fasta, FastaType.GENERIC)
    with pytest.raises(FileNotFoundError):
        fasta.delete_fasta()
