import pytest
from Bio.Seq import Seq as BioSeq

from models.seq import Seq

# FASTA OPERATIONS


def test_fasta_ids(fasta):
    ids = fasta.ids
    assert len(ids) == 3
    assert ids[0] == "seq1"


def test_duplicate_id(fasta_duplicate_id):
    with pytest.raises(ValueError):
        _ = list(fasta_duplicate_id.seqs)


def test_read_nonexistent_fasta(nonexistent_fasta):
    with pytest.raises(FileNotFoundError):
        _ = list(nonexistent_fasta.seqs)


def test_read_empty_fasta(empty_fasta):
    with pytest.raises(ValueError):
        _ = list(empty_fasta.seqs)


def test_fasta_no_id(fasta_no_id):
    with pytest.raises(ValueError):
        _ = list(fasta_no_id.seqs)


def test_get_seqs(fasta):
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
    with pytest.raises(AssertionError):
        _ = fasta.get_seqs()


def test_get_seqs_duplicate_ids(fasta):
    with pytest.raises(AssertionError):
        _ = fasta.get_seqs("seq1", "seq1")


def test_get_seqs_sequence_not_found(fasta):
    with pytest.raises(ValueError):
        _ = fasta.get_seqs("seq10")


def test_add_no_sequences(fasta):
    with pytest.raises(AssertionError):
        fasta.add_seqs()


def test_remove_no_sequences(fasta):
    with pytest.raises(AssertionError):
        fasta.remove_seqs()


def test_remove_seqs(fasta):
    assert len(list(fasta.seqs)) == 3
    ids = fasta.ids
    assert "seq1" in ids
    fasta.remove_seqs("seq1")
    ids = fasta.ids
    assert "seq1" not in ids
    assert "seq2" in ids
    assert "seq3" in ids
    assert len(list(fasta.seqs)) == 2


def test_remove_seqs_not_found(fasta):
    with pytest.raises(ValueError):
        fasta.remove_seqs("seq10")


def test_remove_duplicate_seq_id(fasta):
    with pytest.raises(AssertionError):
        fasta.remove_seqs("seq1", "seq1")


def test_genome_ids_in_generic_fasta(fasta):
    with pytest.raises(AssertionError):
        _ = fasta.genome_ids()


def test_genome_ids(protein_fasta):
    genome_ids = protein_fasta.genome_ids
    assert genome_ids[0] == "genome1"
    assert genome_ids[1] == "genome2"


def test_add_seq(fasta):
    seq = Seq(
        seq=BioSeq("AT"),
        seq_id="seq4",
        seq_description="",
    )
    fasta.add_seqs(seq)
    seqs = list(fasta.seqs)
    assert seqs[3].id == "seq4"
    assert seqs[3].seq == "AT"


def test_n_seqs(fasta):
    assert fasta.n_seqs == 3


# FILE OPERATIONS


def test_move_nonexistent_fasta(nonexistent_fasta, tmp_path):
    with pytest.raises(FileNotFoundError):
        nonexistent_fasta.move_fasta(tmp_path)


def test_move_empty_fasta(empty_fasta, tmp_path):
    with pytest.raises(ValueError):
        empty_fasta.move_fasta(tmp_path)


def test_move_fasta(fasta, tmp_path):
    dir_to = tmp_path / "dir_to"
    dir_to.mkdir()
    original_path = fasta.path
    fasta.move_fasta(dir_to)
    assert fasta.path == dir_to / original_path.name
    assert fasta.path.exists()
    assert not original_path.exists()


def test_move_fasta_to_nonexistent_dir(fasta, tmp_path):
    dir_to = tmp_path / "dir_to"
    with pytest.raises(FileNotFoundError):
        fasta.move_fasta(dir_to)


def test_move_fasta_to_existing_file_path(fasta, tmp_path):
    with pytest.raises(FileExistsError):
        fasta.move_fasta(tmp_path)


def test_rename_nonexistent_fasta(nonexistent_fasta):
    with pytest.raises(FileNotFoundError):
        nonexistent_fasta.rename_fasta("testfile.txt")


def test_rename_fasta(fasta):
    fasta.rename_fasta("testfile.txt")
    assert fasta.path.name == "testfile.txt"
    assert fasta.path.is_file()


def test_rename_empty_filename(fasta):
    with pytest.raises(AssertionError):
        fasta.rename_fasta("")


def test_rename_filename_no_suffix(fasta):
    with pytest.raises(ValueError):
        fasta.rename_fasta("test")


def test_rename_empty_fasta(empty_fasta, fasta):
    other_file_name = fasta.path.name
    with pytest.raises(ValueError):
        empty_fasta.rename_fasta(other_file_name)


def test_rename_fasta_to_existent_file_path(empty_fasta, fasta):
    other_file_name = empty_fasta.path.name
    with pytest.raises(FileExistsError):
        fasta.rename_fasta(other_file_name)


def test_delete_nonexistent_fasta(nonexistent_fasta):
    with pytest.raises(FileNotFoundError):
        nonexistent_fasta.delete_fasta()


def test_delete_fasta(fasta):
    fasta.delete_fasta()
    assert fasta.path.is_file() is False
