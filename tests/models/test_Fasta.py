from pathlib import Path

import pytest
from Bio.Seq import Seq as BioSeq

from models.fasta import Fasta
from models.seq import Seq


@pytest.fixture
def generic_fasta(tmp_path):
    path = tmp_path / "fasta.fasta"
    path.write_text(
        ">generic_seq1 seqdescription\nATCG\n>generic_seq2 seqdescription\nATGC\n>generic_seq3 seqdescription\nCATG\n",
        encoding="utf-8",
    )
    return Fasta(path)


@pytest.fixture
def protein_fasta(tmp_path):
    path = tmp_path / "fasta.fasta"
    path.write_text(
        ">protein_seq1 genome_1 [protein_id=x]\nATCG\n>protein_seq2 genome_2 [protein_id=x]\nATGC\n>protein_seq3 genome_3 [protein_id=x]\nCATG\n",
        encoding="utf-8",
    )
    return Fasta(path)


def test_fasta(generic_fasta):
    seqs = list(generic_fasta.seqs)
    assert len(seqs) == 3
    assert seqs[0].id == "generic_seq1"
    assert seqs[1].id == "generic_seq2"
    assert seqs[2].id == "generic_seq3"
    assert seqs[0].description == "generic_seq1 seqdescription"
    assert seqs[1].description == "generic_seq2 seqdescription"
    assert seqs[2].description == "generic_seq3 seqdescription"
    assert seqs[0].seq == "ATCG"
    assert seqs[1].seq == "ATGC"
    assert seqs[2].seq == "CATG"


def test_fasta_duplicate_id(tmp_path):
    fasta_path = tmp_path / "fasta.fasta"
    fasta_path.write_text(">seq1\nATGC\n>seq1\nATAT\n")
    fasta = Fasta(fasta_path)
    with pytest.raises(ValueError):
        list(fasta.seqs)


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
    fasta = Fasta(path)
    with pytest.raises(ValueError):
        list(fasta.seqs)


def test_fasta_no_seq(tmp_path):
    path = tmp_path / "fasta_no_id.fasta"
    path.write_text(">seq1\n", encoding="utf-8")
    fasta = Fasta(path)
    with pytest.raises(ValueError):
        list(fasta.seqs)


def test_mixed_types_fasta(tmp_path):
    path = tmp_path / "fasta.fasta"
    path.write_text(">seq1 genome [protein_id=z]\nATGC\n>seq2 genome2\nATGG\n")
    fasta = Fasta(path)
    with pytest.raises(ValueError):
        list(fasta.seqs)


def test_genome_ids(protein_fasta):
    genome_ids = protein_fasta.genome_ids
    assert genome_ids[0] == "genome_1"
    assert genome_ids[1] == "genome_2"
    assert genome_ids[2] == "genome_3"


def test_genome_ids_in_generic_fasta(generic_fasta):
    with pytest.raises(AttributeError):
        _ = generic_fasta.genome_ids


def test_n_seqs(generic_fasta):
    assert generic_fasta.n_seqs == 3


def test_get_seqs_right_order(generic_fasta):
    seqs = generic_fasta.get_seqs("generic_seq2", "generic_seq1", "generic_seq3")
    assert len(seqs) == 3
    assert seqs[0].id == "generic_seq2"
    assert seqs[1].id == "generic_seq1"
    assert seqs[2].id == "generic_seq3"


def test_get_seqs_no_ids(generic_fasta):
    with pytest.raises(ValueError):
        _ = generic_fasta.get_seqs()


def test_get_seqs_duplicate_ids(generic_fasta):
    with pytest.raises(ValueError):
        _ = generic_fasta.get_seqs("generic_seq1", "generic_seq1")


def test_get_seqs_sequence_not_found(generic_fasta):
    with pytest.raises(ValueError):
        _ = generic_fasta.get_seqs("seq10")


def test_add_seq(generic_fasta):
    seq = Seq(
        seq=BioSeq("AT"),
        seq_id="seq4",
        seq_description="",
    )
    generic_fasta.add_seqs(seq)
    seqs = list(generic_fasta.seqs)
    assert seqs[3].id == "seq4"
    assert seqs[3].seq == "AT"


def test_add_no_sequences(generic_fasta):
    with pytest.raises(ValueError):
        generic_fasta.add_seqs()


def test_add_seqs_rejects_duplicate_ids(generic_fasta):
    seq1 = Seq(BioSeq("ATGC"), "sequence1", "sequence 1")
    seq2 = Seq(BioSeq("GGCC"), "sequence1", "sequence 2")
    with pytest.raises(ValueError):
        generic_fasta.add_seqs(seq1, seq2)


def test_add_seqs_rejects_existing_ids(generic_fasta):
    seq1 = Seq(BioSeq("ATGC"), "seq1", "sequence 1")
    generic_fasta.add_seqs(seq1)

    seq2 = Seq(BioSeq("GGCC"), "seq1", "sequence 2")
    with pytest.raises(ValueError):
        generic_fasta.add_seqs(seq2)


def test_add_seqs_rejects_different_types(generic_fasta):
    prot_seq = Seq(BioSeq("ATGC"), "seq1", "seq1 genome_id [protein_id=]")

    with pytest.raises(ValueError):
        generic_fasta.add_seqs(prot_seq)


def test_add_seqs_rejects_different_types_prot(protein_fasta):
    prot_seq = Seq(BioSeq("ATGC"), "seq1", "seq1 genome_id")

    with pytest.raises(ValueError):
        protein_fasta.add_seqs(prot_seq)


def test_add_seqs_rejects_different_types_in_same_call(generic_fasta):
    generic_seq = Seq(BioSeq("ATGC"), "seq1", "seq1 genome_1")
    protein_seq = Seq(BioSeq("AUGC"), "seq2", "seq2 genome_id [protein_id=x]")

    with pytest.raises(ValueError):
        generic_fasta.add_seqs(generic_seq, protein_seq)


def test_add_seqs_accepts_same_type(generic_fasta):
    seq1 = Seq(BioSeq("ATGC"), "seq1", "seq1 genome_id")
    seq2 = Seq(BioSeq("GGCC"), "seq2", "seq2 genome_id")

    generic_fasta.add_seqs(seq1, seq2)


def test_add_seqs_same_type_to_existing_fasta(protein_fasta):
    seq1 = Seq(BioSeq("ATGC"), "seq1", "seq1 genome_id [protein_id=x]")
    protein_fasta.add_seqs(seq1)


def test_remove_seqs(generic_fasta):
    assert len(list(generic_fasta.seqs)) == 3
    ids = generic_fasta.ids
    assert "generic_seq1" in ids
    generic_fasta.remove_seqs("generic_seq1")
    ids = generic_fasta.ids
    assert "generic_seq1" not in ids
    assert "generic_seq2" in ids
    assert "generic_seq3" in ids
    assert len(list(generic_fasta.seqs)) == 2


def test_remove_all_seqs(generic_fasta):
    generic_fasta.remove_seqs("generic_seq1", "generic_seq2", "generic_seq3")
    assert not generic_fasta.path.exists()


def test_remove_no_sequences(generic_fasta):
    with pytest.raises(ValueError):
        generic_fasta.remove_seqs()


def test_remove_seqs_not_found(generic_fasta):
    with pytest.raises(ValueError):
        generic_fasta.remove_seqs("seq10")


def test_remove_duplicate_seq_id(generic_fasta):
    with pytest.raises(ValueError):
        generic_fasta.remove_seqs("seq1", "seq1")


def test_move_fasta(generic_fasta, tmp_path):
    dir_to = tmp_path / "dir_to"
    dir_to.mkdir()
    original_path = generic_fasta.path
    generic_fasta.move_fasta(dir_to)
    assert generic_fasta.path == dir_to / original_path.name
    assert generic_fasta.path.exists()
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
    dir_to = tmp_path / "dir_to"
    with pytest.raises(FileNotFoundError):
        generic_fasta.move_fasta(dir_to)


def test_move_fasta_to_existing_file_path(generic_fasta, tmp_path):
    with pytest.raises(FileExistsError):
        generic_fasta.move_fasta(tmp_path)


def test_rename_fasta(generic_fasta):
    generic_fasta.rename_fasta("testfile.txt")
    assert generic_fasta.path.name == "testfile.txt"
    assert generic_fasta.path.is_file()


def test_rename_nonexistent_fasta():
    fasta = Fasta(Path("fasta.fasta"))
    with pytest.raises(FileNotFoundError):
        fasta.rename_fasta("testfile.txt")


def test_rename_empty_filename(generic_fasta):
    with pytest.raises(ValueError):
        generic_fasta.rename_fasta("")


def test_rename_filename_no_suffix(generic_fasta):
    with pytest.raises(ValueError):
        generic_fasta.rename_fasta("test")


def test_rename_empty_fasta(tmp_path):
    fasta_path = tmp_path / "fasta.fasta"
    fasta_path.touch()
    fasta = Fasta(fasta_path)
    with pytest.raises(ValueError):
        fasta.rename_fasta("test")


def test_rename_fasta_to_existent_file_path(generic_fasta, tmp_path):
    other_file_path = tmp_path / "file.file"
    other_file_path.touch()
    other_file_name = other_file_path.name
    with pytest.raises(FileExistsError):
        generic_fasta.rename_fasta(other_file_name)


def test_delete_fasta(generic_fasta):
    generic_fasta.delete_fasta()
    assert generic_fasta.path.exists() is False


def test_delete_nonexistent_fasta():
    fasta = Fasta(Path("fasta.fasta"))
    with pytest.raises(FileNotFoundError):
        fasta.delete_fasta()
