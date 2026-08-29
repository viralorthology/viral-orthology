from pathlib import Path

import pytest

from models.ortholog_group import OrthologGroup


def test_fasta_two_proteins_of_the_same_genome(
    fasta_two_prots_from_the_same_genome,
):
    with pytest.raises(ValueError):
        _ = OrthologGroup(
            fasta_two_prots_from_the_same_genome, Path("test"), Path("test2")
        )


def test_fasta_one_seq(fasta_one_seq):
    with pytest.raises(ValueError):
        _ = OrthologGroup(fasta_one_seq, Path("test"), Path("test2"))


def test_associated_files(valid_ortholog_group):
    og = OrthologGroup(valid_ortholog_group, Path("test"), Path("test2"))
    assert og.alignment_path == og.path.with_suffix(".muscle")
    assert og.hmm_hmmer_path == og.path.with_suffix(".hmm")
    assert og.hmm_hhsuite_path == og.path.with_suffix(".hhm")
    assert og.a2m_path == og.path.with_suffix(".a2m")


def test_delete_associated_files_with_delete_file(valid_ortholog_group):
    og = OrthologGroup(valid_ortholog_group, Path("test"), Path("test2"))
    associated_files = og.associated_files

    for file_ in associated_files:
        assert not file_.is_file()

    for file_ in associated_files:
        file_.touch()

    for file_ in associated_files:
        assert file_.is_file()

    og.delete_fasta()

    for file_ in associated_files:
        assert not file_.is_file()


def test_delete_associated_files_with_rename_file(valid_ortholog_group):
    og = OrthologGroup(valid_ortholog_group, Path("test"), Path("test2"))
    associated_files = og.associated_files

    for file_ in associated_files:
        assert not file_.is_file()

    for file_ in associated_files:
        file_.touch()

    for file_ in associated_files:
        assert file_.is_file()

    og.rename_fasta("test.fasta")

    for file_ in associated_files:
        assert not file_.is_file()


def test_delete_one_seq_in_two_seqs_fasta_annotated_db(valid_ortholog_group):
    annotated_db_path = valid_ortholog_group.with_name("annotated.db")
    orffinder_db_path = valid_ortholog_group.with_name("orffinder.db")
    og = OrthologGroup(valid_ortholog_group, annotated_db_path, orffinder_db_path)
    og.remove_seqs("seq1")

    assert not og.path.is_file()
    assert not orffinder_db_path.is_file()
    assert annotated_db_path.read_text(encoding="utf-8") == ">seq2 genome2\nATCG\n"


def test_delete_one_seq_in_two_seqs_fasta_orffinder_db(
    valid_ortholog_group_with_orffinder_seq,
):
    annotated_db_path = valid_ortholog_group_with_orffinder_seq.with_name(
        "annotated.db"
    )
    orffinder_db_path = valid_ortholog_group_with_orffinder_seq.with_name(
        "orffinder.db"
    )
    og = OrthologGroup(
        valid_ortholog_group_with_orffinder_seq, annotated_db_path, orffinder_db_path
    )
    og.remove_seqs("seq1")

    assert not og.path.is_file()
    assert not annotated_db_path.is_file()
    assert (
        orffinder_db_path.read_text(encoding="utf-8") == ">ORFFINDER1 genome2\nATCG\n"
    )


def test_delete_one_seq(
    valid_ortholog_group_three_seqs,
):
    og = OrthologGroup(valid_ortholog_group_three_seqs, Path("test"), Path("test2"))
    og.remove_seqs("seq1")

    assert og.path.is_file()


def test_delete_missing_seq(valid_ortholog_group):
    annotated_db_path = valid_ortholog_group.with_name("annotated.db")
    orffinder_db_path = valid_ortholog_group.with_name("orffinder.db")
    og = OrthologGroup(valid_ortholog_group, annotated_db_path, orffinder_db_path)
    with pytest.raises(ValueError):
        og.remove_seqs("seq1", "seq10")
