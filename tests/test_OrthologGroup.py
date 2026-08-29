import pytest

from models.ortholog_group import OrthologGroup


def test_fasta_two_proteins_of_the_same_genome(
    fasta_two_prots_from_the_same_genome,
):
    with pytest.raises(ValueError):
        _ = OrthologGroup(fasta_two_prots_from_the_same_genome)


def test_fasta_one_seq(fasta_one_seq):
    with pytest.raises(ValueError):
        _ = OrthologGroup(fasta_one_seq)


def test_associated_files(valid_ortholog_group):
    og = OrthologGroup(valid_ortholog_group)
    assert og.alignment_path == og.path.with_suffix(".muscle")
    assert og.hmm_hmmer_path == og.path.with_suffix(".hmm")
    assert og.hmm_hhsuite_path == og.path.with_suffix(".hhm")
    assert og.a2m_path == og.path.with_suffix(".a2m")


def test_delete_associated_files_with_delete_file(valid_ortholog_group):
    og = OrthologGroup(valid_ortholog_group)
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
    og = OrthologGroup(valid_ortholog_group)
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


def test_delete_one_seq_in_two_seqs_fasta_annotated_db(valid_ortholog_group, context):
    og = OrthologGroup(valid_ortholog_group)
    og.remove_seqs("seq1", ctx=context)

    assert not og.path.is_file()
    assert not context.paths.predicted_prots_db.is_file()
    assert (
        context.paths.annotated_prots_db.read_text(encoding="utf-8")
        == ">seq2 genome2\nATCG\n"
    )


def test_delete_one_seq_in_two_seqs_fasta_predicted_db(
    valid_ortholog_group_with_orffinder_seq, context
):
    og = OrthologGroup(valid_ortholog_group_with_orffinder_seq)
    og.remove_seqs("seq1", ctx=context)

    assert not og.path.is_file()
    assert not context.paths.annotated_prots_db.is_file()
    assert (
        context.paths.predicted_prots_db.read_text(encoding="utf-8")
        == ">ORFFINDER1 genome2\nATCG\n"
    )


def test_delete_one_seq(valid_ortholog_group_three_seqs, context):
    og = OrthologGroup(valid_ortholog_group_three_seqs)
    og.remove_seqs("seq1", ctx=context)

    assert og.path.is_file()


def test_delete_missing_seq(valid_ortholog_group, context):
    og = OrthologGroup(valid_ortholog_group)
    with pytest.raises(ValueError):
        og.remove_seqs("seq1", "seq10", ctx=context)
