import pytest
from Bio.Seq import Seq as BioSeq

from models.ortholog_group import OrthologGroup
from models.seq import Seq


@pytest.fixture
def valid_ortholog_group(tmp_path):
    path = tmp_path / "valid_ortholog_group.fasta"
    path.write_text(
        ">seq1 genome1 [protein_id=123]\nATGC\n>seq2 genome2 [protein_id=234]\nATCG\n",
        encoding="utf-8",
    )
    return path


def test_empty_fasta(tmp_path):
    fasta_path = tmp_path / "fasta.fasta"
    fasta_path.touch()
    with pytest.raises(ValueError):
        _ = OrthologGroup(fasta_path)


def test_fasta_two_proteins_of_the_same_genome(tmp_path):
    path = tmp_path / "fasta_two_prots_from_the_same_genome.fasta"
    path.write_text(
        ">seq1 genome1 [protein_id=123]\nATGC\n>seq2 genome1 [protein_id=234]\nATCG\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError):
        _ = OrthologGroup(path)


def test_fasta_one_seq(tmp_path):
    path = tmp_path / "fasta_one_seq.fasta"
    path.write_text(">seq1 genome1 [protein_id=123]\nATGC\n", encoding="utf-8")

    with pytest.raises(ValueError):
        _ = OrthologGroup(path)


def test_delete_associated_files_with_add_seqs(valid_ortholog_group):
    og = OrthologGroup(valid_ortholog_group)
    associated_files = og.associated_files

    for file_ in associated_files:
        assert not file_.is_file()

    for file_ in associated_files:
        file_.touch()

    for file_ in associated_files:
        assert file_.is_file()

    og.add_seqs(
        Seq(
            BioSeq("ATGC"),
            "seq10",
            "seq10 genome10 [protein_id=123]",
        )
    )

    for file_ in associated_files:
        assert not file_.is_file()


def test_delete_associated_files_with_remove_seqs(
    valid_ortholog_group, context_always_yes
):
    og = OrthologGroup(valid_ortholog_group)
    associated_files = og.associated_files

    for file_ in associated_files:
        assert not file_.is_file()

    for file_ in associated_files:
        file_.touch()

    for file_ in associated_files:
        assert file_.is_file()

    og.remove_seqs("seq1", ctx=context_always_yes)

    for file_ in associated_files:
        assert not file_.is_file()


def test_delete_associated_files_with_delete_fasta(valid_ortholog_group):
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


def test_delete_associated_files_with_rename_fasta(valid_ortholog_group):
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


def test_delete_one_seq_in_two_seqs_fasta_annotated_db(
    valid_ortholog_group, context_always_yes
):
    og = OrthologGroup(valid_ortholog_group)
    og.remove_seqs("seq1", ctx=context_always_yes)

    assert not og.path.is_file()
    assert not context_always_yes.paths.predicted_prots_db.is_file()
    assert (
        context_always_yes.paths.annotated_prots_db.read_text(encoding="utf-8")
        == ">seq2 genome2 [protein_id=234]\nATCG\n"
    )


def test_delete_one_seq_in_two_seqs_fasta_predicted_db(tmp_path, context_always_yes):
    path = tmp_path / "valid_ortholog_group_with_predicted_seq.fasta"
    path.write_text(
        ">seq1 genome1 [protein_id=123]\nATGC\n>PREDICTED_10 genome2 [protein_id=234]\nATCG\n",
        encoding="utf-8",
    )

    og = OrthologGroup(path)
    og.remove_seqs("seq1", ctx=context_always_yes)

    assert not og.path.is_file()
    assert not context_always_yes.paths.annotated_prots_db.is_file()
    assert (
        context_always_yes.paths.predicted_prots_db.read_text(encoding="utf-8")
        == ">PREDICTED_10 genome2 [protein_id=234]\nATCG\n"
    )


def test_delete_one_seq(tmp_path, context_always_yes):
    path = tmp_path / "valid_ortholog_group.fasta"
    path.write_text(
        ">seq1 genome1 [protein_id=123]\nATGC\n>seq2 genome2 [protein_id=234]\nATCG\n>seq3 genome3 [protein_id=345]\nATGG\n",
        encoding="utf-8",
    )
    og = OrthologGroup(path)
    og.remove_seqs("seq1", ctx=context_always_yes)

    assert og.path.is_file()


def test_delete_missing_seq(valid_ortholog_group, context_always_yes):
    og = OrthologGroup(valid_ortholog_group)
    with pytest.raises(ValueError):
        og.remove_seqs("seq1", "seq10", ctx=context_always_yes)


def test_remove_seqs_context_none(valid_ortholog_group):
    og = OrthologGroup(valid_ortholog_group)
    with pytest.raises(ValueError):
        og.remove_seqs("seq1", "seq10")
