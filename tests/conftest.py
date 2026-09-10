from collections.abc import Iterable

import pytest

from cli.args import Args
from cli.ui import UI
from config.context import Context
from config.paths import Paths

# CONTEXT


@pytest.fixture
def ui_always_yes():
    class UIAlwaysYes(UI):
        def show(self, text: str) -> None:
            pass

        def show_error(self, text: str) -> None:
            pass

        def progress_bar[T](self, iterable: Iterable[T]) -> Iterable[T]:
            return iterable

        def ask_yes_no(self, question: str) -> bool:
            return True

    return UIAlwaysYes()


@pytest.fixture
def context_always_yes(tmp_path, ui_always_yes):
    return Context(
        Args(debug=False, assume_yes=False, tool_args={}),
        Paths(tmp_path),
        ui_always_yes,
    )


# FASTA FILES


@pytest.fixture
def fasta(tmp_path):
    path = tmp_path / "fasta.fasta"
    path.write_text(
        ">seq1 seqdescription\nATCG\n>seq2 seqdescription\nATGC\n>seq3 seqdescription\nCATG\n",
        encoding="utf-8",
    )
    return path


@pytest.fixture
def nonexistent_fasta(tmp_path):
    path = tmp_path / "nonexistent_fasta.fasta"
    return path


@pytest.fixture
def empty_fasta(tmp_path):
    path = tmp_path / "empty_fasta.fasta"
    path.write_text("", encoding="utf-8")
    return path


@pytest.fixture
def fasta_duplicate_id(tmp_path):
    path = tmp_path / "fasta_duplicate_id.fasta"
    path.write_text(
        ">seqid seqdescription\nATCG\n>seqid seqdescription\nATCG\n", encoding="utf-8"
    )
    return path


@pytest.fixture
def fasta_no_id(tmp_path):
    path = tmp_path / "fasta_no_id.fasta"
    path.write_text(">\nATCG\n", encoding="utf-8")
    return path


@pytest.fixture
def fasta_no_seqs(tmp_path):
    path = tmp_path / "fasta_no_id.fasta"
    path.write_text(">seq1\n", encoding="utf-8")
    return path


@pytest.fixture
def protein_fasta(tmp_path):
    path = tmp_path / "protein_fasta.fasta"
    path.write_text(
        ">seq1 genome1 [protein_id=123]\nATGC\n>seq2 genome2 [protein_id=234]\nATCG\n",
        encoding="utf-8",
    )
    return path


# ORTHOLOG GROUP FASTA FILES


@pytest.fixture
def fasta_two_prots_from_the_same_genome(tmp_path):
    path = tmp_path / "fasta_two_prots_from_the_same_genome.fasta"
    path.write_text(
        ">seq1 genome1 [protein_id=123]\nATGC\n>seq2 genome1 [protein_id=234]\nATCG\n",
        encoding="utf-8",
    )
    return path


@pytest.fixture
def fasta_one_seq(tmp_path):
    path = tmp_path / "fasta_one_seq.fasta"
    path.write_text(">seq1 genome1 [protein_id=123]\nATGC\n", encoding="utf-8")
    return path


@pytest.fixture
def valid_ortholog_group(tmp_path):
    path = tmp_path / "valid_ortholog_group.fasta"
    path.write_text(
        ">seq1 genome1 [protein_id=123]\nATGC\n>seq2 genome2 [protein_id=234]\nATCG\n",
        encoding="utf-8",
    )
    return path


@pytest.fixture
def valid_ortholog_group_three_seqs(tmp_path):
    path = tmp_path / "valid_ortholog_group.fasta"
    path.write_text(
        ">seq1 genome1 [protein_id=123]\nATGC\n>seq2 genome2 [protein_id=234]\nATCG\n>seq3 genome3 [protein_id=345]\nATGG\n",
        encoding="utf-8",
    )
    return path


@pytest.fixture
def valid_ortholog_group_with_predicted_seq(tmp_path):
    path = tmp_path / "valid_ortholog_group_with_predicted_seq.fasta"
    path.write_text(
        ">seq1 genome1 [protein_id=123]\nATGC\n>PREDICTED_10 genome2 [protein_id=234]\nATCG\n",
        encoding="utf-8",
    )
    return path
