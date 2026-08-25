import pytest
from Bio.Seq import Seq as BioSeq
from Bio.SeqRecord import SeqRecord

from classes.fasta import Fasta
from classes.fasta_type import FastaType
from classes.file import File

# TEXT FILES


@pytest.fixture
def nonexistent_file(tmp_path):
    path = tmp_path / "this_file_does_not_exist.txt"
    return File(path)


@pytest.fixture
def empty_file(tmp_path):
    path = tmp_path / "empty_file.txt"
    path.write_text("", encoding="utf-8")
    return File(path)


@pytest.fixture
def text_file(tmp_path):
    path = tmp_path / "text_file.txt"
    path.write_text("sometext", encoding="utf-8")
    return File(path)


# FASTA FILES


@pytest.fixture
def fasta(tmp_path):
    path = tmp_path / "fasta.fasta"
    path.write_text(
        ">seq1 seqdescription\nATCG\n>seq2 seqdescription\nATGC\n>seq3 seqdescription\nCATG\n",
        encoding="utf-8",
    )
    return Fasta(path, FastaType.GENERIC)


@pytest.fixture
def fasta_duplicate_id(tmp_path):
    path = tmp_path / "fasta_duplicate_id.fasta"
    path.write_text(
        ">seqid seqdescription\nATCG\n>seqid seqdescription\nATCG\n", encoding="utf-8"
    )
    return Fasta(path, FastaType.GENERIC)


@pytest.fixture
def fasta_no_id(tmp_path):
    path = tmp_path / "fasta_no_id.fasta"
    path.write_text(">\nATCG\n", encoding="utf-8")
    return Fasta(path, FastaType.GENERIC)


@pytest.fixture
def nonexistent_fasta(tmp_path):
    path = tmp_path / "nonexistent_fasta.fasta"
    return Fasta(path, FastaType.GENERIC)


@pytest.fixture
def empty_fasta(tmp_path):
    path = tmp_path / "empty_fasta.fasta"
    path.write_text("", encoding="utf-8")
    return Fasta(path, FastaType.GENERIC)


@pytest.fixture
def protein_fasta(tmp_path):
    path = tmp_path / "protein_fasta.fasta"
    path.write_text(">seq1 genome1\nATGC\n>seq2 genome2\nATCG\n", encoding="utf-8")
    return Fasta(path, FastaType.PROTEIN)


# SEQ


@pytest.fixture
def seq_with_genome_id():
    return SeqRecord(BioSeq("ATGC"), id="seq1", description="seq1 genome1")


@pytest.fixture
def seq_without_genome_id():
    return SeqRecord(BioSeq("ATGC"), id="seq1", description="seq1")


# ORTHOLOG GROUP


@pytest.fixture
def fasta_two_prots_from_the_same_genome(tmp_path):
    path = tmp_path / "fasta_two_prots_from_the_same_genome.fasta"
    path.write_text(">seq1 genome1\nATGC\n>seq2 genome1\nATCG\n", encoding="utf-8")
    return path


@pytest.fixture
def fasta_one_seq(tmp_path):
    path = tmp_path / "fasta_one_seq.fasta"
    path.write_text(">seq1 genome1\nATGC\n", encoding="utf-8")
    return path


@pytest.fixture
def valid_ortholog_group(tmp_path):
    path = tmp_path / "valid_ortholog_group.fasta"
    path.write_text(">seq1 genome1\nATGC\n>seq2 genome2\nATCG\n", encoding="utf-8")
    return path


@pytest.fixture
def valid_ortholog_group_three_seqs(tmp_path):
    path = tmp_path / "valid_ortholog_group.fasta"
    path.write_text(
        ">seq1 genome1\nATGC\n>seq2 genome2\nATCG\n>seq3 genome3\nATGG\n",
        encoding="utf-8",
    )
    return path


@pytest.fixture
def valid_ortholog_group_with_orffinder_seq(tmp_path):
    path = tmp_path / "valid_ortholog_group_with_orffinder_seq.fasta"
    path.write_text(
        ">seq1 genome1\nATGC\n>ORFFINDER1 genome2\nATCG\n", encoding="utf-8"
    )
    return path
