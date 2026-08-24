import pytest
from Bio.Seq import Seq as BioSeq
from Bio.SeqRecord import SeqRecord

from classes.fasta_type import FastaType
from classes.seq import Seq


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
    seqs = fasta.get_seqs("seq1", "seq2", "seq3")
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
    seq = Seq(SeqRecord(BioSeq("AT"), id="seq4"), FastaType.GENERIC)
    fasta.add_seqs(seq)
    seqs = list(fasta.seqs)
    assert seqs[3].id == "seq4"
    assert seqs[3].seq == "AT"


def test_n_seqs(fasta):
    assert fasta.n_seqs == 3
