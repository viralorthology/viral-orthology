from models.tmp_fasta import TmpFasta


def test_tmp_fasta_cleans_up_after_context():
    with TmpFasta() as fasta:
        fasta_path = fasta.path
        assert not fasta_path.exists()
        fasta_path.write_text("test")
        assert fasta_path.exists()

    assert not fasta_path.exists()
