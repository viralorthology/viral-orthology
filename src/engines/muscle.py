import utils
from models.ortholog_group import OrthologGroup


def muscle_align(fasta: OrthologGroup) -> None:
    """
    Create a multiple sequence alignment using MUSCLE.

    The -super5 algorithm is used for ortholog groups with more than
    500 sequences; otherwise, the standard alignment algorithm is used.
    """
    align_algorithm = "-super5" if fasta.n_seqs > 500 else "-align"

    utils.run_cmd(
        f"muscle {align_algorithm} {fasta.path} -output {fasta.alignment_path}"
    )
