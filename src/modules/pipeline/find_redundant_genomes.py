import tempfile
from pathlib import Path

from config.context import Context
from engines.blast import make_blast_db, run_blastn
from models.fasta import Fasta


def find_redundant_genomes(ctx: Context) -> None:
    """
    Identify and remove redundant genomes from the first analysis round.

    Redundant genomes are only filtered if at least two genomes remain for
    the first analysis round.
    """
    ctx.ui.show("Searching for redundant genomes...")

    genomes_fasta = Fasta(ctx.paths.genomes_fasta)
    redundant_genome_ids = _get_redundant_genomes(ctx, genomes_fasta)

    n_genomes_for_first_round = genomes_fasta.n_seqs - len(redundant_genome_ids)
    if (  # the pipeline needs at least 2 genomes on first round, else dont filter
        n_genomes_for_first_round >= 2 and redundant_genome_ids
    ):
        ctx.runtime.redundant_genomes = True
        ctx.runtime.redundant_genome_ids = redundant_genome_ids


def _get_redundant_genomes(ctx: Context, genomes_fasta: Fasta) -> set[str]:
    """
    Identify redundant genomes based on pairwise BLASTN similarity.

    Args:
        ctx: Pipeline context used for progress reporting.
        genomes_fasta: Fasta file with all the available genomes.

    Returns:
        The Seq IDs of the redundant genomes.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir_path = Path(tmp_dir)

        genome_fastas = []
        for genome_seq in genomes_fasta.seqs:
            genome_fasta = Fasta(tmp_dir_path / f"{genome_seq.id}.fasta")
            genome_fasta.add_seqs(genome_seq)
            make_blast_db(genome_fasta, "nucl")
            genome_fastas.append(genome_fasta)

        genome_fastas.sort(key=lambda fasta: fasta.path.name)  # Ensure reproducibility
        redundant_genome_ids = _find_redundant_genomes(ctx, genome_fastas)

        return redundant_genome_ids


def _find_redundant_genomes(ctx: Context, genome_fastas: list[Fasta]) -> set[str]:
    """
    Identify redundant genomes based on pairwise sequence similarity.

    Genomes are compared pairwise and a genome is considered redundant when
    it is sufficiently similar to an earlier genome in the input list.

    Args:
        ctx: Pipeline context used for progress reporting.
        genome_fastas: Genomes to compare.

    Returns:
        The Seq IDs of the redundant genomes.
    """
    redundant_genome_ids = set()

    for i, g1 in enumerate(ctx.ui.progress_bar(genome_fastas)):
        if g1.path.stem in redundant_genome_ids:  # genome id = file stem
            continue

        for g2 in genome_fastas[i + 1 :]:
            if g2.path.stem in redundant_genome_ids:
                continue

            if _genomes_are_similar(g1, g2, 90, 95):
                redundant_genome_ids.add(g2.path.stem)

    return redundant_genome_ids


def _genomes_are_similar(
    query_fasta: Fasta,
    subject_fasta: Fasta,
    qcov: int,
    ident: int,
) -> bool:
    """
    Check whether two genomes have a sufficiently similar BLASTN hit.

    Args:
        query_fasta: FASTA file containing the query genome.
        subject_fasta: FASTA file containing the subject genome.
        qcov: Minimum query coverage threshold.
        ident: Minimum sequence identity threshold.

    Returns:
        True if at least one BLASTN hit exceeds both thresholds, otherwise False.
    """
    blastn_hits = run_blastn(query_fasta, subject_fasta, "")

    return any(hit.qcov > qcov and hit.ident > ident for hit in blastn_hits)
