import shutil
import tempfile
from pathlib import Path

import utils
from config.context import Context
from engines.blast import make_blast_db, run_blastn
from models.fasta import Fasta


def find_redundant_genomes(ctx: Context) -> None:
    """
    Identify and remove redundant genomes from the first analysis round.

    Genomes that are sufficiently similar to another genome are identified
    and moved, along with their associated proteome files, to the redundant
    sequences directory. Redundant genomes are only removed if at least two
    genomes remain for the first analysis round.
    """
    ctx.ui.show("Searching for redundant genomes...")

    genome_fastas = utils.get_fastas(ctx.paths.sequences_dir, ".genome")
    genome_fastas.sort(key=lambda fasta: fasta.path.stem)

    redundant_genomes = _get_redundant_genomes(ctx, genome_fastas)

    n_genomes_for_first_round = len(genome_fastas) - len(redundant_genomes)
    if (  # the pipeline needs at least 2 genomes on first round, else dont filter
        n_genomes_for_first_round >= 2 and redundant_genomes
    ):
        # move redundant genomes and proteomes
        ctx.paths.redundant_seqs_dir.mkdir()
        for genome_fasta in redundant_genomes:
            proteome_fasta = Fasta(
                ctx.paths.sequences_dir / f"{genome_fasta.path.stem}.proteome"
            )
            predicted_proteome_fasta = Fasta(
                ctx.paths.sequences_dir / f"{genome_fasta.path.stem}.predicted"
            )
            genome_fasta.move_fasta(ctx.paths.redundant_seqs_dir)
            proteome_fasta.move_fasta(ctx.paths.redundant_seqs_dir)
            predicted_proteome_fasta.move_fasta(ctx.paths.redundant_seqs_dir)


def _get_redundant_genomes(ctx: Context, genome_fastas: list[Fasta]) -> list[Fasta]:
    """
    Identify redundant genomes based on pairwise BLASTN similarity.

    Genome files are copied to a temporary directory where BLAST databases
    are created for the similarity comparisons. The identified redundant
    genomes are mapped back to their original FASTA files.

    Args:
        ctx: Pipeline context used for progress reporting.
        genome_fastas: Genomes to compare.

    Returns:
        The original FASTA files identified as redundant.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir_path = Path(tmp_dir)

        original_fasta_by_tmp_path = {}
        tmp_fastas = []
        for fasta in genome_fastas:
            fasta_tmp_path = tmp_dir_path / fasta.path.name
            original_fasta_by_tmp_path[fasta_tmp_path] = fasta
            shutil.copy(fasta.path, fasta_tmp_path)
            tmp_fasta = Fasta(fasta_tmp_path)
            tmp_fastas.append(tmp_fasta)

        for genome_fasta in tmp_fastas:
            make_blast_db(genome_fasta, "nucl")

        redundant_genomes = _find_redundant_genomes(ctx, tmp_fastas)

        return [
            original_fasta_by_tmp_path[redundant.path]
            for redundant in redundant_genomes
        ]


def _find_redundant_genomes(ctx: Context, genome_fastas: list[Fasta]) -> set[Fasta]:
    """
    Identify redundant genomes based on pairwise sequence similarity.

    Genomes are compared pairwise and a genome is considered redundant when
    it is sufficiently similar to an earlier genome in the input list.

    Args:
        ctx: Pipeline context used for progress reporting.
        genome_fastas: Genomes to compare.

    Returns:
        The genomes identified as redundant.
    """
    redundant_genomes = set()

    for i, g1 in enumerate(ctx.ui.progress_bar(genome_fastas)):
        if g1 in redundant_genomes:
            continue

        for g2 in genome_fastas[i + 1 :]:
            if g2 in redundant_genomes:
                continue

            if _genomes_are_similar(g1, g2, 90, 95):
                redundant_genomes.add(g2)

    return redundant_genomes


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
