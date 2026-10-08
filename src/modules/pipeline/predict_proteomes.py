import logging

import engines
from config.context import Context
from models.fasta import Fasta

logger = logging.getLogger(__name__)


def predict_proteomes(ctx: Context, genome_ids: set[str]) -> None:
    """Predict proteomes from genome sequences and save the resulting proteomes."""
    ctx.ui.show("Predicting proteomes...")

    proteomes_fasta = Fasta(ctx.paths.proteomes_fasta)
    predicted_unique_prots_fasta = Fasta(ctx.paths.predicted_unique_prots_fasta)
    genomes_fasta = Fasta(ctx.paths.genomes_fasta)

    for genome_id in ctx.ui.progress_bar(sorted(genome_ids)):
        genome_seq = genomes_fasta.get_seqs(genome_id)[0]
        proteome = [
            seq for seq in proteomes_fasta.seqs if seq.genome_id == genome_id
        ] or None

        predicted_proteins = engines.predict_proteome(
            genome_seq,
            ctx.args.tool_args["orffinder"],
            ctx.runtime.predicted_proteins_count,
            proteome,
        )

        logger.info(
            "%d proteins predicted for genome %s",
            len(predicted_proteins),
            genome_id,
        )

        predicted_unique_prots_fasta.add_seqs(*predicted_proteins)
        ctx.runtime.predicted_proteins_count += len(predicted_proteins)
