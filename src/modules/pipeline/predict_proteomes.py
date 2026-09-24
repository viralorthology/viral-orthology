import logging

from config.context import Context
from engines import predict_proteome
from models.fasta import Fasta

logger = logging.getLogger(__name__)


def predict_proteomes(ctx: Context) -> None:
    """Predict proteomes from genome sequences and save the resulting proteomes."""
    ctx.ui.show("Predicting proteomes...")

    annotated_proteomes_fasta = Fasta(ctx.paths.proteomes_fasta)
    predicted_proteomes_fasta = Fasta(ctx.paths.predicted_proteomes_fasta)
    genomes_fasta = Fasta(ctx.paths.genomes_fasta)
    for genome_seq in ctx.ui.progress_bar(
        genomes_fasta.seqs, total=genomes_fasta.n_seqs
    ):
        logger.info("Predicting %s proteome", genome_seq.id)

        proteome = [
            seq
            for seq in annotated_proteomes_fasta.seqs
            if seq.genome_id == genome_seq.id
        ] or None

        predicted_proteins = predict_proteome(
            genome_seq,
            ctx.args.tool_args["orffinder"],
            ctx.runtime.predicted_proteins_count,
            proteome,
        )

        logger.info(
            "%d proteins predicted for genome %s",
            len(predicted_proteins),
            genome_seq.id,
        )

        predicted_proteomes_fasta.add_seqs(*predicted_proteins)
        ctx.runtime.predicted_proteins_count += len(predicted_proteins)
