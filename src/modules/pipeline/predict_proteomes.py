import utils
from config.context import Context
from engines import predict_proteome
from models.fasta import Fasta


def predict_proteomes(ctx: Context) -> None:
    """Predict proteomes from genome sequences and save the resulting proteomes."""
    ctx.ui.show("Predicting proteomes...")

    genome_fastas = utils.get_fastas(ctx.paths.sequences_dir, ".genome")
    n_predicted_proteins = ctx.runtime.predicted_proteins_count
    for genome_fasta in ctx.ui.progress_bar(genome_fastas):
        genome_id = genome_fasta.ids[0]
        proteome_fasta = Fasta(ctx.paths.sequences_dir / f"{genome_id}.proteome")
        predicted_proteins = predict_proteome(
            genome_fasta,
            ctx.args.tool_args["orffinder"],
            n_predicted_proteins,
            proteome_fasta,
        )
        predicted_fasta = Fasta(ctx.paths.sequences_dir / f"{genome_id}.predicted")
        predicted_fasta.add_seqs(*predicted_proteins)
        ctx.runtime.predicted_proteins_count += len(predicted_proteins)
