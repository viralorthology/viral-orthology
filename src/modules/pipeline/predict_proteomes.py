from config.context import Context
from engines import predict_proteome
from models.fasta import Fasta


def predict_proteomes(ctx: Context) -> None:
    """Predict proteomes from genome sequences and save the resulting proteomes."""
    ctx.ui.show("Predicting proteomes...")

    predicted_proteomes_fasta = Fasta(ctx.paths.predicted_proteomes_fasta)
    genomes_fasta = Fasta(ctx.paths.genomes_fasta)
    for genome_seq in ctx.ui.progress_bar(
        genomes_fasta.seqs, total=genomes_fasta.n_seqs
    ):
        proteome_fasta = Fasta(
            ctx.paths.proteomes_dir / f"{genome_seq.id}.fasta"
        )  # TODO it may not exist
        predicted_proteins = predict_proteome(
            genome_seq,
            ctx.args.tool_args["orffinder"],
            ctx.runtime.predicted_proteins_count,
            proteome_fasta,
        )
        predicted_proteomes_fasta.add_seqs(*predicted_proteins)
        ctx.runtime.predicted_proteins_count += len(predicted_proteins)
