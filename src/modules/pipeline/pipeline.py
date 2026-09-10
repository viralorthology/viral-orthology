import utils
from config.context import Context
from models.fasta import Fasta
from modules.pipeline.predict_proteomes import predict_proteomes
from modules.pipeline.split_fastas import split_fastas


def run(ctx: Context) -> None:
    ### VALIDATION ###
    utils.check_dependencies(
        "blastn",
        "blastp",
        "ORFfinder",
        "proteinortho",
        "proteinortho_grab_proteins.pl",
        "hmmbuild",
        "muscle",
        "hhmake",
        "hhalign",
    )
    utils.ensure_files_have_content(ctx.paths.genomes_fasta, ctx.paths.proteomes_fasta)
    utils.ensure_paths_do_not_exist(
        ctx.paths.annotated_prots_db,
        ctx.paths.predicted_prots_db,
        ctx.paths.ortholog_groups_dir,
        ctx.paths.sequences_dir,
        ctx.paths.redundant_seqs_dir,
        ctx.paths.paralogs_dir,
    )

    genome_ids = Fasta(ctx.paths.genomes_fasta).ids
    genome_ids_from_proteomes = Fasta(ctx.paths.proteomes_fasta).genome_ids

    if len(genome_ids) < 2:
        raise ValueError(
            "ViralOrthology requires at least two different genome sequences to continue."
        )

    proteomes_without_genome = set(genome_ids_from_proteomes) - set(genome_ids)
    if proteomes_without_genome:
        raise ValueError(
            f"There are proteomes without their corresponding genome sequence: {proteomes_without_genome}"
        )

    genomes_without_proteome = set(genome_ids) - set(genome_ids_from_proteomes)
    if genomes_without_proteome:
        raise ValueError(
            f"There are genomes without their corresponding proteome: {genomes_without_proteome}"
        )

    ### RUN ###
    _preparation_stage(ctx)

    ### MAKE REPORTS ###


def _preparation_stage(ctx: Context) -> None:
    split_fastas(ctx)
    predict_proteomes(ctx)
    # find paralogs
    # find redundant genomes
