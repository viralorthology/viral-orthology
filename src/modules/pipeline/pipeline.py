import utils
from config.context import Context
from models.fasta import Fasta
from modules.pipeline.find_paralogs import find_paralogs
from modules.pipeline.find_redundant_genomes import find_redundant_genomes
from modules.pipeline.make_initial_ortholog_groups import make_initial_ortholog_groups
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
        ctx.paths.annotated_unique_prots_fasta,
        ctx.paths.predicted_unique_prots_fasta,
        ctx.paths.ortholog_groups_dir,
        ctx.paths.proteomes_dir,
        ctx.paths.paralogs_dir,
        ctx.paths.predicted_proteomes_fasta,
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

    ### RUN ###
    _preparation_stage(ctx)
    _initial_stage(ctx)

    ### MAKE REPORTS ###


def _preparation_stage(ctx: Context) -> None:
    split_fastas(ctx)
    predict_proteomes(ctx)
    find_paralogs(ctx)
    find_redundant_genomes(ctx)


def _initial_stage(ctx: Context) -> None:
    make_initial_ortholog_groups(ctx)
