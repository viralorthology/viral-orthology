import logging

import utils
from config.context import Context
from models.fasta import Fasta
from modules.pipeline.blastp_search import blastp_search
from modules.pipeline.clean_ortholog_groups import clean_ortholog_groups
from modules.pipeline.find_paralogs import find_paralogs
from modules.pipeline.find_redundant_genomes import find_redundant_genomes
from modules.pipeline.hmm_search import hmm_search
from modules.pipeline.make_annotated_unique_prots_fasta import (
    make_annotated_unique_prots_fasta,
)
from modules.pipeline.make_initial_ortholog_groups import make_initial_ortholog_groups
from modules.pipeline.merge_by_annotation import merge_by_annotation
from modules.pipeline.predict_proteomes import predict_proteomes
from modules.pipeline.rename_og_fastas_by_annotation import (
    rename_og_fastas_by_annotation,
)

logger = logging.getLogger(__name__)


def run(ctx: Context) -> None:
    logger.info("Running pipeline")

    ### PRE-RUN ###
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

    logger.info("Validation completed successfully")

    ### RUN ###
    _preparation_stage(ctx)
    _initial_stage(ctx)
    _enrichment_stage(ctx)


def _preparation_stage(ctx: Context) -> None:
    find_paralogs(ctx)
    find_redundant_genomes(ctx)


def _initial_stage(ctx: Context) -> None:
    predict_proteomes(ctx, ctx.runtime.active_genome_ids)
    make_initial_ortholog_groups(ctx)
    rename_og_fastas_by_annotation(ctx)
    clean_ortholog_groups(ctx)
    make_annotated_unique_prots_fasta(ctx)


def _enrichment_stage(ctx: Context) -> None:
    merge_by_annotation(ctx)
    hmm_search(ctx)
    blastp_search(ctx, ctx.paths.annotated_unique_prots_fasta)
    hmm_search(ctx)
    blastp_search(ctx, ctx.paths.predicted_unique_prots_fasta)
    hmm_search(ctx)
