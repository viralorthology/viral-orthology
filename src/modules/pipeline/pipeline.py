import logging
import sys

import utils
from config.context import Context
from models.fasta import Fasta
from modules.pipeline.blastp_search import blastp_search
from modules.pipeline.clean_ortholog_groups import clean_ortholog_groups
from modules.pipeline.find_paralogs import find_paralogs
from modules.pipeline.find_redundant_genomes import find_redundant_genomes
from modules.pipeline.hmm_search import hmm_search
from modules.pipeline.make_initial_ortholog_groups import make_initial_ortholog_groups
from modules.pipeline.manage_unique_prots_fasta import (
    add_redundant_genome_prots_unique_prots_fasta,
    make_annotated_unique_prots_fasta,
)
from modules.pipeline.merge_by_annotation import merge_by_annotation
from modules.pipeline.merge_by_hmm import merge_by_hmm
from modules.pipeline.paralogs_hmm_search import paralogs_hmm_search
from modules.pipeline.predict_proteomes import predict_proteomes
from modules.pipeline.rename_og_fastas_by_annotation import (
    rename_og_fastas_by_annotation,
)

logger = logging.getLogger(__name__)


def run(ctx: Context) -> None:
    logger.info("Running pipeline: %s", (" ").join(sys.argv))
    logger.info("Pipeline parameters:")
    for tool, params in ctx.args.tool_args.items():
        logger.info(" %s: %s", tool, params if params else "tool defaults")
    logger.info(" assume yes: %s", ctx.args.assume_yes)

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
    dataset_hash = utils.get_dataset_hash(genome_ids)
    logger.info("Starting pipeline on dataset %s", dataset_hash)

    ### RUN ###
    # redundancy reduction stage
    find_redundant_genomes(ctx)

    # initial stage
    find_paralogs(ctx, ctx.runtime.active_genome_ids)
    predict_proteomes(ctx, ctx.runtime.active_genome_ids)
    make_initial_ortholog_groups(ctx)
    rename_og_fastas_by_annotation(ctx)
    clean_ortholog_groups(ctx)
    merge_by_annotation(ctx)
    merge_by_hmm(ctx)
    make_annotated_unique_prots_fasta(ctx)

    # enrichment stage
    hmm_search(ctx)
    blastp_search(ctx, ctx.paths.annotated_unique_prots_fasta)
    hmm_search(ctx)
    blastp_search(ctx, ctx.paths.predicted_unique_prots_fasta)
    hmm_search(ctx)

    if ctx.runtime.redundant_genome_ids:
        ctx.runtime.active_genome_ids.update(ctx.runtime.redundant_genome_ids)
        add_redundant_genome_prots_unique_prots_fasta(ctx)
        find_paralogs(ctx, ctx.runtime.redundant_genome_ids)
        predict_proteomes(ctx, ctx.runtime.redundant_genome_ids)

        # enrichment stage
        hmm_search(ctx)
        blastp_search(ctx, ctx.paths.annotated_unique_prots_fasta)
        hmm_search(ctx)
        blastp_search(ctx, ctx.paths.predicted_unique_prots_fasta)
        hmm_search(ctx)

    # final stage
    merge_by_hmm(ctx)
    paralogs_hmm_search(ctx)
    rename_og_fastas_by_annotation(ctx)
    merge_by_annotation(ctx)

    ### POST-RUN ###
    rename_og_fastas_by_annotation(ctx)
    ctx.paths.output_dir.mkdir(exist_ok=True)
    ctx.ui.show("Writing reports...")
    utils.make_og_report(ctx, dataset_hash)
    ctx.paths.predicted_unique_prots_fasta.unlink()

    logger.info("Pipeline ended successfully")
