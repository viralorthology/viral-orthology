import logging

import utils
from config.context import Context
from models.fasta import Fasta
from modules.pipeline.hmm_search import search_and_add_prots_with_hmm

logger = logging.getLogger(__name__)


def run(ctx: Context) -> None:
    # PRE-RUN
    utils.check_dependencies("hmmbuild", "muscle", "hmmsearch")
    utils.ensure_files_have_content(
        ctx.paths.genomes_fasta,
        ctx.paths.annotated_unique_prots_fasta,
    )
    utils.ensure_dir_has_files(ctx.paths.ortholog_groups_dir, ".fasta")

    genomes_fasta = Fasta(ctx.paths.genomes_fasta)

    # LOG
    dataset_hash = utils.get_dataset_hash(genomes_fasta.ids)
    utils.module_start_logging(
        logger, "HMM-based search module", dataset_hash, ctx.args.tool_args
    )

    # RUN
    genome_count = genomes_fasta.n_seqs
    min_n_proteins = genome_count - 1  # TODO accept arg from argv

    ogs = utils.get_ortholog_groups(
        ctx.paths.ortholog_groups_dir
    )  # TODO finish if no selected ogs at start
    annotated_unique_prots_fasta = Fasta(ctx.paths.annotated_unique_prots_fasta)

    while True:
        selected_ogs = [og for og in ogs if min_n_proteins <= og.n_seqs < genome_count]
        if not selected_ogs:
            break

        annotated_prots_added = search_and_add_prots_with_hmm(
            ctx.ui,
            selected_ogs,
            annotated_unique_prots_fasta,
            ctx.args.tool_args["hmm-search"],
        )

        if not annotated_prots_added:
            break

    # POST-RUN
    ctx.paths.output_dir.mkdir(exist_ok=True)
    ctx.ui.show("Writing reports...")
    utils.make_og_report(ctx, dataset_hash)

    logger.info("HMM-based search module completed succesfully")
