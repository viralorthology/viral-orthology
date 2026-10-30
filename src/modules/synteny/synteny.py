import logging
from collections import Counter
from pathlib import Path

import engines
import utils
from config.context import Context
from models.fasta import Fasta
from models.ortholog_group import OrthologGroup
from models.seq import Seq
from models.tmp_fasta import TmpFasta
from modules.synteny.synteny_windows import (
    SyntenyWindow,
    classify_sws,
    get_protein_pairs,
)

SYNTENY_MIN_CONSERVATION_PERC = 0.9  # TODO pass user perc

logger = logging.getLogger(__name__)


def run(ctx: Context) -> None:
    # PRE-RUN
    utils.check_dependencies("blastp", "makeblastdb")
    utils.ensure_files_have_content(
        ctx.paths.genomes_fasta,
        ctx.paths.proteomes_fasta,
        ctx.paths.annotated_unique_prots_fasta,
    )
    utils.ensure_dir_has_files(ctx.paths.ortholog_groups_dir, ".fasta")

    genome_ids = Fasta(ctx.paths.genomes_fasta).ids

    # LOG
    dataset_hash = utils.get_dataset_hash(genome_ids)
    utils.module_start_logging(
        logger, "Synteny module", dataset_hash, ctx.args.tool_args
    )

    # RUN
    proteomes_fasta = Fasta(ctx.paths.proteomes_fasta)
    annotated_unique_proteins_fasta = Fasta(ctx.paths.annotated_unique_prots_fasta)
    ogs = utils.get_ortholog_groups(ctx.paths.ortholog_groups_dir)
    non_paralog_ids = set(annotated_unique_proteins_fasta.ids).union(
        {seq.id for og in ogs for seq in og.seqs}
    )

    protein_pairs = get_protein_pairs(proteomes_fasta, non_paralog_ids)

    complete_sws, incomplete_sws, discarded_sws = classify_sws(protein_pairs, ogs)

    conserved_sws = _get_conserved_sws(complete_sws, len(genome_ids))
    if not conserved_sws or not incomplete_sws:
        logger.info("No conserved or incomplete synteny window found")
        return

    logger.info("Conserved synteny windows found:")
    for sw in conserved_sws:
        strands = ["direct" if strand else "reverse" for strand in sw.strands]
        logger.info("%s %s", (" ").join(sw.prots), (" ").join(strands))

    added_prot_ids: set[str] = set()
    for incomplete_sw in incomplete_sws:
        if set(incomplete_sw.prots).intersection(added_prot_ids):
            continue

        compatible_conserved_sws = [
            sw for sw in conserved_sws if _conserved_sw_is_compatible(incomplete_sw, sw)
        ]
        if not compatible_conserved_sws:
            continue

        assert len(compatible_conserved_sws) == 1
        conserved_sw = next(iter(compatible_conserved_sws))

        og, potential_new_ortholog = _get_og_and_potential_new_protein(
            set(conserved_sw.prots),
            set(incomplete_sw.prots),
            annotated_unique_proteins_fasta,
            ctx.paths.ortholog_groups_dir,
        )
        if potential_new_ortholog.genome_id in og.genome_ids:
            continue

        if _has_blastp_hit(
            og, potential_new_ortholog, ctx.args.tool_args["blastp-search"]
        ):
            og.add_seqs(potential_new_ortholog)
            added_prot_ids.add(potential_new_ortholog.id)
            logger.info(
                "%s added to %s based on synteny conservation",
                potential_new_ortholog.id,
                og.path.name,
            )

            # update discarded sws
            for sw in discarded_sws.copy():
                if potential_new_ortholog.id in sw.prots:
                    sw.prots[sw.prots.index(potential_new_ortholog.id)] = og.path.name
                    incomplete_sws.append(sw)

                    assert sw in discarded_sws
                    discarded_sws.remove(sw)

    if added_prot_ids:
        annotated_unique_proteins_fasta.remove_seqs(*added_prot_ids)

    # POST-RUN
    ctx.paths.output_dir.mkdir(exist_ok=True)
    ctx.ui.show("Writing reports...")
    utils.make_og_report(ctx, dataset_hash)

    logger.info("Synteny module completed succesfully")


def _get_conserved_sws(
    complete_sws: list[SyntenyWindow], n_genomes: int
) -> list[SyntenyWindow]:
    """
    Identify synteny windows conserved across genomes.

    Groups compatible synteny windows and selects representative
    windows whose compatibility counts meet the minimum conservation
    threshold without reaching the total number of genomes.

    Args:
        complete_sws: Complete synteny windows to compare.
        n_genomes: Total number of genomes in the dataset.

    Returns:
        Representative synteny windows meeting the conservation criteria.
    """
    grouped_sws = set()
    compatibility_counter: Counter[SyntenyWindow] = Counter()
    for i, sw1 in enumerate(complete_sws):
        if sw1 in grouped_sws:
            continue
        compatibility_counter[sw1] += 1

        for sw2 in complete_sws[i + 1 :]:
            if sw2 in grouped_sws:
                continue

            if _sws_are_compatible(sw1, sw2):
                compatibility_counter[sw1] += 1
                grouped_sws.add(sw2)

    convserved_sws = []
    for sw, n in compatibility_counter.most_common():
        if n < n_genomes and n >= n_genomes * SYNTENY_MIN_CONSERVATION_PERC:
            convserved_sws.append(sw)

    return convserved_sws


def _sws_are_compatible(
    sw1: SyntenyWindow,
    sw2: SyntenyWindow,
) -> bool:
    """
    Check whether two synteny windows have equivalent gene order.

    Windows are compatible if their protein identifiers and strand
    orientations match in the same order or in reverse order, with
    strand orientations inverted in the reverse configuration.
    """
    same_sw = sw1.prots == sw2.prots and sw1.strands == sw2.strands
    reverse_sw = sw1.prots == sw2.prots[::-1] and sw1.strands == [
        not strand for strand in sw2.strands[::-1]
    ]

    return same_sw or reverse_sw


def _conserved_sw_is_compatible(
    incomplete_sw: SyntenyWindow,
    conserved_sw: SyntenyWindow,
) -> bool:
    """
    Check whether an incomplete window matches a conserved window.

    Determines whether the windows share an ortholog-group identifier
    and whether their strand orientations are consistent with the
    relative position of that shared group.
    """
    og_in_common = set(incomplete_sw.prots).intersection(set(conserved_sw.prots))

    if not og_in_common:
        return False

    assert len(og_in_common) == 1

    og_path_name = next(iter(og_in_common))

    same_order = conserved_sw.prots.index(og_path_name) == incomplete_sw.prots.index(
        og_path_name
    )

    if same_order:
        return conserved_sw.strands == incomplete_sw.strands
    return conserved_sw.strands == [
        not strand for strand in incomplete_sw.strands[::-1]
    ]


def _get_og_and_potential_new_protein(
    conserved_sw_prots: set[str],
    incomplete_sw_prots: set[str],
    unique_prots_fasta: Fasta,
    ogs_dir_path: Path,
) -> tuple[OrthologGroup, Seq]:
    """
    Retrieve the ortholog group and candidate protein from two windows.

    Identifies the ortholog group represented only in the conserved
    window and the candidate protein represented only in the incomplete
    window. Loads the ortholog group and retrieves the candidate sequence
    from the unique-proteins FASTA file.
    """
    og_path_name = next(iter(conserved_sw_prots.difference(incomplete_sw_prots)))
    potential_new_prot_id = next(
        iter(incomplete_sw_prots.difference(conserved_sw_prots))
    )

    potential_new_prot = unique_prots_fasta.get_seqs(potential_new_prot_id)[0]
    og = OrthologGroup(ogs_dir_path / og_path_name)

    return og, potential_new_prot


def _has_blastp_hit(
    og: OrthologGroup, potential_new_protein: Seq, blastp_params: str
) -> bool:
    """Check whether a candidate protein has a BLASTP hit against an ortholog group."""
    with TmpFasta() as blast_db:
        blast_db.add_seqs(*og.seqs)
        engines.make_blast_db(blast_db, "prot")

        with TmpFasta() as query_fasta:
            query_fasta.add_seqs(potential_new_protein)

            blast_hits = engines.blastp_search(query_fasta, blast_db, blastp_params)

    return bool(blast_hits)
