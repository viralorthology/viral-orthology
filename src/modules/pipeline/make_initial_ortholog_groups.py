import shutil
import tempfile
from pathlib import Path

import engines
import utils
from config.context import Context
from models.fasta import Fasta
from models.seq import Seq


def make_initial_ortholog_groups(ctx: Context) -> None:
    """
    Create initial ortholog groups from the non-redundant proteomes.

    Proteomes are compared with ProteinOrtho to generate ortholog groups.
    Groups containing paralogs are filtered by retaining the paralog with
    the best BLASTP hit against the other proteins in the group. Groups
    containing proteins from only one genome are deleted.
    """
    ctx.ui.show("Creating initial ortholog groups...")

    proteome_fastas = utils.get_fastas(ctx.paths.proteomes_dir, ".proteome")
    non_redundant_proteome_fastas = [
        fasta
        for fasta in proteome_fastas
        if fasta.path.stem  # file stem == genome id
        not in ctx.runtime.redundant_genome_ids
    ]

    with tempfile.TemporaryDirectory() as tmp_path:
        # copy proteomes to temp dir
        tmp_dir = Path(tmp_path)

        tmp_proteome_fastas = []
        for fasta in non_redundant_proteome_fastas:
            shutil.copy(fasta.path, tmp_dir / fasta.path.name)
            tmp_proteome_fastas.append(Fasta(tmp_dir / fasta.path.name))

        # run proteinortho
        _run_proteinortho(
            tmp_proteome_fastas, ctx.args.tool_args["proteinortho"], tmp_dir
        )

        # remove paralogs from ortholog groups
        ortholog_groups = utils.get_fastas(tmp_dir, ".fasta")
        assert ortholog_groups
        deleted_ogs = set()
        for og in ortholog_groups:
            og_was_deleted = _remove_paralogs(og, ctx.paths.paralogs_dir)
            if og_was_deleted:
                deleted_ogs.add(og)

        # move ortholog groups
        ctx.paths.ortholog_groups_dir.mkdir()
        for og in ortholog_groups:
            if og in deleted_ogs:
                continue
            og.move_fasta(ctx.paths.ortholog_groups_dir)


def _run_proteinortho(proteome_fastas: list[Fasta], params: str, cwd: Path) -> None:
    """
    Run ProteinOrtho and extract the resulting ortholog group proteins.

    Args:
        proteome_fastas: Proteome FASTA files to compare.
        params: Additional parameters to pass to ProteinOrtho.
        cwd: Working directory for ProteinOrtho and its output files.
    """
    fasta_paths = (" ").join(str(fasta.path) for fasta in proteome_fastas)

    utils.run_cmd(
        f"proteinortho {fasta_paths} {params}",
        cwd=cwd,
    )
    utils.run_cmd(
        f"proteinortho_grab_proteins.pl -tofiles myproject.proteinortho.tsv -exact {fasta_paths}",
        cwd=cwd,
    )


def _remove_paralogs(og: Fasta, paralogs_dir_path: Path) -> bool:
    """
    Remove paralogs from an ortholog group and save them in the paralogs dir.

    For each genome containing multiple proteins in the group, the paralog
    with the best BLASTP hit against proteins from the other genomes is
    retained in the group. The remaining paralogs are removed from the
    group and saved to the paralogs directory.
    An ortholog group containing proteins from only one genome is deleted.

    Args:
        og: Ortholog group FASTA file to process.
        paralogs_dir_path: Directory where removed paralog sequences are saved.

    Returns:
        True if the ortholog group was deleted because it contained proteins
        from only one genome, otherwise False.
    """
    if len(set(og.genome_ids)) == 1:  # all seqs from the same genome
        og.delete_fasta()
        return True

    for genome_id in _get_genome_ids_with_paralogs(og.genome_ids):
        paralog_seqs, other_og_seqs = _filter_seqs_for_evaluation(og, genome_id)
        best_hit_seq, other_paralog_seqs = _evaluate_paralogs(
            paralog_seqs, other_og_seqs
        )

        # remove paralogs from og
        og.remove_seqs(*[seq.id for seq in other_paralog_seqs])
        assert og.path.exists()

        # move paralogs to paralogs dir
        paralogs_fasta = Fasta(paralogs_dir_path / f"{best_hit_seq.id}.fasta")
        paralogs_fasta.add_seqs(*other_paralog_seqs)

    return False


def _get_genome_ids_with_paralogs(genome_ids: list[str]) -> list[str]:
    """
    Find genome IDs that occur more than once.

    Args:
        genome_ids: List of genome IDs associated with the sequences in an
            orthologous group.

    Returns:
        Sorted list of genome IDs that have multiple sequences.
    """
    assert len(set(genome_ids)) > 1

    return sorted(  # Ensure reproducibility
        [genome_id for genome_id in set(genome_ids) if genome_ids.count(genome_id) > 1]
    )


def _filter_seqs_for_evaluation(
    og: Fasta, genome_id: str
) -> tuple[list[Seq], list[Seq]]:
    """
    Separate paralogs from the remaining proteins in an ortholog group.

    Args:
        og: Ortholog group containing the sequences to separate.
        genome_id: Genome ID whose paralog sequences should be used as
            BLASTP queries.

    Returns:
        - The paralogous sequences from genome_id
        - The sequences from all other genomes in the ortholog group.
    """
    paralog_seqs = []
    other_og_seqs = []

    for seq in og.seqs:
        if seq.genome_id == genome_id:
            paralog_seqs.append(seq)
        else:
            other_og_seqs.append(seq)

    assert len(paralog_seqs) > 1
    assert other_og_seqs

    return paralog_seqs, other_og_seqs


def _evaluate_paralogs(
    paralog_seqs: list[Seq], other_og_seqs: list[Seq]
) -> tuple[Seq, list[Seq]]:
    """
    Identify the best-supported paralog using BLASTP.

    The paralog sequences are used as BLASTP queries against a database
    containing the proteins from other genomes in the ortholog group.
    The query producing the best hit is retained, while the remaining
    paralogs are returned for removal.

    Args:
        paralog_seqs: Paralogous sequences from the same genome.
        other_og_seqs: Sequences from other genomes in the ortholog
            group, used to build the BLASTP database.

    Returns:
        - The best-supported paralog
        - The remaining paralog sequences.
    """
    with tempfile.TemporaryDirectory() as tmp_path:
        tmp_dir = Path(tmp_path)

        assert other_og_seqs
        fasta_blastp_db = Fasta(tmp_dir / "blastp_db.fasta")
        fasta_blastp_db.add_seqs(*other_og_seqs)
        engines.make_blast_db(fasta_blastp_db, "prot")

        assert len(paralog_seqs) > 1
        fasta_query = Fasta(tmp_dir / "blastp_query.fasta")
        fasta_query.add_seqs(*paralog_seqs)

        hits = engines.blastp_search(fasta_query, fasta_blastp_db, "-evalue 100000")
        assert hits
        best_query_id = hits[0].query_id

        best_seq = None
        other_seqs = []
        for seq in paralog_seqs:
            if seq.id == best_query_id:
                assert best_seq is None
                best_seq = seq
            else:
                other_seqs.append(seq)
        assert best_seq is not None
        assert other_seqs

        return best_seq, other_seqs
