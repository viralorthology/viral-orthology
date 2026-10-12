from dataclasses import dataclass

import utils
from engines.muscle import muscle_align
from models.fasta import Fasta
from models.ortholog_group import OrthologGroup


@dataclass(frozen=True)
class HMMHit:
    evalue: float
    seq_id: str
    genome_id: str


def hmm_search(
    query_fasta: OrthologGroup, db_fasta: Fasta, params: str = ""
) -> list[HMMHit]:
    """
    Search a protein database using an HMM built from an ortholog group.

    Args:
        query_fasta: Ortholog group used to build the query HMM.
        db_fasta: FASTA file containing the protein sequences to search.
        params: Additional command-line parameters to pass to hmmsearch.

    Returns:
        A list of HMMER hits sorted by E-value.
    """
    if not query_fasta.alignment_path.is_file():
        muscle_align(query_fasta)
        _build_hmm_hmmer(query_fasta)

    hits = _search_hmm_hmmer(query_fasta, db_fasta, params)

    return hits


def _build_hmm_hmmer(fasta: OrthologGroup) -> None:
    """
    Build an HMM from a multiple sequence alignment using HMMER.
    """
    utils.run_cmd(f"hmmbuild {fasta.hmm_hmmer_path} {fasta.alignment_path}")


def _search_hmm_hmmer(
    query_fasta: OrthologGroup, db_fasta: Fasta, search_params: str
) -> list[HMMHit]:
    """
    Search a protein database using an HMM with hmmsearch.

    Returns:
        A list of HMMER hits sorted by E-value.
    """
    output = utils.run_cmd(
        f"hmmsearch {search_params} --noali {query_fasta.hmm_hmmer_path} {db_fasta.path}"
    )
    return _parse_hmmsearch_output(output)


def _parse_hmmsearch_output(hmmsearch_output: str) -> list[HMMHit]:
    """
    Parse HMMER hmmsearch output into HMMHit objects.

    Args:
        hmmsearch_output: Raw output produced by hmmsearch.

    Returns:
        A list of HMMER hits sorted by E-value.
    """

    def _parse_hmmsearch_hit(line: str) -> HMMHit:
        fields = line.split()

        return HMMHit(
            evalue=float(fields[0]),
            seq_id=fields[8],
            genome_id=fields[9],
        )

    hits = []
    read = False

    for line in hmmsearch_output.splitlines():
        line = line.strip()
        if line.startswith("E-value  score  bias"):
            read = True
            continue

        if read:
            if (
                line.startswith(
                    (
                        "------ inclusion threshold ------",
                        "Domain annotation for each sequence",
                    )
                )
                or not line
            ):
                break

            if line.startswith("-"):
                continue

            hit = _parse_hmmsearch_hit(line)
            hits.append(hit)

    return sorted(
        hits, key=lambda hit: (hit.evalue, hit.seq_id)
    )  # Ensure reproducibility
