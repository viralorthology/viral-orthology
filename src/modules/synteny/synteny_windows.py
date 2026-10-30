import statistics
from collections import defaultdict
from typing import Literal

from models.fasta import Fasta
from models.ortholog_group import OrthologGroup

SWCategory = Literal["complete", "incomplete", "discarded"]


class SyntenyWindow:
    """
    Represent a pair of adjacent proteins and their strand orientations.

    Attributes:
        prots: Protein or ortholog-group identifiers in genomic order.
        strands: Strand orientations corresponding to the proteins. (direct strand = True)
    """

    def __init__(self, prot1: str, prot2: str, prot1_strand: bool, prot2_strand: bool):
        self.prots = [prot1, prot2]
        self.strands = [prot1_strand, prot2_strand]

    def __repr__(self) -> str:
        return f"SyntenyWindow({self.prots} {self.strands})"


def get_protein_pairs(
    proteomes_fasta: Fasta, non_paralog_ids: set[str]
) -> list[tuple[str, str, bool, bool]]:
    """
    Identify adjacent non-paralogous protein pairs in each genome.

    Groups proteins by genome, sorts them by their genomic midpoint, and
    generates pairs of consecutive proteins. Pairs containing proteins
    outside the set of non-paralogous identifiers are discarded.

    Returns:
        A list of tuples containing the two protein identifiers and their
        respective strand orientations.
    """
    seq_strands = {}
    seq_mid_positions = {}
    seqs_by_genome_id = defaultdict(list)
    for seq in proteomes_fasta.seqs:
        seqs_by_genome_id[seq.genome_id].append(seq)
        strand, mid_position = _get_seq_data(seq.description)
        seq_strands[seq.id] = strand
        seq_mid_positions[seq.id] = mid_position

    # get contiguos protein pairs
    protein_id_pairs = []
    for genome_proteins in seqs_by_genome_id.values():
        sorted_protein_ids = [
            seq.id
            for seq in sorted(
                genome_proteins, key=lambda seq: seq_mid_positions[seq.id]
            )
        ]

        genome_protein_id_pairs = _get_pairs_of_contiguos_proteins(sorted_protein_ids)
        protein_id_pairs.extend(genome_protein_id_pairs)

    # filter protein pairs with paralogs
    filtered_protein_pairs = []
    for pair in protein_id_pairs:
        if all(prot_id in non_paralog_ids for prot_id in pair):
            prot1_id, prot2_id = pair
            filtered_protein_pairs.append(
                (prot1_id, prot2_id, seq_strands[prot1_id], seq_strands[prot2_id])
            )

    return filtered_protein_pairs


def classify_sws(
    protein_pairs: list[tuple[str, str, bool, bool]], ogs: list[OrthologGroup]
) -> tuple[list[SyntenyWindow], list[SyntenyWindow], list[SyntenyWindow]]:
    """
    Classify protein pairs according to their ortholog-group membership.

    Replaces protein identifiers belonging to ortholog groups with their
    corresponding group names and classifies each pair according to how
    many proteins belong to an ortholog group.

    Returns:
        A tuple containing three lists of synteny windows:
            - complete
            - incomplete
            - discarded
    """
    og_path_name_by_prot_id = {}
    seqs_in_og_ids = set()
    for og in ogs:
        for seq in og.seqs:
            og_path_name_by_prot_id[seq.id] = og.path.name
            seqs_in_og_ids.add(seq.id)

    # split synteny windows into categories
    synteny_windows: dict[SWCategory, list[SyntenyWindow]] = {
        "complete": [],
        "incomplete": [],
        "discarded": [],
    }

    for protein_pair in protein_pairs:
        prot1_id, prot2_id, prot1_strand, prot2_strand = protein_pair
        prot1 = og_path_name_by_prot_id.get(prot1_id, prot1_id)
        prot2 = og_path_name_by_prot_id.get(prot2_id, prot2_id)

        sw = SyntenyWindow(prot1, prot2, prot1_strand, prot2_strand)

        sw_category = _classify_sw(prot1_id, prot2_id, seqs_in_og_ids)
        synteny_windows[sw_category].append(sw)

    return (
        synteny_windows["complete"],
        synteny_windows["incomplete"],
        synteny_windows["discarded"],
    )


def _classify_sw(prot1_id: str, prot2_id: str, seqs_in_og_ids: set[str]) -> SWCategory:
    """
    Determine the category of a protein pair from ortholog-group membership.

    Returns:
        "complete" if both proteins belong to ortholog groups, "incomplete"
        if exactly one does, or "discarded" if neither does.
    """
    assert seqs_in_og_ids
    og_protein_count = sum(id_ in seqs_in_og_ids for id_ in (prot1_id, prot2_id))

    match og_protein_count:
        case 2:
            return "complete"
        case 1:
            return "incomplete"
        case 0:
            return "discarded"
        case _:
            raise ValueError(f"Unexpected count: {og_protein_count}")


def _get_pairs_of_contiguos_proteins(
    sorted_protein_ids: list[str],
) -> list[tuple[str, str]]:
    """Generate pairs of consecutive protein identifiers."""
    pairs = []
    first, last = 0, 2
    while True:
        gene_pair = sorted_protein_ids[first:last]
        if len(gene_pair) < 2:
            break

        pairs.append((gene_pair[0], gene_pair[1]))

        first += 1
        last += 1

    return pairs


def _get_seq_data(seq_description: str) -> tuple[bool, float]:
    assert "[location=" in seq_description

    seq_location_str = seq_description.split("[location=")[1].split("]")[0]

    return _get_seq_strand(seq_location_str), _get_seq_mid_position(seq_location_str)


def _get_seq_strand(seq_location_str: str) -> bool:
    """
    Determine the strand orientation from a sequence description.

    Returns:
        True if the sequence is on the forward strand, or False if it is
        on the complementary strand.
    """
    return "complement" not in seq_location_str


def _get_seq_mid_position(seq_location_str: str) -> float:
    """
    Calculate the midpoint of a sequence's genomic location interval.

    Returns:
        The arithmetic mean of the interval's start and end coordinates.
    """
    clean_location = ("").join([c for c in seq_location_str if c in "0123456789."])
    split_location = clean_location.split(".")
    location = [int(split_location[0]), int(split_location[-1])]

    return statistics.mean(location)
