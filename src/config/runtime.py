from dataclasses import dataclass


@dataclass
class Runtime:
    redundant_genomes: bool = False
    redundant_genome_ids: set[str] | None = None
    predicted_proteins_count: int = 0
    active_genomes: list[str] | None = None
