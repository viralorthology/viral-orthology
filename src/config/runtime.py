from dataclasses import dataclass


@dataclass
class Runtime:
    redundant_genomes: bool = False
    predicted_proteins_count: int = 0
    active_genomes: list[str] | None = None
