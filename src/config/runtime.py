from dataclasses import dataclass, field


@dataclass
class Runtime:
    redundant_genomes: bool = False
    redundant_genome_ids: set[str] = field(default_factory=set)
    predicted_proteins_count: int = 0
