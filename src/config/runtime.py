from dataclasses import dataclass, field


@dataclass
class Runtime:
    redundant_genome_ids: set[str] = field(default_factory=set)
    active_genome_ids: set[str] = field(default_factory=set)
    predicted_proteins_count: int = 0
    paralog_ids: set[str] = field(default_factory=set)
