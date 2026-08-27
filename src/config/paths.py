from dataclasses import dataclass, field
from pathlib import Path


@dataclass(frozen=True)
class Paths:
    """Filesystem paths used throughout the pipeline."""

    base: Path
    output_dir: Path = field(init=False)
    proteomes_dir: Path = field(init=False)
    genomes_dir: Path = field(init=False)
    paralogs_dir: Path = field(init=False)
    ortholog_groups_dir: Path = field(init=False)

    annotated_prots_db: Path = field(init=False)
    predicted_prots_db: Path = field(init=False)
    genomes_fasta: Path = field(init=False)
    proteomes_fasta: Path = field(init=False)
    orfeomes_fasta: Path = field(init=False)
    ids_txt: Path = field(init=False)

    def __post_init__(self) -> None:
        if not isinstance(self.base, Path):
            raise TypeError(
                f"base must be a pathlib.Path, got {type(self.base).__name__}"
            )
        base = self.base.resolve()

        object.__setattr__(self, "base", base)
        object.__setattr__(self, "output_dir", base / "output")
        object.__setattr__(self, "proteomes_dir", base / "proteomes")
        object.__setattr__(self, "genomes_dir", base / "genomes")
        object.__setattr__(self, "paralogs_dir", base / "paralogs")
        object.__setattr__(self, "ortholog_groups_dir", base / "ortholog_groups")
        object.__setattr__(self, "annotated_prots_db", base / "gb_unique_prots.fasta")
        object.__setattr__(
            self, "predicted_prots_db", base / "predicted_unique_prots.fasta"
        )
        object.__setattr__(self, "genomes_fasta", base / "genomes.fasta")
        object.__setattr__(self, "proteomes_fasta", base / "proteomes.fasta")
        object.__setattr__(self, "orfeomes_fasta", base / "orfeomes.fasta")
        object.__setattr__(self, "ids_txt", base / "ids.txt")
