from dataclasses import dataclass

from cli.args import Args
from cli.ui import UI
from config.paths import Paths


@dataclass(frozen=True)
class Context:
    args: Args
    paths: Paths
    ui: UI
