import shlex
import shutil
import subprocess
from pathlib import Path


def run_cmd(cmd_str: str, cwd: Path | None = None) -> str:
    """
    Run a command and return its standard output.

    Args:
        cmd_str: The command to run.
        cwd: The working directory in which to run the command.

    Raises:
        subprocess.CalledProcessError: if the command exits with a non-zero status
        FileNotFoundError: if the command is not found
    """
    output = subprocess.run(
        shlex.split(cmd_str), capture_output=True, text=True, check=True, cwd=cwd
    )
    return output.stdout


def check_dependencies(*cmds: str) -> None:
    """
    Check that the given dependencies are available.

    Raises:
        FileNotFoundError: if the dependency is not found.
    """
    for cmd in cmds:
        if not shutil.which(cmd):
            raise FileNotFoundError(f"{cmd} not found.")
