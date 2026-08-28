import shlex
import shutil
import subprocess


def run_cmd(cmd_str: str) -> str:
    """
    Run a command and return its standard output.

    Raises:
        subprocess.CalledProcessError: if the command exits with a non-zero status
        FileNotFoundError: if the command is not found
    """
    output = subprocess.run(
        shlex.split(cmd_str), capture_output=True, text=True, check=True
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
