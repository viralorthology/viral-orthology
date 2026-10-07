import logging
from pathlib import Path


def configure_logging(base_path: Path, dev_level: bool) -> None:
    """
    Configure application logging.

    Args:
        base_path: Directory where the log file will be created.
        dev_level: If True, enables DEBUG-level logging; otherwise,
            logging starts at INFO level.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG if dev_level else logging.INFO)

    if dev_level:
        format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    else:
        format = "%(asctime)s - %(message)s"
    formatter = logging.Formatter(
        format,
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    log_file = logging.FileHandler(base_path / "viralorthology.log")
    log_file.setFormatter(formatter)
    logger.addHandler(log_file)

    if dev_level:
        console = logging.StreamHandler()
        console.setLevel(logging.INFO)
        console.setFormatter(
            logging.Formatter(
                "%(asctime)s - %(message)s",
                datefmt="%H:%M:%S",
            )
        )
        logger.addHandler(console)
