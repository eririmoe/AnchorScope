from __future__ import annotations

import logging
from pathlib import Path

_our_handlers: list[logging.Handler] = []


def setup_logging(log_file: Path | None = None, level: str = "INFO", log_format: str | None = None) -> None:
    """Configure root logging for CLI execution."""
    global _our_handlers

    if log_format is None:
        log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

    root_logger = logging.getLogger()
    root_logger.setLevel(getattr(logging, level.upper()))

    for handler in _our_handlers[:]:
        root_logger.removeHandler(handler)
        handler.close()
        _our_handlers.remove(handler)

    formatter = logging.Formatter(log_format, datefmt="%Y-%m-%d %H:%M:%S")

    console_handler = logging.StreamHandler()
    console_handler.setLevel(getattr(logging, level.upper()))
    console_handler.setFormatter(formatter)
    root_logger.addHandler(console_handler)
    _our_handlers.append(console_handler)

    if log_file is not None:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file, encoding="utf-8")
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)
        _our_handlers.append(file_handler)
