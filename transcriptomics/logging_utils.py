"""Minimal, dependency-free logging setup.

Uses ``rich`` for colour if it happens to be installed, but falls back cleanly to
the standard library so the tool runs with zero third-party packages.
"""
from __future__ import annotations

import logging
import sys

_ROOT = "transcriptomics"


def setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    root = logging.getLogger(_ROOT)
    root.handlers.clear()
    root.setLevel(level)
    root.propagate = False

    try:  # pragma: no cover - cosmetic only
        from rich.logging import RichHandler

        handler: logging.Handler = RichHandler(
            show_time=True, show_path=False, rich_tracebacks=True, markup=False
        )
        handler.setFormatter(logging.Formatter("%(message)s", datefmt="[%H:%M:%S]"))
    except Exception:
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(
            logging.Formatter("%(asctime)s %(levelname)-7s %(message)s", datefmt="%H:%M:%S")
        )

    root.addHandler(handler)


def get_logger(name: str = _ROOT) -> logging.Logger:
    if name == _ROOT or name.startswith(_ROOT + "."):
        return logging.getLogger(name)
    return logging.getLogger(f"{_ROOT}.{name.rsplit('.', 1)[-1]}")
