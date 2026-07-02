"""Endpoints, defaults and small environment helpers."""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

# --- Remote services -------------------------------------------------------
NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"

TOOL_NAME = "transcriptomics"
USER_AGENT = f"{TOOL_NAME}/0.1.0 (python-urllib)"

# --- HTTP behaviour --------------------------------------------------------
HTTP_TIMEOUT = 60      # seconds per request
HTTP_RETRIES = 3       # attempts before giving up


def ncbi_email() -> Optional[str]:
    """Optional contact email; NCBI E-utilities etiquette asks for one."""
    return os.environ.get("NCBI_EMAIL")


def ncbi_api_key() -> Optional[str]:
    """Optional NCBI API key — raises the rate limit from 3 to 10 req/s."""
    return os.environ.get("NCBI_API_KEY")


def detect_threads(reserve: int = 2, cap: Optional[int] = None) -> int:
    """Pick a sensible default thread count, leaving headroom for the OS."""
    cpu = os.cpu_count() or 1
    threads = max(1, cpu - reserve)
    if cap is not None:
        threads = min(threads, cap)
    return threads


@dataclass
class FetchConfig:
    """Resolved options for a single ``fetch`` invocation."""

    outdir: Path
    method: str = "auto"          # auto | ena | sra
    threads: int = 1
    max_spots: Optional[int] = None
    overwrite: bool = False
    dry_run: bool = False
    keep_intermediates: bool = False
