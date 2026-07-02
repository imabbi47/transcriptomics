"""Shared data structures for the ingest stage.

``Run`` is the canonical description of a single sequencing run, merged from ENA
and/or NCBI metadata. ``RunFiles`` pairs a run with the FASTQ files we actually
wrote to disk.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import List, Optional


class LibraryLayout(str, Enum):
    SINGLE = "SINGLE"
    PAIRED = "PAIRED"
    UNKNOWN = "UNKNOWN"


def _safe(name: str) -> str:
    """Make a string safe to use as a file/sample name."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip()).strip("._-")
    return cleaned


@dataclass
class Run:
    """One sequencing run plus enough metadata to download and label it."""

    run_accession: str
    experiment: Optional[str] = None
    sample_accession: Optional[str] = None
    study_accession: Optional[str] = None
    sample_title: Optional[str] = None
    sample_alias: Optional[str] = None
    organism: Optional[str] = None
    library_strategy: Optional[str] = None
    library_layout: LibraryLayout = LibraryLayout.UNKNOWN
    platform: Optional[str] = None
    instrument: Optional[str] = None
    spots: Optional[int] = None
    bases: Optional[int] = None
    # Direct-download info from ENA (semicolon lists, split out). May be empty,
    # in which case the downloader falls back to the SRA Toolkit.
    fastq_urls: List[str] = field(default_factory=list)
    fastq_md5: List[str] = field(default_factory=list)
    fastq_bytes: List[int] = field(default_factory=list)
    # Set when the run was reached via a GEO accession.
    geo_accession: Optional[str] = None

    @property
    def has_ena_fastq(self) -> bool:
        return any(self.fastq_urls)

    @property
    def is_paired(self) -> bool:
        if self.library_layout is LibraryLayout.PAIRED:
            return True
        if self.library_layout is LibraryLayout.SINGLE:
            return False
        # Unknown layout: infer from the number of ENA FASTQ files.
        paired_files = [u for u in self.fastq_urls if u.endswith(("_1.fastq.gz", "_2.fastq.gz"))]
        return len(paired_files) >= 2

    @property
    def total_bytes(self) -> int:
        return sum(self.fastq_bytes)

    def friendly_name(self) -> str:
        """A human-meaningful sample label (GEO/alias/title) for the samplesheet.

        Files on disk are named by ``run_accession`` for uniqueness; this is the
        biological label used in the ``sample`` column.
        """
        for candidate in (self.geo_accession, self.sample_alias, self.sample_title):
            if candidate and candidate.strip():
                safe = _safe(candidate)
                if safe:
                    return safe
        return self.run_accession


@dataclass
class RunFiles:
    """A run together with the FASTQ files written for it."""

    run: Run
    fastq_1: Path
    fastq_2: Optional[Path] = None

    @property
    def is_paired(self) -> bool:
        return self.fastq_2 is not None
