"""Write the artifacts downstream steps consume.

* ``samplesheet.csv`` — the input contract for QC / trimming / quantification.
  Columns mirror common RNA-seq conventions: ``sample, run_accession, fastq_1,
  fastq_2, strandedness, organism, layout``.
* ``runs.json`` — full resolved metadata for provenance and debugging.
"""
from __future__ import annotations

import csv
import json
from dataclasses import asdict
from pathlib import Path
from typing import List, Optional

from .models import Run, RunFiles

SAMPLESHEET_COLUMNS = [
    "sample",
    "run_accession",
    "fastq_1",
    "fastq_2",
    "strandedness",
    "organism",
    "layout",
]


def _row(run: Run, fastq_1: str, fastq_2: str) -> List[str]:
    return [
        run.friendly_name(),
        run.run_accession,
        fastq_1,
        fastq_2,
        "auto",  # let the quantifier auto-detect strandedness later
        run.organism or "",
        run.library_layout.value,
    ]


def write_samplesheet(run_files: List[RunFiles], path: Path) -> Path:
    """Write a samplesheet describing FASTQ files that exist on disk."""
    path = Path(path)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(SAMPLESHEET_COLUMNS)
        for item in run_files:
            writer.writerow(
                _row(
                    item.run,
                    str(item.fastq_1.resolve()),
                    str(item.fastq_2.resolve()) if item.fastq_2 else "",
                )
            )
    return path


def write_planned_samplesheet(runs: List[Run], path: Path) -> Path:
    """Dry-run variant: predict file names without requiring them to exist."""
    path = Path(path)
    outdir = path.parent.resolve()
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(SAMPLESHEET_COLUMNS)
        for run in runs:
            acc = run.run_accession
            fastq_1 = str(outdir / f"{acc}_R1.fastq.gz")
            fastq_2 = str(outdir / f"{acc}_R2.fastq.gz") if run.is_paired else ""
            writer.writerow(_row(run, fastq_1, fastq_2))
    return path


def write_runs_json(runs: List[Run], path: Path) -> Path:
    """Dump full run metadata as JSON (enums rendered as their values)."""
    path = Path(path)

    def _clean(run: Run) -> dict:
        data = asdict(run)
        data["library_layout"] = run.library_layout.value
        return data

    path.write_text(json.dumps([_clean(run) for run in runs], indent=2))
    return path
