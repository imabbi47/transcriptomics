"""transcriptomics — accession-driven RNA-seq preprocessing.

Stage 0 (ingest): resolve an SRA / ENA / DDBJ / GEO accession into sequencing runs
and download them as analysis-ready FASTQ, emitting a samplesheet that downstream
steps (QC, trimming, quantification) consume directly.
"""
from __future__ import annotations

__version__ = "0.1.0"

from .accession import Accession, AccessionType, classify, parse_many

__all__ = [
    "__version__",
    "Accession",
    "AccessionType",
    "classify",
    "parse_many",
]
