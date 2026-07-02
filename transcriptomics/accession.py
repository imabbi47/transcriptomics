"""Accession parsing and classification.

Recognises the accession families used by NCBI SRA, EMBL-EBI ENA, DDBJ and GEO so
the rest of the pipeline can decide *how* to resolve an identifier into sequencing
runs:

    GEO series  (GSE…)  -> many samples -> many runs
    GEO sample  (GSM…)  -> one sample   -> one or more runs
    study       (SRP/ERP/DRP, PRJ…)     -> many runs
    experiment  (SRX/ERX/DRX)           -> one or more runs
    sample      (SRS/ERS/DRS, SAM…)     -> one or more runs
    run         (SRR/ERR/DRR)           -> the unit we actually download

Only classification lives here (pure, dependency-free, easy to unit test). The
network calls that turn these into concrete runs live in ``metadata.py``.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from enum import Enum


class AccessionType(str, Enum):
    """The family an accession belongs to."""

    GEO_SERIES = "geo_series"   # GSE…
    GEO_SAMPLE = "geo_sample"   # GSM…
    RUN = "run"                 # SRR / ERR / DRR
    EXPERIMENT = "experiment"   # SRX / ERX / DRX
    SAMPLE = "sample"           # SRS / ERS / DRS
    STUDY = "study"             # SRP / ERP / DRP
    BIOPROJECT = "bioproject"   # PRJNA / PRJEB / PRJDB
    BIOSAMPLE = "biosample"     # SAMN / SAMEA / SAMD
    UNKNOWN = "unknown"


# Ordered most-specific first; the first match wins.
_PATTERNS: list[tuple[AccessionType, "re.Pattern[str]"]] = [
    (AccessionType.GEO_SERIES, re.compile(r"^GSE\d+$", re.IGNORECASE)),
    (AccessionType.GEO_SAMPLE, re.compile(r"^GSM\d+$", re.IGNORECASE)),
    (AccessionType.RUN, re.compile(r"^[SED]RR\d+$", re.IGNORECASE)),
    (AccessionType.EXPERIMENT, re.compile(r"^[SED]RX\d+$", re.IGNORECASE)),
    (AccessionType.SAMPLE, re.compile(r"^[SED]RS\d+$", re.IGNORECASE)),
    (AccessionType.STUDY, re.compile(r"^[SED]RP\d+$", re.IGNORECASE)),
    (AccessionType.BIOPROJECT, re.compile(r"^PRJ(NA|EB|DB)\d+$", re.IGNORECASE)),
    (AccessionType.BIOSAMPLE, re.compile(r"^SAM[NED][A-Z]?\d+$", re.IGNORECASE)),
]


@dataclass(frozen=True)
class Accession:
    """A user-supplied identifier together with its detected family."""

    raw: str
    type: AccessionType

    @property
    def value(self) -> str:
        """Normalised (upper-cased) accession; APIs are case-insensitive but tidy."""
        return self.raw.upper()

    @property
    def is_geo(self) -> bool:
        return self.type in (AccessionType.GEO_SERIES, AccessionType.GEO_SAMPLE)

    @property
    def is_run(self) -> bool:
        return self.type is AccessionType.RUN

    @property
    def is_resolvable(self) -> bool:
        """False only for identifiers we don't recognise at all."""
        return self.type is not AccessionType.UNKNOWN

    def __str__(self) -> str:
        return f"{self.value} ({self.type.value})"


def classify(accession: str) -> Accession:
    """Classify a single accession string into an :class:`Accession`.

    Unrecognised input is returned as ``AccessionType.UNKNOWN`` rather than
    raising, so callers can collect and report every bad token at once.
    """
    raw = accession.strip()
    for acc_type, pattern in _PATTERNS:
        if pattern.match(raw):
            return Accession(raw=raw, type=acc_type)
    return Accession(raw=raw, type=AccessionType.UNKNOWN)


def parse_many(accessions: list[str]) -> list[Accession]:
    """Classify many accessions, skipping blank tokens."""
    return [classify(a) for a in accessions if a and a.strip()]
