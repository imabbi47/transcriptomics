"""Resolve an accession into concrete sequencing runs.

Two complementary sources are used:

* **ENA portal API** (``filereport``) — accepts INSDC accessions (run/experiment/
  sample/study/bioproject for SRA, ENA and DDBJ) and, crucially, returns *direct
  FASTQ download URLs* in one call.
* **NCBI E-utilities** (``esearch`` + ``efetch`` runinfo) — the reliable way to
  turn a **GEO** accession (GSE/GSM) into its underlying SRA runs, and a fallback
  when ENA has no record.

Strategy: INSDC accessions go straight to ENA; GEO goes through NCBI and is then
enriched with ENA FASTQ URLs. Everything degrades gracefully — if no URLs are
found, the downloader falls back to the SRA Toolkit.
"""
from __future__ import annotations

import csv
import io
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Optional

from . import config
from .accession import Accession, AccessionType
from .logging_utils import get_logger
from .models import LibraryLayout, Run

log = get_logger(__name__)


class MetadataError(RuntimeError):
    """Raised when an accession cannot be resolved into runs."""


_ENA_FIELDS = [
    "run_accession",
    "experiment_accession",
    "sample_accession",
    "study_accession",
    "sample_title",
    "sample_alias",
    "scientific_name",
    "library_strategy",
    "library_layout",
    "instrument_platform",
    "instrument_model",
    "read_count",
    "base_count",
    "fastq_ftp",
    "fastq_md5",
    "fastq_bytes",
]


# --------------------------------------------------------------------------- #
# Low-level HTTP                                                               #
# --------------------------------------------------------------------------- #
def _http_get(url: str, timeout: Optional[int] = None, retries: Optional[int] = None) -> str:
    timeout = timeout or config.HTTP_TIMEOUT
    retries = retries or config.HTTP_RETRIES
    last_error: Optional[Exception] = None
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": config.USER_AGENT})
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read().decode("utf-8", errors="replace")
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as error:
            last_error = error
            log.debug("GET failed (%d/%d) %s: %s", attempt, retries, url, error)
            if attempt < retries:
                time.sleep(min(2 ** attempt, 8))
    raise MetadataError(f"HTTP request failed after {retries} attempts: {url} ({last_error})")


def _to_int(value: Optional[str]) -> Optional[int]:
    if value is None:
        return None
    value = value.strip()
    return int(value) if value.isdigit() else None


def _xml_tag(xml: str, tag: str) -> Optional[str]:
    match = re.search(rf"<{tag}>(.*?)</{tag}>", xml, re.DOTALL)
    return match.group(1).strip() if match else None


# --------------------------------------------------------------------------- #
# ENA                                                                         #
# --------------------------------------------------------------------------- #
def _ena_filereport(accession: str) -> List[Dict[str, str]]:
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": ",".join(_ENA_FIELDS),
        "format": "tsv",
        "limit": "0",  # 0 = no limit, return every run
    }
    url = f"{config.ENA_FILEREPORT}?{urllib.parse.urlencode(params)}"
    text = _http_get(url)
    if not text.strip():
        return []
    reader = csv.DictReader(io.StringIO(text), delimiter="\t")
    return [row for row in reader if row.get("run_accession")]


def _run_from_ena_row(row: Dict[str, str]) -> Run:
    urls = ["https://" + part for part in (row.get("fastq_ftp") or "").split(";") if part]
    md5 = [part for part in (row.get("fastq_md5") or "").split(";") if part]
    sizes = [int(part) for part in (row.get("fastq_bytes") or "").split(";") if part.isdigit()]
    layout = (row.get("library_layout") or "").upper()
    return Run(
        run_accession=row.get("run_accession", ""),
        experiment=row.get("experiment_accession") or None,
        sample_accession=row.get("sample_accession") or None,
        study_accession=row.get("study_accession") or None,
        sample_title=row.get("sample_title") or None,
        sample_alias=row.get("sample_alias") or None,
        organism=row.get("scientific_name") or None,
        library_strategy=row.get("library_strategy") or None,
        library_layout=LibraryLayout(layout) if layout in ("SINGLE", "PAIRED") else LibraryLayout.UNKNOWN,
        platform=row.get("instrument_platform") or None,
        instrument=row.get("instrument_model") or None,
        spots=_to_int(row.get("read_count")),
        bases=_to_int(row.get("base_count")),
        fastq_urls=urls,
        fastq_md5=md5,
        fastq_bytes=sizes,
    )


# --------------------------------------------------------------------------- #
# NCBI E-utilities                                                             #
# --------------------------------------------------------------------------- #
def _eutils_common() -> Dict[str, str]:
    params = {"tool": config.TOOL_NAME}
    email = config.ncbi_email()
    api_key = config.ncbi_api_key()
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    return params


def _ncbi_runinfo(term: str) -> List[Dict[str, str]]:
    """Search the SRA database for ``term`` and fetch the runinfo CSV.

    Works for GEO accessions (GSE/GSM) because NCBI cross-indexes them in SRA,
    as well as for any INSDC accession.
    """
    search = {**_eutils_common(), "db": "sra", "term": term, "usehistory": "y", "retmax": "0"}
    search_xml = _http_get(f"{config.NCBI_EUTILS}/esearch.fcgi?{urllib.parse.urlencode(search)}")

    count = _xml_tag(search_xml, "Count")
    web_env = _xml_tag(search_xml, "WebEnv")
    query_key = _xml_tag(search_xml, "QueryKey")
    if not count or count == "0" or not web_env or not query_key:
        return []

    fetch = {
        **_eutils_common(),
        "db": "sra",
        "query_key": query_key,
        "WebEnv": web_env,
        "rettype": "runinfo",
        "retmode": "text",
    }
    runinfo = _http_get(f"{config.NCBI_EUTILS}/efetch.fcgi?{urllib.parse.urlencode(fetch)}")
    reader = csv.DictReader(io.StringIO(runinfo))
    return [row for row in reader if row.get("Run")]


def _run_from_ncbi_row(row: Dict[str, str]) -> Run:
    layout = (row.get("LibraryLayout") or "").upper()
    return Run(
        run_accession=row.get("Run", ""),
        experiment=row.get("Experiment") or None,
        sample_accession=row.get("Sample") or None,
        study_accession=row.get("SRAStudy") or None,
        sample_title=row.get("SampleName") or None,
        sample_alias=row.get("SampleName") or None,
        organism=row.get("ScientificName") or None,
        library_strategy=row.get("LibraryStrategy") or None,
        library_layout=LibraryLayout(layout) if layout in ("SINGLE", "PAIRED") else LibraryLayout.UNKNOWN,
        platform=row.get("Platform") or None,
        instrument=row.get("Model") or None,
        spots=_to_int(row.get("spots")),
        bases=_to_int(row.get("bases")),
    )


def _geo_to_sra_runinfo(geo_accession: str) -> List[Dict[str, str]]:
    """GEO (GSE/GSM) → SRA runs via the ``gds`` database and ``elink``.

    GEO accessions are **not** searchable in the SRA Entrez database directly
    (a bare ``GSE…`` term returns "PhraseNotFound"), so we take the documented
    route instead:

    1. find the GEO record(s) in ``gds`` with ``<accession>[ACCN]``,
    2. ``elink`` that set from ``gds`` to ``sra`` on the history server,
    3. ``efetch`` the runinfo CSV for the linked runs.
    """
    search = {
        **_eutils_common(),
        "db": "gds",
        "term": f"{geo_accession}[ACCN]",
        "usehistory": "y",
        "retmax": "0",
    }
    search_xml = _http_get(f"{config.NCBI_EUTILS}/esearch.fcgi?{urllib.parse.urlencode(search)}")
    count = _xml_tag(search_xml, "Count")
    web_env = _xml_tag(search_xml, "WebEnv")
    query_key = _xml_tag(search_xml, "QueryKey")
    if not count or count == "0" or not web_env or not query_key:
        return []

    link = {
        **_eutils_common(),
        "dbfrom": "gds",
        "db": "sra",
        "query_key": query_key,
        "WebEnv": web_env,
        "cmd": "neighbor_history",
    }
    link_xml = _http_get(f"{config.NCBI_EUTILS}/elink.fcgi?{urllib.parse.urlencode(link)}")
    linked_key = _xml_tag(link_xml, "QueryKey")
    linked_env = _xml_tag(link_xml, "WebEnv") or web_env
    if not linked_key:
        return []

    fetch = {
        **_eutils_common(),
        "db": "sra",
        "query_key": linked_key,
        "WebEnv": linked_env,
        "rettype": "runinfo",
        "retmode": "text",
    }
    runinfo = _http_get(f"{config.NCBI_EUTILS}/efetch.fcgi?{urllib.parse.urlencode(fetch)}")
    reader = csv.DictReader(io.StringIO(runinfo))
    return [row for row in reader if row.get("Run")]


def _geo_sample_to_sra_accession(gsm: str) -> Optional[str]:
    """A GEO sample (GSM) maps to exactly one SRA experiment (SRX).

    Read the sample's ``esummary`` from ``gds`` and return the SRA accession from
    its ``ExtRelations``. This lets us resolve just *this* sample's runs, rather
    than the whole parent study that a ``gds → sra`` elink would over-select.
    """
    search = {**_eutils_common(), "db": "gds", "term": f"{gsm}[ACCN]", "retmax": "10"}
    search_xml = _http_get(f"{config.NCBI_EUTILS}/esearch.fcgi?{urllib.parse.urlencode(search)}")
    uids = re.findall(r"<Id>(\d+)</Id>", search_xml)
    # GEO sample UIDs in gds start with '3'; prefer that over series/platform.
    uid = next((u for u in uids if u.startswith("3")), uids[0] if uids else None)
    if not uid:
        return None

    summary = {**_eutils_common(), "db": "gds", "id": uid}
    summary_xml = _http_get(f"{config.NCBI_EUTILS}/esummary.fcgi?{urllib.parse.urlencode(summary)}")
    match = re.search(
        r'<Item Name="RelationType"[^>]*>SRA</Item>\s*'
        r'<Item Name="TargetObject"[^>]*>([^<]+)</Item>',
        summary_xml,
    )
    return match.group(1).strip() if match else None


def _enrich_with_ena_fastq(runs: List[Run]) -> None:
    """Fill in ENA FASTQ URLs for runs that don't have them yet (e.g. from GEO)."""
    pending = [run for run in runs if not run.fastq_urls and run.run_accession]
    if not pending:
        return

    def fill(run: Run) -> None:
        try:
            rows = _ena_filereport(run.run_accession)
        except MetadataError as error:
            log.debug("ENA enrichment failed for %s: %s", run.run_accession, error)
            return
        if not rows:
            return
        enriched = _run_from_ena_row(rows[0])
        run.fastq_urls = enriched.fastq_urls
        run.fastq_md5 = enriched.fastq_md5
        run.fastq_bytes = enriched.fastq_bytes
        if run.library_layout is LibraryLayout.UNKNOWN:
            run.library_layout = enriched.library_layout

    log.info("Looking up FASTQ URLs for %d run(s) via ENA…", len(pending))
    with ThreadPoolExecutor(max_workers=min(8, len(pending))) as pool:
        list(pool.map(fill, pending))


# --------------------------------------------------------------------------- #
# Public API                                                                  #
# --------------------------------------------------------------------------- #
def resolve_runs(accession: Accession) -> List[Run]:
    """Turn an :class:`Accession` into a sorted list of :class:`Run`."""
    if not accession.is_resolvable:
        raise MetadataError(
            f"Unrecognised accession: {accession.raw!r}. "
            "Expected something like SRR…, SRX…, SRP…, GSE… or GSM…"
        )

    if accession.type is AccessionType.GEO_SAMPLE:
        log.info("Resolving GEO sample %s via NCBI (gds → SRA experiment)…", accession.value)
        srx = _geo_sample_to_sra_accession(accession.value)
        if not srx:
            raise MetadataError(f"Could not find an SRA experiment linked to {accession.value}.")
        log.info("%s → %s", accession.value, srx)
        runs = [_run_from_ena_row(row) for row in _ena_filereport(srx)]
        if not runs:
            runs = [_run_from_ncbi_row(row) for row in _ncbi_runinfo(srx)]
            _enrich_with_ena_fastq(runs)
        for run in runs:
            run.geo_accession = accession.value
    elif accession.is_geo:
        log.info("Resolving GEO series %s via NCBI (gds → sra)…", accession.value)
        runs = [_run_from_ncbi_row(row) for row in _geo_to_sra_runinfo(accession.value)]
        _enrich_with_ena_fastq(runs)
    else:
        log.info("Resolving %s via ENA…", accession.value)
        runs = [_run_from_ena_row(row) for row in _ena_filereport(accession.value)]
        if not runs:
            log.info("ENA returned no records for %s; trying NCBI…", accession.value)
            runs = [_run_from_ncbi_row(row) for row in _ncbi_runinfo(accession.value)]
            _enrich_with_ena_fastq(runs)

    runs = [run for run in runs if run.run_accession]
    runs.sort(key=lambda run: run.run_accession)

    if not runs:
        raise MetadataError(
            f"No sequencing runs found for {accession.value}. It may be an "
            "array-only GEO record, embargoed, or not yet public."
        )
    log.info("Resolved %s → %d run(s).", accession.value, len(runs))
    return runs


def resolve_many(accessions: List[Accession]) -> List[Run]:
    """Resolve several accessions, de-duplicating runs by accession."""
    seen: Dict[str, Run] = {}
    for accession in accessions:
        for run in resolve_runs(accession):
            seen.setdefault(run.run_accession, run)
    return sorted(seen.values(), key=lambda run: run.run_accession)
