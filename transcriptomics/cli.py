"""Command-line interface for the ingest stage.

    transcriptomics fetch   SRR000001 GSE12345 [-o data] [--dry-run] …
    transcriptomics resolve SRP123456            # metadata only, no download

Stdlib-only (argparse), so it runs without installing anything beyond the package.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Sequence

from . import __version__, config
from .accession import Accession, parse_many
from .download import DownloadError, download_run
from .logging_utils import get_logger, setup_logging
from .metadata import MetadataError, resolve_runs
from .models import Run, RunFiles
from .samplesheet import write_planned_samplesheet, write_runs_json, write_samplesheet

log = get_logger("cli")


# --------------------------------------------------------------------------- #
# Presentation helpers                                                         #
# --------------------------------------------------------------------------- #
def _human_bytes(num: int) -> str:
    if not num:
        return "?"
    value = float(num)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024 or unit == "TB":
            return f"{value:.0f}{unit}" if unit == "B" else f"{value:.1f}{unit}"
        value /= 1024
    return f"{value:.1f}TB"


def _print_runs_table(runs: Sequence[Run]) -> None:
    header = ("RUN", "SAMPLE", "LAYOUT", "ORGANISM", "SIZE", "SOURCE")
    rows = [
        (
            run.run_accession,
            run.friendly_name()[:22],
            run.library_layout.value,
            (run.organism or "")[:20],
            _human_bytes(run.total_bytes),
            "ENA" if run.has_ena_fastq else "SRA",
        )
        for run in runs
    ]
    widths = [max(len(str(col)), *(len(str(r[i])) for r in rows)) if rows else len(col)
              for i, col in enumerate(header)]
    line = "  ".join(str(col).ljust(widths[i]) for i, col in enumerate(header))
    print(line)
    print("  ".join("-" * widths[i] for i in range(len(header))))
    for row in rows:
        print("  ".join(str(col).ljust(widths[i]) for i, col in enumerate(row)))
    total = sum(run.total_bytes for run in runs)
    if total:
        print(f"\nTotal download size (ENA estimate): ~{_human_bytes(total)} across {len(runs)} run(s).")


def _resolve_all(accessions: List[Accession]) -> List[Run]:
    runs: dict = {}
    for accession in accessions:
        for run in resolve_runs(accession):
            runs.setdefault(run.run_accession, run)
    return sorted(runs.values(), key=lambda run: run.run_accession)


# --------------------------------------------------------------------------- #
# Commands                                                                     #
# --------------------------------------------------------------------------- #
def cmd_resolve(args: argparse.Namespace) -> int:
    accessions = parse_many(args.accessions)
    unknown = [a.raw for a in accessions if not a.is_resolvable]
    if unknown:
        log.error("Unrecognised accession(s): %s", ", ".join(unknown))
        return 2
    runs = _resolve_all(accessions)
    _print_runs_table(runs)
    return 0


def cmd_fetch(args: argparse.Namespace) -> int:
    accessions = parse_many(args.accessions)
    unknown = [a.raw for a in accessions if not a.is_resolvable]
    if unknown:
        log.error("Unrecognised accession(s): %s", ", ".join(unknown))
        return 2

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    log.info("Resolving %d accession(s)…", len(accessions))
    runs = _resolve_all(accessions)
    _print_runs_table(runs)

    write_runs_json(runs, outdir / "runs.json")

    if args.dry_run:
        sheet = write_planned_samplesheet(runs, outdir / "samplesheet.csv")
        log.info("Dry run — wrote metadata only:")
        log.info("  %s", (outdir / "runs.json"))
        log.info("  %s  (planned)", sheet)
        return 0

    threads = args.threads or config.detect_threads()
    log.info("Downloading %d run(s) to %s using %d thread(s)…", len(runs), outdir, threads)

    downloaded: List[RunFiles] = []
    failures: List[str] = []
    for index, run in enumerate(runs, start=1):
        log.info("[%d/%d] %s", index, len(runs), run.run_accession)
        try:
            downloaded.append(
                download_run(
                    run,
                    outdir,
                    method=args.method,
                    threads=threads,
                    max_spots=args.max_spots,
                    overwrite=args.overwrite,
                    keep_intermediates=args.keep_intermediates,
                )
            )
        except DownloadError as error:
            log.error("  failed: %s", error)
            failures.append(run.run_accession)

    if downloaded:
        sheet = write_samplesheet(downloaded, outdir / "samplesheet.csv")
        log.info("Wrote samplesheet: %s", sheet)

    log.info("Done: %d/%d run(s) downloaded.", len(downloaded), len(runs))
    if failures:
        log.warning("Failed: %s", ", ".join(failures))
        return 1
    return 0


# --------------------------------------------------------------------------- #
# Parser                                                                       #
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="transcriptomics",
        description="Fetch SRA/ENA/DDBJ/GEO sequencing runs as analysis-ready FASTQ.",
        epilog=(
            "pipeline stages (some need extras: pip install -e '.[analysis]'):\n"
            "  quantify     Salmon quantification of FASTQ (needs the salmon binary)\n"
            "  de           differential expression (pyDESeq2)\n"
            "  enrich       GO / KEGG / Hallmark enrichment (gseapy)\n"
            "  report       bundle DE + enrichment into one HTML file\n"
            "  geo-design   parse a GEO series matrix into a design table\n"
            "run 'transcriptomics <stage> --help' for stage options."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("-v", "--verbose", action="store_true", help="debug logging")
    sub = parser.add_subparsers(dest="command", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("accessions", nargs="+", metavar="ACCESSION",
                        help="one or more of SRR/SRX/SRP/SRS, ERR…, DRR…, PRJ…, GSE…, GSM…")

    p_resolve = sub.add_parser("resolve", parents=[common],
                               help="resolve accessions to runs and print a table (no download)")
    p_resolve.set_defaults(func=cmd_resolve)

    p_fetch = sub.add_parser("fetch", parents=[common],
                             help="resolve and download runs, then write a samplesheet")
    p_fetch.add_argument("-o", "--outdir", default="data", help="output directory (default: ./data)")
    p_fetch.add_argument("--method", choices=("auto", "ena", "sra"), default="auto",
                         help="download route (default: auto — ENA if available, else SRA Toolkit)")
    p_fetch.add_argument("--threads", type=int, default=0, help="threads for fasterq-dump (default: auto)")
    p_fetch.add_argument("--max-spots", type=int, default=None,
                         help="download only the first N reads (quick test runs; SRA route)")
    p_fetch.add_argument("--overwrite", action="store_true", help="re-download even if files exist")
    p_fetch.add_argument("--keep-intermediates", action="store_true",
                         help="keep .sra/work files instead of cleaning up")
    p_fetch.add_argument("--dry-run", action="store_true",
                         help="resolve and write metadata only; download nothing")
    p_fetch.set_defaults(func=cmd_fetch)
    return parser


_STAGES = {"quantify": "quantify", "de": "de", "enrich": "enrich",
           "report": "report", "geo-design": "geo_design"}


def _run_stage(stage: str, rest: Sequence[str]) -> int:
    """Delegate to a pipeline-stage module (lazy import keeps the core dep-free)."""
    import importlib

    saved = sys.argv
    sys.argv = [f"transcriptomics {stage}", *rest]
    try:
        module = importlib.import_module(f"transcriptomics.{_STAGES[stage]}")
        return module.main() or 0
    except ModuleNotFoundError as exc:
        log.error("stage '%s' needs extra packages (missing: %s). "
                  "Install with:  pip install -e '.[analysis]'", stage, exc.name)
        return 1
    finally:
        sys.argv = saved


def main(argv: Optional[Sequence[str]] = None) -> int:
    raw = list(sys.argv[1:] if argv is None else argv)
    if raw and raw[0] in _STAGES:
        return _run_stage(raw[0], raw[1:])
    args = build_parser().parse_args(argv)
    setup_logging(getattr(args, "verbose", False))
    try:
        return args.func(args)
    except MetadataError as error:
        log.error("%s", error)
        return 1
    except KeyboardInterrupt:  # pragma: no cover
        log.error("Interrupted.")
        return 130


if __name__ == "__main__":
    sys.exit(main())
