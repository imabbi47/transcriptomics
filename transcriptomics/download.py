"""Download a run's reads as analysis-ready, gzipped FASTQ.

Two routes, chosen automatically:

1. **ENA direct FASTQ** (preferred) — files are already gzipped FASTQ, so we just
   fetch them with ``curl``/``wget`` and verify the MD5. Fast, no conversion.
2. **SRA Toolkit fallback** — ``prefetch`` + ``fasterq-dump`` (or ``fastq-dump``
   when ``--max-spots`` is set for quick test runs), then gzip with ``pigz``/``gzip``.

Output files are named ``<run_accession>_R1.fastq.gz`` (+ ``_R2`` when paired) so
downstream steps have a predictable layout.
"""
from __future__ import annotations

import hashlib
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional, Tuple

from .logging_utils import get_logger
from .models import Run, RunFiles

log = get_logger(__name__)


class DownloadError(RuntimeError):
    """Raised when a run cannot be downloaded."""


# --------------------------------------------------------------------------- #
# Small process / file helpers                                                #
# --------------------------------------------------------------------------- #
def _which(*names: str) -> Optional[str]:
    for name in names:
        found = shutil.which(name)
        if found:
            return found
    return None


def _run_cmd(cmd: List[str]) -> str:
    log.debug("exec: %s", " ".join(cmd))
    proc = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    )
    if proc.returncode != 0:
        tail = (proc.stdout or "")[-2000:]
        raise DownloadError(f"command failed ({proc.returncode}): {' '.join(cmd)}\n{tail}")
    return proc.stdout or ""


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _gzip(path: Path) -> Path:
    """Gzip ``path`` in place (preferring pigz), returning the .gz path."""
    target = path.with_suffix(path.suffix + ".gz")
    pigz = _which("pigz")
    if pigz:
        _run_cmd([pigz, "-f", str(path)])
    elif _which("gzip"):
        _run_cmd(["gzip", "-f", str(path)])
    else:  # pure-python fallback
        import gzip as _gz

        with path.open("rb") as src, _gz.open(target, "wb") as dst:
            shutil.copyfileobj(src, dst)
        path.unlink(missing_ok=True)
    return target


# --------------------------------------------------------------------------- #
# ENA direct download                                                         #
# --------------------------------------------------------------------------- #
def _pick_fastq_pair(urls: List[str]) -> Tuple[Optional[str], Optional[str]]:
    """From ENA's URL list, choose the R1 (and R2) we actually want.

    ENA sometimes also lists a bare ``<run>.fastq.gz`` alongside ``_1``/``_2``;
    we prefer the explicitly paired files when present.
    """
    read1 = next((u for u in urls if u.endswith("_1.fastq.gz")), None)
    read2 = next((u for u in urls if u.endswith("_2.fastq.gz")), None)
    if read1 and read2:
        return read1, read2
    singles = [u for u in urls if not u.endswith(("_1.fastq.gz", "_2.fastq.gz"))]
    if singles:
        return singles[0], None
    if read1:
        return read1, None
    return (urls[0] if urls else None), None


def _http_download(url: str, dest: Path, expected_md5: Optional[str], overwrite: bool) -> Path:
    if dest.exists() and not overwrite:
        if expected_md5 and _md5(dest) != expected_md5:
            log.warning("%s exists but MD5 differs — re-downloading.", dest.name)
        else:
            log.info("✓ %s already present.", dest.name)
            return dest

    tmp = dest.with_name(dest.name + ".part")
    tmp.unlink(missing_ok=True)
    curl = _which("curl")
    wget = _which("wget")
    if curl:
        # --speed-limit/--speed-time: abort if the transfer stalls below
        # 10 KB/s for 60 s, so --retry can recover instead of hanging forever.
        _run_cmd([curl, "-L", "--fail", "--silent", "--show-error",
                  "--retry", "5", "--retry-delay", "5",
                  "--speed-limit", "10240", "--speed-time", "60",
                  "-o", str(tmp), url])
    elif wget:
        _run_cmd([wget, "--quiet", "--tries=3", "-O", str(tmp), url])
    else:
        import urllib.request

        log.debug("curl/wget not found; using urllib for %s", url)
        urllib.request.urlretrieve(url, tmp)  # noqa: S310 - URL is from ENA

    if expected_md5:
        actual = _md5(tmp)
        if actual != expected_md5:
            tmp.unlink(missing_ok=True)
            raise DownloadError(f"MD5 mismatch for {dest.name}: {actual} != {expected_md5}")

    tmp.replace(dest)
    log.info("✓ downloaded %s", dest.name)
    return dest


def _download_ena(run: Run, outdir: Path, overwrite: bool) -> RunFiles:
    url1, url2 = _pick_fastq_pair(run.fastq_urls)
    if not url1:
        raise DownloadError(f"{run.run_accession}: no usable ENA FASTQ URL")

    md5_by_url = dict(zip(run.fastq_urls, run.fastq_md5)) if run.fastq_md5 else {}
    acc = run.run_accession

    dest1 = _http_download(url1, outdir / f"{acc}_R1.fastq.gz", md5_by_url.get(url1), overwrite)
    dest2: Optional[Path] = None
    if url2:
        dest2 = _http_download(url2, outdir / f"{acc}_R2.fastq.gz", md5_by_url.get(url2), overwrite)
    return RunFiles(run=run, fastq_1=dest1, fastq_2=dest2)


# --------------------------------------------------------------------------- #
# SRA Toolkit fallback                                                        #
# --------------------------------------------------------------------------- #
def _download_sra(
    run: Run, outdir: Path, threads: int, max_spots: Optional[int], overwrite: bool, keep: bool
) -> RunFiles:
    acc = run.run_accession
    final1 = outdir / f"{acc}_R1.fastq.gz"
    final2 = outdir / f"{acc}_R2.fastq.gz"
    if final1.exists() and not overwrite:
        log.info("✓ %s already present.", final1.name)
        return RunFiles(run, final1, final2 if final2.exists() else None)

    if not _which("fasterq-dump", "fastq-dump"):
        raise DownloadError(
            "SRA Toolkit not found and no ENA FASTQ available for "
            f"{acc}. Install it with:  conda install -c bioconda sra-tools"
        )

    workdir = outdir / f".{acc}.work"
    workdir.mkdir(parents=True, exist_ok=True)
    try:
        prefetch = _which("prefetch")
        if prefetch:
            log.info("prefetch %s …", acc)
            _run_cmd([prefetch, acc, "-O", str(workdir), "--max-size", "100g"])

        if max_spots and _which("fastq-dump"):
            log.info("fastq-dump %s (first %d spots) …", acc, max_spots)
            _run_cmd([_which("fastq-dump"), "-X", str(max_spots), "--split-files",
                      "--gzip", "-O", str(workdir), acc])
            produced = sorted(workdir.glob(f"{acc}*.fastq.gz"))
        else:
            log.info("fasterq-dump %s (%d threads) …", acc, threads)
            _run_cmd([_which("fasterq-dump"), acc, "-O", str(workdir),
                      "--split-files", "-e", str(threads), "-f"])
            produced = [_gzip(p) for p in sorted(workdir.glob(f"{acc}*.fastq"))]

        produced = [p for p in produced if not p.name.endswith("_3.fastq.gz")]
        if not produced:
            raise DownloadError(f"{acc}: SRA Toolkit produced no FASTQ files")

        read1 = next((p for p in produced if p.name.endswith("_1.fastq.gz")), None)
        read2 = next((p for p in produced if p.name.endswith("_2.fastq.gz")), None)
        if read1 and read2:
            read1.replace(final1)
            read2.replace(final2)
            return RunFiles(run, final1, final2)
        single = read1 or produced[0]
        single.replace(final1)
        return RunFiles(run, final1, None)
    finally:
        if not keep:
            shutil.rmtree(workdir, ignore_errors=True)


# --------------------------------------------------------------------------- #
# Public API                                                                  #
# --------------------------------------------------------------------------- #
def download_run(
    run: Run,
    outdir: Path,
    method: str = "auto",
    threads: int = 1,
    max_spots: Optional[int] = None,
    overwrite: bool = False,
    keep_intermediates: bool = False,
) -> RunFiles:
    """Download a single run, returning the FASTQ files written."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    use_ena = method == "ena" or (method == "auto" and run.has_ena_fastq)
    if method == "ena" and not run.has_ena_fastq:
        raise DownloadError(f"{run.run_accession}: ENA has no FASTQ; try --method sra")

    if use_ena:
        log.info("Downloading %s from ENA…", run.run_accession)
        return _download_ena(run, outdir, overwrite)

    log.info("Downloading %s via SRA Toolkit…", run.run_accession)
    return _download_sra(run, outdir, threads, max_spots, overwrite, keep_intermediates)
