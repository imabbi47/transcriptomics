#!/usr/bin/env python3
"""Alignment-based quantification: HISAT2 (or STAR) -> featureCounts.

Actual splice-aware genome alignment + gene-level counting, as an alternative to
the alignment-free Salmon quantify stage. HISAT2 is the default because it is
memory-light (~4-8 GB for human); STAR is more accurate but needs ~30 GB RAM.

    python align.py --aligner hisat2 --index <hisat2_prefix> --gtf annotation.gtf \
        --r1 R1.fq.gz [--r2 R2.fq.gz] --outdir out [--strandedness 0|1|2]

Emits counts.txt (featureCounts) + a 2-column gene/count table ready for `de`.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess


def _which(explicit, *names):
    if explicit:
        return explicit
    for n in names:
        p = shutil.which(n)
        if p:
            return p
    return None


def _run(cmd):
    print("[align] running:", os.path.basename(cmd[0]), " ".join(cmd[1:3]), "...")
    # HISAT2/STAR ship as python/perl wrapper scripts; ensure the interpreters
    # sitting alongside the binary (e.g. a conda env's bin) are found on PATH.
    env = os.environ.copy()
    bindir = os.path.dirname(os.path.abspath(cmd[0]))
    env["PATH"] = bindir + os.pathsep + env.get("PATH", "")
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)
    if proc.returncode != 0:
        raise SystemExit(f"[align] {os.path.basename(cmd[0])} failed ({proc.returncode}):\n{proc.stdout[-1800:]}")
    return proc.stdout


def align_hisat2(hisat2, index, r1, r2, sam, threads):
    cmd = [hisat2, "-p", str(threads), "-x", index, "-S", sam, "--new-summary"]
    cmd += (["-1", r1, "-2", r2] if r2 else ["-U", r1])
    out = _run(cmd)
    tail = "\n".join(ln for ln in out.splitlines() if "%" in ln or "aligned" in ln.lower())
    if tail:
        print("[align] HISAT2 summary:\n" + tail)
    return sam


def align_star(star, genome_dir, r1, r2, outdir, threads):
    prefix = os.path.join(outdir, "star_")
    cmd = [star, "--runThreadN", str(threads), "--genomeDir", genome_dir, "--readFilesIn"]
    cmd += ([r1, r2] if r2 else [r1])
    if r1.endswith(".gz"):
        cmd += ["--readFilesCommand", "zcat"]
    cmd += ["--outSAMtype", "BAM", "Unsorted", "--outFileNamePrefix", prefix]
    _run(cmd)
    final = prefix + "Log.final.out"
    if os.path.exists(final):
        print("[align] STAR summary:")
        for line in open(final):
            if "%" in line or "Uniquely mapped" in line or "input reads" in line:
                print("   " + line.strip())
    return prefix + "Aligned.out.bam"


def count(featurecounts, gtf, aln, out, paired, strandedness, threads):
    cmd = [featurecounts, "-T", str(threads), "-a", gtf, "-o", out, "-s", str(strandedness)]
    if paired:
        cmd += ["-p", "--countReadPairs"]
    cmd += [aln]
    _run(cmd)


def summarize(counts_path):
    genes = nonzero = 0
    assigned = 0.0
    top = []
    with open(counts_path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("Geneid"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            gid, c = parts[0], float(parts[-1])
            genes += 1
            assigned += c
            if c > 0:
                nonzero += 1
            top.append((c, gid))
    print(f"\n[align] genes in annotation : {genes:,}")
    print(f"[align] genes with counts   : {nonzero:,}")
    print(f"[align] reads assigned      : {assigned:,.0f}")
    print("[align] top genes by count:")
    for c, g in sorted(top, reverse=True)[:10]:
        print(f"   {g:<26} {c:,.0f}")
    summ = counts_path + ".summary"
    if os.path.exists(summ):
        print("[align] featureCounts assignment:")
        for line in open(summ):
            k, *v = line.rstrip("\n").split("\t")
            if v and v[-1] not in ("", "0") and k != "Status":
                print(f"   {k}: {v[-1]}")
    # gene/count 2-column table for the de stage
    two_col = counts_path.replace(".txt", "_genecounts.tsv")
    with open(counts_path) as fh, open(two_col, "w") as out:
        out.write("gene_id\tcount\n")
        for line in fh:
            if line.startswith(("#", "Geneid")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                out.write(f"{parts[0]}\t{parts[-1]}\n")
    print(f"[align] wrote {two_col}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--aligner", choices=["hisat2", "star"], default="hisat2")
    ap.add_argument("--index", required=True, help="HISAT2 index prefix, or STAR genomeDir")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", default=None)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--strandedness", type=int, default=0, choices=[0, 1, 2],
                    help="featureCounts -s: 0 unstranded, 1 stranded, 2 reverse")
    ap.add_argument("--hisat2", default=None, help="path to hisat2 (else PATH)")
    ap.add_argument("--star", default=None, help="path to STAR (else PATH)")
    ap.add_argument("--featurecounts", default=None, help="path to featureCounts (else PATH)")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    fc = _which(a.featurecounts, "featureCounts")
    if not fc:
        raise SystemExit("[align] featureCounts not found — install:  conda install -c bioconda subread")

    if a.aligner == "hisat2":
        hisat2 = _which(a.hisat2, "hisat2")
        if not hisat2:
            raise SystemExit("[align] hisat2 not found — install:  conda install -c bioconda hisat2")
        aln = align_hisat2(hisat2, a.index, a.r1, a.r2, os.path.join(a.outdir, "aln.sam"), a.threads)
    else:
        star = _which(a.star, "STAR")
        if not star:
            raise SystemExit("[align] STAR not found — install:  conda install -c bioconda star")
        aln = align_star(star, a.index, a.r1, a.r2, a.outdir, a.threads)

    counts_out = os.path.join(a.outdir, "counts.txt")
    count(fc, a.gtf, aln, counts_out, bool(a.r2), a.strandedness, a.threads)
    summarize(counts_out)

    if a.aligner == "hisat2" and os.path.exists(aln):
        os.remove(aln)  # SAM is large; drop it after counting
    print(f"[align] done -> {counts_out}")


if __name__ == "__main__":
    main()
