#!/usr/bin/env python3
"""Quantify a sample against a Salmon index and summarize the result.

Prototype of the pipeline's quantify stage: runs `salmon quant` (auto library
type, selective alignment) then reports mapping rate + top transcripts. Stdlib.

    python quantify.py --salmon <bin> --index <dir> --r1 R1.fq.gz [--r2 R2.fq.gz] --outdir out
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess


def gencode_symbol(name):
    # GENCODE FASTA headers: ENST|ENSG|OTTHUMG|OTTHUMT|tx_name|gene_name|len|...
    parts = name.split("|")
    return parts[5] if len(parts) > 5 else name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--salmon", required=True)
    ap.add_argument("--index", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", default=None)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    cmd = [args.salmon, "quant", "-i", args.index, "-l", "A",
           "-p", str(args.threads), "--validateMappings", "-o", args.outdir]
    cmd += (["-1", args.r1, "-2", args.r2] if args.r2 else ["-r", args.r1])
    print(f"[quant] {' '.join(os.path.basename(c) for c in cmd[:3])} ...")
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(proc.stdout[-1200:])
    if proc.returncode != 0:
        raise SystemExit(f"salmon quant failed ({proc.returncode})")

    meta_path = os.path.join(args.outdir, "aux_info", "meta_info.json")
    meta = json.load(open(meta_path)) if os.path.exists(meta_path) else {}
    print("\n[quant] === mapping summary ===")
    print(f"[quant] reads processed : {meta.get('num_processed', 0):,}")
    print(f"[quant] reads mapped    : {meta.get('num_mapped', 0):,}")
    if "percent_mapped" in meta:
        print(f"[quant] mapping rate    : {meta['percent_mapped']:.2f}%")
    print(f"[quant] library type    : {meta.get('library_types', ['?'])}")

    rows = []
    with open(os.path.join(args.outdir, "quant.sf")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rows.append((row["Name"], float(row["TPM"]), float(row["NumReads"])))
    nonzero = [r for r in rows if r[1] > 0]
    print(f"\n[quant] transcripts in index : {len(rows):,}")
    print(f"[quant] transcripts TPM>0    : {len(nonzero):,}")
    print("\n[quant] top 12 transcripts by TPM:")
    for name, tpm, nreads in sorted(rows, key=lambda x: x[1], reverse=True)[:12]:
        print(f"   {gencode_symbol(name):<14} TPM={tpm:9.1f}  reads={nreads:8.0f}")


if __name__ == "__main__":
    main()
