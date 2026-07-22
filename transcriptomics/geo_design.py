#!/usr/bin/env python3
"""Parse a GEO series_matrix.txt.gz into a design table.

Summarizes each sample characteristic (so you can pick the comparison) and
writes a design CSV keyed by sample title. Stdlib only.

    python geo_design.py --matrix GSE..._series_matrix.txt.gz --out design_full.csv
"""
from __future__ import annotations

import argparse
import collections
import csv
import gzip
import re


def read_lines(path):
    opener = gzip.open(path, "rt", errors="replace") if path.endswith(".gz") else open(path)
    with opener as fh:
        for line in fh:
            yield line.rstrip("\n")


def split_row(line):
    return [c.strip().strip('"') for c in line.split("\t")[1:]]


def parse(path):
    titles, gsms, sources = [], [], []
    chars = collections.OrderedDict()  # "idx:key" -> list of values (per sample)
    ci = 0
    for line in read_lines(path):
        if line.startswith("!Sample_title"):
            titles = split_row(line)
        elif line.startswith("!Sample_geo_accession"):
            gsms = split_row(line)
        elif line.startswith("!Sample_source_name_ch1"):
            sources = split_row(line)
        elif line.startswith("!Sample_characteristics_ch1"):
            vals = split_row(line)
            keys_seen = {v.split(":", 1)[0].strip().lower() for v in vals if ":" in v}
            key = vals[0].split(":")[0].strip() if vals and ":" in vals[0] else f"char{ci}"
            key = re.sub(r"\s+", "_", key.lower())
            if len(keys_seen) > 1:
                print(f"[geo] WARNING: characteristics row {ci:02d} mixes keys across samples "
                      f"{sorted(keys_seen)} -> column '{ci:02d}_{key}' may be mislabelled; "
                      "verify the design before running DE")
            chars[f"{ci:02d}_{key}"] = [v.split(":", 1)[1].strip() if ":" in v else v for v in vals]
            ci += 1
    return titles, gsms, sources, chars


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", required=True)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    titles, gsms, sources, chars = parse(args.matrix)
    n = len(titles)
    print(f"[geo] samples: {n}")
    print(f"[geo] title examples: {', '.join(titles[:6])}")
    if sources:
        print(f"[geo] source examples: {', '.join(sorted(set(sources))[:4])}")

    print("\n[geo] characteristics (value distribution):")
    for key, vals in chars.items():
        counts = collections.Counter(vals)
        # show compactly; useful comparison variables have few levels
        summary = ", ".join(f"{v}={c}" for v, c in counts.most_common(8))
        more = "" if len(counts) <= 8 else f"  (+{len(counts)-8} more)"
        print(f"   {key} [{len(counts)} levels]: {summary}{more}")

    if args.out:
        keys = list(chars.keys())
        with open(args.out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["sample", "gsm", "source"] + keys)
            for i in range(n):
                row = [titles[i],
                       gsms[i] if i < len(gsms) else "",
                       sources[i] if i < len(sources) else ""]
                row += [chars[k][i] if i < len(chars[k]) else "" for k in keys]
                w.writerow(row)
        print(f"\n[geo] wrote design table -> {args.out}")


if __name__ == "__main__":
    main()
