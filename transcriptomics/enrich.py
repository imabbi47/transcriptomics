#!/usr/bin/env python3
"""Pathway over-representation on DE results via gseapy / Enrichr.

Takes the significant up- and down-regulated genes and tests them against
MSigDB Hallmark, KEGG and GO. Symbols are uppercased so mouse symbols query the
standard (human) Enrichr libraries as orthologs (harmless for human symbols).

    python enrich.py --results de_results.csv --outdir results --label "A vs B"
"""
from __future__ import annotations

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import gseapy as gp

LIBRARIES = ["MSigDB_Hallmark_2020", "KEGG_2021_Human", "GO_Biological_Process_2023"]


def enrich(symbols, tag, outdir, label):
    genes = sorted({str(s).upper() for s in symbols if isinstance(s, str) and s.strip()})
    genes = [g for g in genes if not g.startswith("ENSG")]  # drop unmapped Ensembl IDs
    if len(genes) < 10:
        print(f"[{tag}] only {len(genes)} usable genes — skipping")
        return None
    print(f"[{tag}] {len(genes)} genes -> Enrichr")
    try:
        enr = gp.enrichr(gene_list=genes, gene_sets=LIBRARIES, organism="human", outdir=None)
    except Exception as exc:  # noqa: BLE001
        print(f"[{tag}] enrichr failed: {exc!r}")
        return None
    d = enr.results.sort_values("Adjusted P-value")
    d.to_csv(os.path.join(outdir, f"enrichr_{tag}.csv"), index=False)
    print(f"\n=== {tag}-regulated ({label}) — top pathways ===")
    print(d.head(12)[["Gene_set", "Term", "Adjusted P-value", "Overlap"]].to_string(index=False))
    return d


def barplot(d, tag, outdir, label):
    if d is None or d.empty:
        return
    top = d.head(12).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.5, 5))
    ax.barh(range(len(top)),
            -np.log10(top["Adjusted P-value"].clip(lower=1e-300)),
            color="#c0392b" if tag == "up" else "#2c82c9")
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels([str(t)[:50] for t in top["Term"]], fontsize=8)
    ax.set_xlabel("-log10 adjusted p")
    ax.set_title(f"Enriched pathways — {tag}-regulated ({label})")
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, f"enrichment_{tag}.png"), dpi=130)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--label", default="contrast")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    res = pd.read_csv(args.results, index_col=0).dropna(subset=["padj"])
    sig = res[res["padj"] < args.alpha]
    up = sig[sig["log2FoldChange"] > 0]["symbol"]
    down = sig[sig["log2FoldChange"] < 0]["symbol"]

    for genes, tag in ((up, "up"), (down, "down")):
        barplot(enrich(genes, tag, args.outdir, args.label), tag, args.outdir, args.label)
    print(f"\n[enrich] wrote enrichr_*.csv + enrichment_*.png to {args.outdir}/")


if __name__ == "__main__":
    main()
