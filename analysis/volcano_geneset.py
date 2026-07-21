#!/usr/bin/env python3
"""Volcano plot with a gene set highlighted over the full DE background."""
from __future__ import annotations

import argparse

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--de", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--title", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--lfc", type=float, default=1.0)
    ap.add_argument("--nlabel", type=int, default=12)
    ap.add_argument("--ymax", type=float, default=0.0, help="y-axis cap; 0 = auto from gene set")
    args = ap.parse_args()

    gset = {ln.strip().upper() for ln in open(args.genes) if ln.strip()}
    de = pd.read_csv(args.de).dropna(subset=["padj", "log2FoldChange"]).copy()
    de["nlp"] = -np.log10(de["padj"].clip(lower=1e-300))
    de["inset"] = de["symbol"].astype(str).str.upper().isin(gset)

    bg = de[~de["inset"]]
    st = de[de["inset"]]
    sig = st[st["padj"] < args.alpha]
    nsig = st[st["padj"] >= args.alpha]
    up, dn = sig[sig["log2FoldChange"] > 0], sig[sig["log2FoldChange"] < 0]

    fig, ax = plt.subplots(figsize=(7.2, 6))
    ax.scatter(bg["log2FoldChange"], bg["nlp"], s=4, c="#d9dce1", alpha=.5, linewidths=0, label="all genes")
    ax.scatter(nsig["log2FoldChange"], nsig["nlp"], s=16, c="#f0a24b", alpha=.8, linewidths=0, label="gene set (ns)")
    ax.scatter(up["log2FoldChange"], up["nlp"], s=24, c="#c0392b", edgecolors="white", linewidths=.3, label="gene set up (sig)")
    ax.scatter(dn["log2FoldChange"], dn["nlp"], s=24, c="#2c6fb0", edgecolors="white", linewidths=.3, label="gene set down (sig)")

    ax.axhline(-np.log10(args.alpha), ls="--", c="grey", lw=.7)
    ax.axvline(args.lfc, ls="--", c="grey", lw=.7)
    ax.axvline(-args.lfc, ls="--", c="grey", lw=.7)

    for _, r in sig.sort_values("padj").head(args.nlabel).iterrows():
        ax.annotate(str(r["symbol"]), (r["log2FoldChange"], r["nlp"]),
                    fontsize=7, xytext=(3, 3), textcoords="offset points")

    xmax = min(10.0, max(3.0, float(np.ceil(np.percentile(np.abs(de["log2FoldChange"]), 99)))))
    ax.set_xlim(-xmax, xmax)
    # zoom y to the gene set so it isn't dwarfed by a few ultra-significant background genes
    set_max = float(st["nlp"].max()) if len(st) else float(de["nlp"].max())
    ymax = args.ymax if args.ymax > 0 else max(5.0, float(np.ceil(set_max * 1.15)))
    ax.set_ylim(-ymax * 0.03, ymax)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title(args.title, fontsize=10)
    ax.legend(fontsize=7, loc="upper right", framealpha=.9)
    ax.text(.02, .98, f"gene set: {len(st)} found · {len(sig)} significant",
            transform=ax.transAxes, va="top", fontsize=8, color="#333")
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"[volcano] {args.out}  (set={len(st)}, sig up={len(up)}, down={len(dn)})")


if __name__ == "__main__":
    main()
