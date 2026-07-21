#!/usr/bin/env python3
"""Bundle analysis results into a multi-sheet .xlsx for plotting in R / Excel.

Sheets: Info, DE_results, PCA, Enrichment_up, Enrichment_down.

    python export_excel.py --results DIR --counts counts.tsv.gz --design design.csv --out out.xlsx
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd


def load_counts(path):
    sep = "\t" if path.endswith((".tsv.gz", ".txt.gz", ".tsv", ".txt")) else ","
    df = pd.read_csv(path, sep=sep, index_col=0)
    return df.apply(pd.to_numeric, errors="coerce").dropna(how="all")


def compute_pca(counts, design, sample_col, n_top=500):
    """PCA of samples on log2-CPM of the most variable genes."""
    from sklearn.decomposition import PCA

    samples = [s for s in design[sample_col].astype(str) if s in counts.columns]
    X = counts[samples]
    cpm = X / X.sum(axis=0) * 1e6
    logcpm = np.log2(cpm + 1)
    top = logcpm.loc[logcpm.var(axis=1).sort_values(ascending=False).index[:n_top]]
    mat = top.T.values
    pca = PCA(n_components=2).fit(mat - mat.mean(axis=0))
    coords = pca.transform(mat - mat.mean(axis=0))
    pct = pca.explained_variance_ratio_ * 100
    out = pd.DataFrame({sample_col: samples, "PC1": coords[:, 0], "PC2": coords[:, 1]})
    out = out.merge(design.astype({sample_col: str}), on=sample_col, how="left")
    return out, pct


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True, help="dir with de_results.csv + enrichr_*.csv")
    ap.add_argument("--counts", required=True)
    ap.add_argument("--design", required=True)
    ap.add_argument("--sample-col", default="sample")
    ap.add_argument("--top", type=int, default=50, help="max enrichment rows per direction")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    de = pd.read_csv(os.path.join(args.results, "de_results.csv"))
    design = pd.read_csv(args.design)

    def enr(name):
        p = os.path.join(args.results, name)
        return pd.read_csv(p).head(args.top) if os.path.exists(p) else pd.DataFrame()

    up_df, dn_df = enr("enrichr_up.csv"), enr("enrichr_down.csv")

    pca_df, pct = pd.DataFrame(), (float("nan"), float("nan"))
    try:
        pca_df, pct = compute_pca(load_counts(args.counts), design, args.sample_col)
    except Exception as exc:  # never let PCA block the export
        print(f"[xlsx] PCA skipped ({exc})")

    info = pd.DataFrame({
        "sheet": ["DE_results", "PCA", "Enrichment_up", "Enrichment_down"],
        "contents": [
            "Per-gene differential expression (log2FoldChange, padj, baseMean, symbol)",
            f"Sample PCA coordinates (PC1={pct[0]:.1f}% var, PC2={pct[1]:.1f}% var) + design",
            "Enriched pathways among UP-regulated genes",
            "Enriched pathways among DOWN-regulated genes",
        ],
        "make_this_graph_in_R": [
            "Volcano: log2FoldChange vs -log10(padj) | MA: log2(baseMean) vs log2FoldChange",
            "PCA scatter: PC1 vs PC2, colour by the design factor",
            "Bar chart: -log10(Adjusted P-value) per Term",
            "Bar chart: -log10(Adjusted P-value) per Term",
        ],
    })

    with pd.ExcelWriter(args.out, engine="openpyxl") as xl:
        info.to_excel(xl, sheet_name="Info", index=False)
        de.to_excel(xl, sheet_name="DE_results", index=False)
        if not pca_df.empty:
            pca_df.to_excel(xl, sheet_name="PCA", index=False)
        if not up_df.empty:
            up_df.to_excel(xl, sheet_name="Enrichment_up", index=False)
        if not dn_df.empty:
            dn_df.to_excel(xl, sheet_name="Enrichment_down", index=False)

    print(f"[xlsx] wrote {args.out} ({os.path.getsize(args.out)/1024:.0f} KB) — "
          f"{de.shape[0]:,} genes, PCA {len(pca_df)} samples")


if __name__ == "__main__":
    main()
