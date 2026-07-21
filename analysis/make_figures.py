#!/usr/bin/env python3
"""Regenerate report figures at high resolution from existing result files.

Reuses de_results.csv, the count matrix, the design and the enrichr tables — so
no DESeq2/Enrichr re-run is needed. Produces <prefix>_PCA/_volcano/_enrichment.png.
"""
from __future__ import annotations

import argparse
import os

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def _load_counts(p):
    sep = "\t" if p.endswith((".tsv.gz", ".txt.gz", ".tsv", ".txt")) else ","
    return pd.read_csv(p, sep=sep, index_col=0).apply(pd.to_numeric, errors="coerce").dropna(how="all")


def fig_pca(counts_p, design_p, factor, sample_col, title, out, dpi):
    from sklearn.decomposition import PCA
    counts = _load_counts(counts_p)
    design = pd.read_csv(design_p)
    samples = [s for s in design[sample_col].astype(str) if s in counts.columns]
    X = counts[samples]
    lg = np.log2(X / X.sum(0) * 1e6 + 1)
    top = lg.loc[lg.var(1).sort_values(ascending=False).index[:500]]
    m = top.T.values - top.T.values.mean(0)
    pca = PCA(2).fit(m)
    co = pca.transform(m)
    pct = pca.explained_variance_ratio_ * 100
    d = design.copy()
    d[sample_col] = d[sample_col].astype(str)
    grp = d.set_index(sample_col).loc[samples, factor].astype(str).values
    fig, ax = plt.subplots(figsize=(6.6, 5.3))
    for g in sorted(set(grp)):
        mk = grp == g
        ax.scatter(co[mk, 0], co[mk, 1], s=48, alpha=.85, edgecolors="white", linewidths=.4, label=g)
    ax.set_xlabel(f"PC1 ({pct[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({pct[1]:.1f}%)")
    ax.set_title(title)
    ax.legend(fontsize=8, title=factor, frameon=True)
    fig.tight_layout()
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    print("  ", out)


def fig_volcano(de_p, genes_p, title, out, dpi, alpha=.05, lfc=1., nlabel=12):
    gset = {ln.strip().upper() for ln in open(genes_p) if ln.strip()}
    de = pd.read_csv(de_p).dropna(subset=["padj", "log2FoldChange"]).copy()
    de["nlp"] = -np.log10(de["padj"].clip(lower=1e-300))
    de["inset"] = de["symbol"].astype(str).str.upper().isin(gset)
    bg, st = de[~de.inset], de[de.inset]
    sig = st[st.padj < alpha]
    nsig = st[st.padj >= alpha]
    up, dn = sig[sig.log2FoldChange > 0], sig[sig.log2FoldChange < 0]
    fig, ax = plt.subplots(figsize=(7.2, 6))
    ax.scatter(bg.log2FoldChange, bg.nlp, s=5, c="#d9dce1", alpha=.5, linewidths=0, label="all genes")
    ax.scatter(nsig.log2FoldChange, nsig.nlp, s=18, c="#f0a24b", alpha=.8, linewidths=0, label="gene set (ns)")
    ax.scatter(up.log2FoldChange, up.nlp, s=26, c="#c0392b", edgecolors="white", linewidths=.3, label="set up (sig)")
    ax.scatter(dn.log2FoldChange, dn.nlp, s=26, c="#2c6fb0", edgecolors="white", linewidths=.3, label="set down (sig)")
    ax.axhline(-np.log10(alpha), ls="--", c="grey", lw=.7)
    ax.axvline(lfc, ls="--", c="grey", lw=.7)
    ax.axvline(-lfc, ls="--", c="grey", lw=.7)
    for _, r in sig.sort_values("padj").head(nlabel).iterrows():
        ax.annotate(str(r["symbol"]), (r.log2FoldChange, r.nlp), fontsize=7, xytext=(3, 3), textcoords="offset points")
    xm = min(10., max(3., float(np.ceil(np.percentile(np.abs(de.log2FoldChange), 99)))))
    ax.set_xlim(-xm, xm)
    sm = float(st.nlp.max()) if len(st) else float(de.nlp.max())
    ym = max(5., float(np.ceil(sm * 1.15)))
    ax.set_ylim(-ym * .03, ym)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p-value")
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    print("  ", out)


def fig_enrichment(up_p, dn_p, title, out, dpi, n=8):
    fig, axs = plt.subplots(2, 1, figsize=(7.6, 8))
    for ax, path, ttl, color in [(axs[0], up_p, "Up-regulated", "#c0392b"),
                                 (axs[1], dn_p, "Down-regulated", "#2c6fb0")]:
        if os.path.exists(path):
            e = pd.read_csv(path).sort_values("Adjusted P-value").head(n).iloc[::-1]
            y = list(range(len(e)))
            ax.barh(y, -np.log10(e["Adjusted P-value"].clip(lower=1e-300)), color=color, alpha=.85)
            ax.set_yticks(y)
            ax.set_yticklabels([str(t)[:44] for t in e["Term"]], fontsize=7)
            ax.set_xlabel("-log10 adjusted p")
            ax.set_title(ttl, fontsize=10, loc="left")
        else:
            ax.axis("off")
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    print("  ", out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--de", required=True)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--design", required=True)
    ap.add_argument("--factor", required=True)
    ap.add_argument("--genes", required=True)
    ap.add_argument("--prefix", required=True)
    ap.add_argument("--title", required=True)
    ap.add_argument("--sample-col", default="sample")
    ap.add_argument("--dpi", type=int, default=300)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    edir = os.path.dirname(a.de)
    print(f"[figures] {a.prefix} @ {a.dpi} dpi:")
    fig_pca(a.counts, a.design, a.factor, a.sample_col, a.title + " — PCA",
            os.path.join(a.outdir, a.prefix + "_PCA.png"), a.dpi)
    fig_volcano(a.de, a.genes, a.title + " — cilia gene set",
                os.path.join(a.outdir, a.prefix + "_volcano.png"), a.dpi)
    fig_enrichment(os.path.join(edir, "enrichr_up.csv"), os.path.join(edir, "enrichr_down.csv"),
                   a.title + " — pathway enrichment",
                   os.path.join(a.outdir, a.prefix + "_enrichment.png"), a.dpi)


if __name__ == "__main__":
    main()
