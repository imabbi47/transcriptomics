#!/usr/bin/env python3
"""Differential expression from a GEO-style count matrix using pyDESeq2.

Reads a genes x samples count table + a design CSV, runs DESeq2 with an optional
blocking covariate, and writes a results table plus PCA / volcano / MA plots.

    python de.py --counts counts.tsv.gz --design design.csv \
        --factor treatment --ref Control --alt Sildenafil \
        --covariate cell_line --outdir results
"""
from __future__ import annotations

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def load_counts(path):
    """Read a genes x samples table -> samples x genes integer DataFrame."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    numeric = df.apply(pd.to_numeric, errors="coerce")
    coerced = int(numeric.isna().sum().sum() - df.isna().sum().sum())
    if coerced:
        print(f"[de] WARNING: {coerced} non-numeric count value(s) coerced to 0 — "
              "check the matrix for stray text columns/rows before trusting results")
    df = numeric.dropna(how="all").fillna(0)
    counts = df.T.round(0).astype(int)
    counts.index = [str(i) for i in counts.index]
    return counts


def build_dds(counts, meta, design_str, factors):
    from pydeseq2.dds import DeseqDataSet
    errs = []
    for kwargs in ({"design": design_str}, {"design_factors": factors}):
        try:
            return DeseqDataSet(counts=counts, metadata=meta, **kwargs)
        except TypeError as exc:
            errs.append(repr(exc))
    raise RuntimeError("DeseqDataSet construction failed: " + " | ".join(errs))


def extract_symbols(gene_ids):
    """Best-effort gene symbols from row IDs.

    Handles IDs that embed a symbol (e.g. ENSMUSG..._Gnai3) and bare Ensembl IDs
    (e.g. ENSG...), the latter mapped via mygene when available.
    """
    ids = [str(g) for g in gene_ids]
    if any("_" in g for g in ids[:200]):
        return [g.split("_", 1)[-1] for g in ids]
    head = ids[0] if ids else ""
    if head.startswith(("ENSG", "ENSMUSG", "ENST", "ENSMUST")):
        species = "human" if head.startswith(("ENSG", "ENST")) else "mouse"
        try:
            import mygene
            clean = [g.split(".")[0] for g in ids]
            hits = mygene.MyGeneInfo().querymany(
                clean, scopes="ensembl.gene", fields="symbol", species=species, verbose=False)
            mapping = {}
            for h in hits:
                q, s = h.get("query"), h.get("symbol")
                if q and s and q not in mapping:
                    mapping[q] = s
            print(f"[de] mapped {len(mapping)}/{len(set(clean))} Ensembl IDs to symbols")
            return [mapping.get(c, orig) for c, orig in zip(clean, ids)]
        except Exception as exc:  # noqa: BLE001
            print(f"[de] symbol mapping skipped: {exc!r}")
    return ids


def pca_plot(dds, meta, args):
    from sklearn.decomposition import PCA
    normed = np.asarray(dds.layers["normed_counts"], dtype=float)
    logn = np.log2(normed + 1.0)
    var = logn.var(axis=0)
    idx = np.argsort(var)[-min(1000, logn.shape[1]):]
    X = logn[:, idx]
    X = X - X.mean(axis=0, keepdims=True)
    pcs = PCA(n_components=2).fit_transform(X)
    samples = [str(s) for s in dds.obs_names]
    fac = meta.loc[samples, args.factor].astype(str).values
    cov_first = next((c.strip() for c in (args.covariate.split(",") if args.covariate else []) if c.strip()), None)
    cov = (meta.loc[samples, cov_first].astype(str).values
           if cov_first and meta[cov_first].nunique() <= 7 else None)

    fig, ax = plt.subplots(figsize=(6.6, 5))
    fac_levels = sorted(set(fac))
    palette = list(plt.cm.Set1.colors)
    color = {lv: palette[i % len(palette)] for i, lv in enumerate(fac_levels)}
    markers = ["o", "s", "^", "D", "v", "P", "X"]
    cov_levels = sorted(set(cov)) if cov is not None else [None]
    mark = {lv: markers[i % len(markers)] for i, lv in enumerate(cov_levels)}

    for i, _ in enumerate(samples):
        ax.scatter(pcs[i, 0], pcs[i, 1], color=color[fac[i]],
                   marker=mark[cov[i]] if cov is not None else "o",
                   s=95, edgecolor="k", linewidth=0.4)

    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], marker="o", color="w", markerfacecolor=color[lv],
                      markersize=10, label=lv) for lv in fac_levels]
    if cov is not None:
        handles += [Line2D([0], [0], marker=mark[lv], color="k", markerfacecolor="w",
                           markersize=9, label=lv) for lv in cov_levels]
    ax.legend(handles=handles, fontsize=8, frameon=False)
    ax.set_xlabel("PC1"), ax.set_ylabel("PC2")
    ax.set_title("PCA (top variable genes)")
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "pca.png"), dpi=130)
    plt.close(fig)


def volcano_plot(res, args):
    d = res.dropna(subset=["padj", "log2FoldChange"]).copy()
    d["nlp"] = -np.log10(d["padj"].clip(lower=1e-300))
    sig = d["padj"] < args.alpha
    fig, ax = plt.subplots(figsize=(6.6, 5))
    ax.scatter(d.loc[~sig, "log2FoldChange"], d.loc[~sig, "nlp"], s=6, c="#bfbfbf", alpha=0.5)
    ax.scatter(d.loc[sig, "log2FoldChange"], d.loc[sig, "nlp"], s=9, c="#c0392b", alpha=0.7)
    for _, r in d.sort_values("padj").head(10).iterrows():
        ax.annotate(str(r["symbol"]), (r["log2FoldChange"], r["nlp"]), fontsize=7)
    ax.axhline(-np.log10(args.alpha), ls="--", c="k", lw=0.6)
    ax.axvline(0, ls=":", c="k", lw=0.5)
    ax.set_xlabel(f"log2 fold-change ({args.alt} vs {args.ref})")
    ax.set_ylabel("-log10 adjusted p")
    ax.set_title("Volcano")
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "volcano.png"), dpi=130)
    plt.close(fig)


def ma_plot(res, args):
    d = res.dropna(subset=["log2FoldChange", "baseMean"]).copy()
    d = d[d["baseMean"] > 0]
    sig = d["padj"].notna() & (d["padj"] < args.alpha)
    fig, ax = plt.subplots(figsize=(6.6, 5))
    ax.scatter(d.loc[~sig, "baseMean"], d.loc[~sig, "log2FoldChange"], s=6, c="#bfbfbf", alpha=0.5)
    ax.scatter(d.loc[sig, "baseMean"], d.loc[sig, "log2FoldChange"], s=9, c="#2c82c9", alpha=0.7)
    ax.set_xscale("log")
    ax.axhline(0, c="k", lw=0.6)
    ax.set_xlabel("mean normalized count")
    ax.set_ylabel("log2 fold-change")
    ax.set_title("MA plot")
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "ma.png"), dpi=130)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--design", required=True)
    ap.add_argument("--factor", default="treatment")
    ap.add_argument("--ref", required=True)
    ap.add_argument("--alt", required=True)
    ap.add_argument("--covariate", default=None)
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--alpha", type=float, default=0.05)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    counts = load_counts(args.counts)
    meta = pd.read_csv(args.design, index_col=0)
    meta.index = [str(i) for i in meta.index]

    common = [s for s in meta.index if s in counts.index]
    counts = counts.loc[common]
    meta = meta.loc[common]
    print(f"[de] aligned: {counts.shape[0]} samples x {counts.shape[1]} genes")
    covariates = [c.strip() for c in (args.covariate.split(",") if args.covariate else []) if c.strip()]
    print(f"[de] {args.factor}: {meta[args.factor].value_counts().to_dict()}")
    for cv in covariates:
        print(f"[de] {cv}: {meta[cv].value_counts().to_dict()}")

    factor_levels = set(meta[args.factor].astype(str))
    missing = [x for x in (args.ref, args.alt) if x not in factor_levels]
    if missing:
        raise SystemExit(f"[de] --ref/--alt {missing} not found in '{args.factor}' "
                         f"(levels present: {sorted(factor_levels)})")
    if not common:
        raise SystemExit("[de] no samples shared between the count matrix and the design "
                         "— check that the design 'sample' column matches the count columns")

    keep = counts.sum(axis=0) >= 10
    counts = counts.loc[:, keep]
    print(f"[de] {int(keep.sum())} genes pass total-count>=10 filter")

    factors = covariates + [args.factor]
    meta = meta.copy()
    for f in factors:
        meta[f] = meta[f].astype(str)  # treat every design term as categorical
    design_str = "~" + " + ".join(factors)
    print(f"[de] design {design_str} | contrast {args.factor}: {args.alt} vs {args.ref}")

    dds = build_dds(counts, meta, design_str, factors)
    dds.deseq2()

    from pydeseq2.ds import DeseqStats
    ds = DeseqStats(dds, contrast=[args.factor, args.alt, args.ref])
    ds.summary()
    res = ds.results_df.copy()
    res["symbol"] = extract_symbols(list(res.index))
    res = res.sort_values("padj")
    res.to_csv(os.path.join(args.outdir, "de_results.csv"))

    tested = int(res["padj"].notna().sum())
    sig = res[res["padj"].notna() & (res["padj"] < args.alpha)]
    up = sig[sig["log2FoldChange"] > 0]
    down = sig[sig["log2FoldChange"] < 0]
    print(f"\n[de] === {args.alt} vs {args.ref} (blocking on {args.covariate}) ===")
    print(f"[de] genes tested: {tested:,}")
    print(f"[de] significant padj<{args.alpha}: {len(sig):,}  (up={len(up):,}, down={len(down):,})")
    print("\n[de] Top 15 genes by adjusted p-value:")
    show = res.head(15)[["symbol", "baseMean", "log2FoldChange", "padj"]].copy()
    show["baseMean"] = show["baseMean"].round(0)
    show["log2FoldChange"] = show["log2FoldChange"].round(2)
    print(show.to_string())

    for fn in (pca_plot, volcano_plot, ma_plot):
        try:
            fn(dds, meta, args) if fn is pca_plot else fn(res, args)
        except Exception as exc:  # noqa: BLE001
            print(f"[de] {fn.__name__} failed: {exc!r}")

    print(f"\n[de] wrote de_results.csv + pca.png + volcano.png + ma.png to {args.outdir}/")


if __name__ == "__main__":
    main()
