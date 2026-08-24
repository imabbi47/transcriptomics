#!/usr/bin/env python3
"""Aggregate per-sample Salmon quant.sf files into a gene x sample count matrix.

Maps transcript -> gene symbol via tx2gene.tsv (transcript_id<TAB>symbol), summing
Salmon's NumReads per gene (tximport-style, no length correction). Versions are
stripped from transcript IDs on both sides so GENCODE-versioned names always match.

    python aggregate.py <quant_root> <tx2gene.tsv> <out_counts.tsv>
"""
import os
import sys
import glob
import pandas as pd

quant_root, tx2gene_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]

t2g = pd.read_csv(tx2gene_path, sep="\t", header=None, names=["tx", "gene"]).dropna()
t2g["txbase"] = t2g["tx"].astype(str).str.split(".").str[0]
tx2gene = dict(zip(t2g["txbase"], t2g["gene"]))
print(f"[agg] tx2gene: {len(tx2gene)} transcripts -> {t2g['gene'].nunique()} genes")

cols = {}
for sf in sorted(glob.glob(os.path.join(quant_root, "*", "quant.sf"))):
    sample = os.path.basename(os.path.dirname(sf))
    q = pd.read_csv(sf, sep="\t")
    q["txbase"] = q["Name"].astype(str).str.split(".").str[0]
    q["gene"] = q["txbase"].map(tx2gene)
    unmapped = int(q["gene"].isna().sum())
    g = q.dropna(subset=["gene"]).groupby("gene")["NumReads"].sum()
    cols[sample] = g
    print(f"[agg] {sample}: {len(q)} tx, {unmapped} unmapped, "
          f"{int(g.gt(0).sum())} genes>0, {q['NumReads'].sum():,.0f} reads assigned")

if not cols:
    sys.exit("[agg] no quant.sf files found under " + quant_root)

mat = pd.DataFrame(cols).fillna(0.0)
mat.index.name = "gene"
mat = mat.round(0).astype(int)
mat = mat.reindex(sorted(mat.columns), axis=1)
mat.to_csv(out_path, sep="\t")
print(f"[agg] wrote {out_path}: {mat.shape[0]} genes x {mat.shape[1]} samples")
print("[agg] samples:", list(mat.columns))
