#!/usr/bin/env python3
"""Look up a gene set across one or more DE result tables (matched by symbol).

    python lookup_geneset.py --genes genes.txt --out out.xlsx \
        --de LABEL1:results/a/de_results.csv --de LABEL2:results/b/de_results.csv

Matching is case-insensitive on the `symbol` column, so human ALL-CAPS symbols
also hit their same-spelling mouse orthologs (e.g. BBS1 -> Bbs1). Orthologs with
divergent names (e.g. TP53 vs Trp53) won't match — that needs a homology map.
"""
from __future__ import annotations

import argparse

import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", required=True, help="text file, one symbol per line")
    ap.add_argument("--de", action="append", required=True, metavar="LABEL:PATH")
    ap.add_argument("--alpha", type=float, default=0.05)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.genes) as fh:
        raw = [g.strip() for g in fh if g.strip()]
    seen = {}
    for g in raw:
        seen.setdefault(g.upper(), g)  # dedupe on upper key, keep original spelling
    q = pd.DataFrame({"gene": list(seen.values()), "KEY": list(seen.keys())})

    summary = []
    for spec in args.de:
        label, path = spec.split(":", 1)
        de = pd.read_csv(path)
        de = de[de["symbol"].notna()].copy()
        de["KEY"] = de["symbol"].astype(str).str.upper()
        de["_p"] = de["padj"].fillna(1.0)
        de = de.sort_values(["KEY", "_p"]).drop_duplicates("KEY", keep="first").set_index("KEY")
        q[f"{label}__found"] = q["KEY"].isin(de.index)
        q[f"{label}__log2FC"] = q["KEY"].map(de["log2FoldChange"]).round(3)
        q[f"{label}__padj"] = q["KEY"].map(de["padj"])
        q[f"{label}__sig"] = q[f"{label}__padj"] < args.alpha
        summary.append({
            "dataset": label,
            "genes_in_set": len(q),
            "found_in_data": int(q[f"{label}__found"].sum()),
            f"significant_padj<{args.alpha}": int(q[f"{label}__sig"].sum()),
        })

    found_cols = [c for c in q.columns if c.endswith("__found")]
    sig_cols = [c for c in q.columns if c.endswith("__sig")]
    q["found_in_all"] = q[found_cols].all(axis=1)
    q["sig_in_all"] = q[sig_cols].all(axis=1)

    summ = pd.DataFrame(summary)
    out = q.drop(columns=["KEY"])
    with pd.ExcelWriter(args.out, engine="openpyxl") as xl:
        summ.to_excel(xl, sheet_name="summary", index=False)
        out.to_excel(xl, sheet_name="geneset_in_data", index=False)
        out[q["sig_in_all"]].to_excel(xl, sheet_name="significant_in_all", index=False)

    print(summ.to_string(index=False))
    print(f"\nfound in ALL datasets   : {int(q['found_in_all'].sum())}")
    print(f"significant in ALL (padj<{args.alpha}): {int(q['sig_in_all'].sum())}")
    print(f"[xlsx] wrote {args.out}")


if __name__ == "__main__":
    main()
