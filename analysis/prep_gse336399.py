#!/usr/bin/env python3
"""Build the Cancer-vs-Benign design for GSE336399 from a geo-design table.

    python prep_gse336399.py [design_full.csv] [design.csv]

Keeps only Cancer + Benign samples (drops Metastatic / Not Collected), and
carries sequencing_batch + gender as covariates. Prints a batch x status
crosstab so you can spot batches confounded with the contrast.
"""
from __future__ import annotations

import sys

import pandas as pd

src = sys.argv[1] if len(sys.argv) > 1 else "results/gse336399/design_full.csv"
dst = sys.argv[2] if len(sys.argv) > 2 else "results/gse336399/design.csv"

d = pd.read_csv(src).rename(columns={
    "03_cancer_status": "cancer_status",
    "05_gender": "gender",
    "07_sequencing_batch": "sequencing_batch",
})
d = d[d["cancer_status"].isin(["Cancer", "Benign"])].copy()
d["sequencing_batch"] = "b" + d["sequencing_batch"].astype(str)  # keep categorical

out = d[["sample", "cancer_status", "sequencing_batch", "gender"]]
out.to_csv(dst, index=False)

print("cancer_status:", out["cancer_status"].value_counts().to_dict())
print("gender:", out["gender"].value_counts().to_dict())
print("\nbatch x cancer_status (confounding check):")
ct = pd.crosstab(d["sequencing_batch"], d["cancer_status"])
print(ct.to_string())
single = ct[(ct == 0).any(axis=1)]
print(f"\nbatches with only one group: {list(single.index) if len(single) else 'none'}")
print(f"wrote {dst} ({len(out)} samples)")
