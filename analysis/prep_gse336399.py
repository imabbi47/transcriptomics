#!/usr/bin/env python3
"""Build the Cancer-vs-Benign design for GSE336399 and check confounding."""
from __future__ import annotations

import pandas as pd

d = pd.read_csv("/tmp/GSE336399_design_full.csv")
d = d.rename(columns={
    "03_cancer_status": "cancer_status",
    "05_gender": "gender",
    "07_sequencing_batch": "sequencing_batch",
})

d = d[d["cancer_status"].isin(["Cancer", "Benign"])].copy()
d["sequencing_batch"] = "b" + d["sequencing_batch"].astype(str)  # keep categorical

out = d[["sample", "cancer_status", "sequencing_batch", "gender"]]
out.to_csv("/tmp/GSE336399_design.csv", index=False)

print("cancer_status:", out["cancer_status"].value_counts().to_dict())
print("gender:", out["gender"].value_counts().to_dict())
print("\nbatch x cancer_status (checking confounding):")
ct = pd.crosstab(d["sequencing_batch"], d["cancer_status"])
print(ct.to_string())
# a batch is usable only if it has both groups; flag single-group batches
single = ct[(ct == 0).any(axis=1)]
print(f"\nbatches with only one group (drop-if-confounded): {list(single.index) if len(single) else 'none'}")
print(f"wrote /tmp/GSE336399_design.csv ({len(out)} samples)")
