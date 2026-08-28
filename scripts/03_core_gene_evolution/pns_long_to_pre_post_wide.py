#!/usr/bin/env python3
"""
Convert pNS_by_gene long table to one-row-per-DCC-gene wide format.

Input columns should include:
  DCC, stage, gene, pNS, observed_syn_raw, observed_nonsyn,
  expected_syn, expected_nonsyn, total_mut_raw, syn_set_to_1

Output:
  DCC, gene, post_pNS, post_observed_syn_raw, ..., pre_pNS, ...
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd


if len(sys.argv) not in (2, 3):
    print("Usage: python pns_long_to_pre_post_wide.py input.pNS_by_gene.csv [output.csv]")
    sys.exit(1)


input_csv = Path(sys.argv[1])
output_csv = Path(sys.argv[2]) if len(sys.argv) == 3 else input_csv.with_name(
    input_csv.stem + "_pre_post_wide.csv"
)

df = pd.read_csv(input_csv)

required = {"DCC", "stage", "gene"}
missing = required - set(df.columns)
if missing:
    raise SystemExit("Missing required columns: " + ", ".join(sorted(missing)))

df["stage"] = df["stage"].astype(str).str.strip()
df = df[df["stage"].isin(["pre", "post"])].copy()

value_cols = [c for c in df.columns if c not in ["DCC", "stage", "gene"]]

wide = df.pivot_table(
    index=["DCC", "gene"],
    columns="stage",
    values=value_cols,
    aggfunc="first",
    dropna=False,
)

wide.columns = [f"{stage}_{metric}" for metric, stage in wide.columns]
wide = wide.reset_index()

ordered_cols = ["DCC", "gene"]
for stage in ["post", "pre"]:
    for metric in value_cols:
        col = f"{stage}_{metric}"
        if col in wide.columns:
            ordered_cols.append(col)

remaining = [c for c in wide.columns if c not in ordered_cols]
wide = wide[ordered_cols + remaining]

count_like = [
    c for c in wide.columns
    if any(x in c for x in [
        "observed_syn_raw", "observed_nonsyn", "total_mut_raw"
    ])
]
for col in count_like:
    wide[col] = wide[col].astype("Int64")

wide.to_csv(output_csv, index=False, na_rep="NA")
print(output_csv)
