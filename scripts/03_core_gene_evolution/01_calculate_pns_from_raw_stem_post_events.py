#!/usr/bin/env python3
"""Step 01: calculate mutation-only gene-level pN/pS before and after DCC emergence."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd


if len(sys.argv) != 3:
    print(
        """
Usage:
python calculate_pns_from_raw_stem_post_events.py input_dir output_dir

Input:
  input_dir should contain DCC*.raw_stem_post_events.csv files.

Output:
  <DCC>.pNS_by_gene.csv
  all_DCC_pNS_by_gene_merged.csv
  DCC1-7_total.pNS_by_gene.csv

Notes:
  - Only rows with region == mutation are used.
  - Recombination rows are excluded because pNS is not meaningful for r events.
  - pNS follows pNS_setsynto1.py: if observed synonymous count is 0,
    observed synonymous is set to 1 for the pNS denominator.
"""
    )
    sys.exit(1)


input_dir = Path(sys.argv[1])
output_dir = Path(sys.argv[2])
output_dir.mkdir(parents=True, exist_ok=True)


MUT_PROB = {
    "A": {"C": 0.297619047619, "G": 0.678571428571, "T": 0.0238095238095},
    "C": {"A": 0.116666666667, "G": 0.238888888889, "T": 0.644444444444},
    "G": {"A": 0.619565217391, "C": 0.239130434783, "T": 0.141304347826},
    "T": {"A": 0.0285714285714, "C": 0.7, "G": 0.271428571429},
}

CODON_TABLE = {
    "AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N", "ACA": "T",
    "ACC": "T", "ACG": "T", "ACT": "T", "AGA": "R", "AGC": "S",
    "AGG": "R", "AGT": "S", "ATA": "I", "ATC": "I", "ATG": "M",
    "ATT": "I", "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P", "CGA": "R",
    "CGC": "R", "CGG": "R", "CGT": "R", "CTA": "L", "CTC": "L",
    "CTG": "L", "CTT": "L", "GAA": "E", "GAC": "D", "GAG": "E",
    "GAT": "D", "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G", "GTA": "V",
    "GTC": "V", "GTG": "V", "GTT": "V", "TAA": "STOP",
    "TAC": "Y", "TAG": "STOP", "TAT": "Y", "TCA": "S",
    "TCC": "S", "TCG": "S", "TCT": "S", "TGA": "STOP",
    "TGC": "C", "TGG": "W", "TGT": "C", "TTA": "L", "TTC": "F",
    "TTG": "L", "TTT": "F",
}


def expected_nonsyn_probability(codon):
    codon = str(codon).upper().replace("U", "T")
    if codon not in CODON_TABLE or CODON_TABLE[codon] == "STOP":
        return np.nan

    expected = 0.0
    for pos, ref_base in enumerate(codon):
        if ref_base not in MUT_PROB:
            return np.nan

        for alt_base, prob in MUT_PROB[ref_base].items():
            alt_codon = list(codon)
            alt_codon[pos] = alt_base
            alt_codon = "".join(alt_codon)

            if alt_codon not in CODON_TABLE:
                continue

            is_nonsyn = CODON_TABLE[codon] != CODON_TABLE[alt_codon]
            expected += (1.0 / 3.0) * prob * int(is_nonsyn)

    return expected


NSY_EXPECTATION = {
    codon: expected_nonsyn_probability(codon)
    for codon in CODON_TABLE
}


def clean_raw_events(df):
    required = {"DCC", "stage", "gene", "region", "mutation_type", "ref_codon"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError("Missing required columns: " + ", ".join(sorted(missing)))

    out = df.copy()
    out["DCC"] = out["DCC"].astype(str).str.strip()
    out["stage"] = out["stage"].astype(str).str.strip()
    out["gene"] = out["gene"].astype(str).str.strip()
    out["region"] = out["region"].astype(str).str.strip().str.lower()
    out["mutation_type"] = out["mutation_type"].astype(str).str.strip().str.lower()
    out["ref_codon"] = out["ref_codon"].astype(str).str.upper().str.replace("U", "T")

    out = out[
        (out["region"] == "mutation")
        & (out["mutation_type"].isin(["synonymous", "nonsynonymous"]))
        & (out["ref_codon"].isin(NSY_EXPECTATION))
    ].copy()
    out["expected_nonsyn_site"] = out["ref_codon"].map(NSY_EXPECTATION)
    out = out[out["expected_nonsyn_site"].notna()].copy()
    return out


def calculate_pns(df, group_cols):
    rows = []

    for keys, sub in df.groupby(group_cols, dropna=False):
        if not isinstance(keys, tuple):
            keys = (keys,)

        observed_nonsyn = int((sub["mutation_type"] == "nonsynonymous").sum())
        observed_syn_raw = int((sub["mutation_type"] == "synonymous").sum())
        total_mut_raw = observed_nonsyn + observed_syn_raw

        observed_syn_for_pns = observed_syn_raw
        syn_set_to_1 = False
        if observed_syn_for_pns == 0:
            observed_syn_for_pns = 1
            syn_set_to_1 = True

        expected_nonsyn = float(sub["expected_nonsyn_site"].sum())
        expected_syn = float((observed_nonsyn + observed_syn_for_pns) - expected_nonsyn)

        pns = np.nan
        if (
            expected_nonsyn > 0
            and expected_syn > 0
            and observed_syn_for_pns > 0
        ):
            pns = (observed_nonsyn / expected_nonsyn) / (
                observed_syn_for_pns / expected_syn
            )
            if pns == 0 or np.isinf(pns):
                pns = np.nan

        row = dict(zip(group_cols, keys))
        row.update(
            {
                "pNS": pns,
                "observed_syn_raw": observed_syn_raw,
                "observed_syn_for_pNS": observed_syn_for_pns,
                "observed_nonsyn": observed_nonsyn,
                "expected_syn": expected_syn,
                "expected_nonsyn": expected_nonsyn,
                "total_mut_raw": total_mut_raw,
                "syn_set_to_1": syn_set_to_1,
            }
        )
        rows.append(row)

    out = pd.DataFrame(rows)
    return out


raw_files = sorted(input_dir.glob("DCC*.raw_stem_post_events.csv"))
raw_files = [x for x in raw_files if not x.name.startswith("DCC1-7")]

if not raw_files:
    sys.exit(f"Error: no DCC*.raw_stem_post_events.csv files found in {input_dir}")

print(f"[INFO] raw DCC event files: {len(raw_files)}")

all_clean = []
per_dcc_files = []

for raw_file in raw_files:
    print(f"[INFO] reading {raw_file.name}")
    raw = pd.read_csv(raw_file)
    clean = clean_raw_events(raw)
    all_clean.append(clean)

    dcc = str(clean["DCC"].iloc[0]) if clean.shape[0] > 0 else raw_file.name.split(".")[0]
    pns = calculate_pns(clean, ["DCC", "stage", "gene"])
    pns = pns.sort_values(["DCC", "stage", "pNS", "total_mut_raw"], ascending=[True, True, False, False])

    out_file = output_dir / f"{dcc}.pNS_by_gene.csv"
    pns.to_csv(out_file, index=False, na_rep="NA")
    per_dcc_files.append(out_file)
    print(f"[INFO] {dcc}: mutation rows used={clean.shape[0]}, pNS rows={pns.shape[0]}")

all_clean_df = pd.concat(all_clean, ignore_index=True)

merged = calculate_pns(all_clean_df, ["DCC", "stage", "gene"])
merged = merged.sort_values(["DCC", "stage", "gene"])
merged_file = output_dir / "all_DCC_pNS_by_gene_merged.csv"
merged.to_csv(merged_file, index=False, na_rep="NA")

total = calculate_pns(all_clean_df, ["stage", "gene"])
total.insert(0, "DCC", "DCC1-7_total")
total = total.sort_values(["stage", "pNS", "total_mut_raw"], ascending=[True, False, False])
total_file = output_dir / "DCC1-7_total.pNS_by_gene.csv"
total.to_csv(total_file, index=False, na_rep="NA")

print("\n==============================")
print("Finished")
print("==============================")
print("Output:")
for x in per_dcc_files:
    print(x)
print(merged_file)
print(total_file)
