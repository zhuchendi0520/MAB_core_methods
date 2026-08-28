#!/usr/bin/env python3

import sys
import pandas as pd

# =========================================================
# Usage
# =========================================================
if len(sys.argv) != 3:
    print("""
Usage:
python DCC_gene_category_count.py gene_presence_absence.Rtab sample_group.tsv
""")
    sys.exit(1)

rtab_file = sys.argv[1]
group_file = sys.argv[2]

# =========================================================
# Load data
# =========================================================
print("Loading Rtab...")

rtab = pd.read_csv(
    rtab_file,
    sep="\t",
    index_col=0
)

rtab = rtab.apply(pd.to_numeric)

print(f"Genes: {rtab.shape[0]}")
print(f"Genomes: {rtab.shape[1]}")

# =========================================================
# Load metadata
# =========================================================
meta = pd.read_csv(group_file, sep="\t")

meta = meta[meta["sample"].isin(rtab.columns)]

# =========================================================
# ONLY keep DCC1-7
# =========================================================
dcc_groups = [f"DCC{i}" for i in range(1, 8)]

meta = meta[meta["group"].isin(dcc_groups)]

group_dict = {
    g: meta.loc[meta["group"] == g, "sample"].tolist()
    for g in dcc_groups
}

print("\nDetected DCC groups:")
for g in dcc_groups:
    print(f"{g}: {len(group_dict[g])}")

# =========================================================
# Gene category thresholds
# =========================================================
# Core genes      (99% <= strains <= 100%)
# Soft core genes (95% <= strains < 99%)
# Shell genes     (15% <= strains < 95%)
# Cloud genes     (0% <= strains < 15%)

CORE_CUTOFF = 0.99
SOFT_CORE_CUTOFF = 0.95
SHELL_CUTOFF = 0.15

# =========================================================
# Calculate gene prevalence within each DCC
# =========================================================
print("\nCalculating gene prevalence...")

results = []

for g in dcc_groups:

    strains = group_dict[g]

    if len(strains) == 0:
        continue

    sub = rtab[strains]

    prevalence = sub.sum(axis=1) / len(strains)

    # =====================================================
    # Count categories
    # =====================================================
    core = ((prevalence >= CORE_CUTOFF) &
            (prevalence <= 1.0)).sum()

    soft_core = ((prevalence >= SOFT_CORE_CUTOFF) &
                 (prevalence < CORE_CUTOFF)).sum()

    shell = ((prevalence >= SHELL_CUTOFF) &
             (prevalence < SOFT_CORE_CUTOFF)).sum()

    cloud = ((prevalence > 0) &
             (prevalence < SHELL_CUTOFF)).sum()

    results.append({
        "DCC": g,
        "Core": core,
        "Soft_core": soft_core,
        "Shell": shell,
        "Cloud": cloud,
        "Total_genes": core + soft_core + shell + cloud
    })

# =========================================================
# Save results
# =========================================================
final_df = pd.DataFrame(results)

final_df.to_csv(
    "DCC_gene_category_counts.csv",
    index=False
)

print("\n==============================")
print("Finished")
print("==============================")

print("\nOutput:")
print("DCC_gene_category_counts.csv")

print("\nSummary:")
print(final_df)