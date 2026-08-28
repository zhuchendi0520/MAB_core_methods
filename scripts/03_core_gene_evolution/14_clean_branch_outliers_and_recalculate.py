#!/usr/bin/env python3
"""
Clean cns_recalculated r/m and pNS results by removing DCC-specific outlier branches.

The outlier table must contain at least:
  DCC,node

For each <DCC>.raw_stem_post_events.csv, rows are removed only when:
  stage == "post" AND (DCC, node) is present in the outlier table

The script then rebuilds:
  - cleaned <DCC>.raw_stem_post_events.csv
  - cleaned <DCC>.gene_stem_post_summary.csv
  - all_DCC_gene_stem_post_summary_wide.csv
  - DCC1-7_total.gene_stem_post_summary.csv
  - pNS outputs in pns/
"""

from __future__ import annotations

import math
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def usage() -> None:
    print(
        "Usage: python clean_cns_recalculated_by_outliers.py "
        "cns_recalculated_dir outlier.csv output_dir "
        "calculate_pns_from_raw_stem_post_events.py "
        "add_rm_pre_post_fisher_stats_csv.py "
        "pns_long_to_pre_post_wide.py"
    )


if len(sys.argv) != 7:
    usage()
    sys.exit(1)

input_dir = Path(sys.argv[1])
outlier_csv = Path(sys.argv[2])
output_dir = Path(sys.argv[3])
pns_script = Path(sys.argv[4])
rm_stats_script = Path(sys.argv[5])
pns_wide_script = Path(sys.argv[6])

output_dir.mkdir(parents=True, exist_ok=True)
pns_dir = output_dir / "pns"
pns_dir.mkdir(parents=True, exist_ok=True)


def clean_text(x):
    return str(x).strip()


def safe_ratio(num, den):
    return num / den if den else np.nan


def count_stage(sub, dcc, stage):
    columns = [
        "DCC", "gene", "stage",
        "rSNP", "mSNP", "r_syn", "r_nonsyn", "m_syn", "m_nonsyn",
        "synonymous", "nonsynonymous", "total_events", "rm_raw", "NS_S_raw",
    ]
    if sub.shape[0] == 0:
        return pd.DataFrame(columns=columns)

    tmp = sub.copy()
    tmp["is_recombination"] = (tmp["region"] == "recombination").astype(int)
    tmp["is_mutation"] = (tmp["region"] == "mutation").astype(int)
    tmp["is_synonymous"] = (tmp["mutation_type"] == "synonymous").astype(int)
    tmp["is_nonsynonymous"] = (tmp["mutation_type"] == "nonsynonymous").astype(int)
    tmp["r_syn_event"] = (
        (tmp["region"] == "recombination")
        & (tmp["mutation_type"] == "synonymous")
    ).astype(int)
    tmp["r_nonsyn_event"] = (
        (tmp["region"] == "recombination")
        & (tmp["mutation_type"] == "nonsynonymous")
    ).astype(int)
    tmp["m_syn_event"] = (
        (tmp["region"] == "mutation")
        & (tmp["mutation_type"] == "synonymous")
    ).astype(int)
    tmp["m_nonsyn_event"] = (
        (tmp["region"] == "mutation")
        & (tmp["mutation_type"] == "nonsynonymous")
    ).astype(int)

    out = (
        tmp.groupby("gene")
        .agg(
            rSNP=("is_recombination", "sum"),
            mSNP=("is_mutation", "sum"),
            r_syn=("r_syn_event", "sum"),
            r_nonsyn=("r_nonsyn_event", "sum"),
            m_syn=("m_syn_event", "sum"),
            m_nonsyn=("m_nonsyn_event", "sum"),
            synonymous=("is_synonymous", "sum"),
            nonsynonymous=("is_nonsynonymous", "sum"),
            total_events=("region", "size"),
        )
        .reset_index()
    )
    out["DCC"] = dcc
    out["stage"] = stage
    out["rm_raw"] = np.where(out["mSNP"] > 0, out["rSNP"] / out["mSNP"], np.nan)
    out["NS_S_raw"] = np.where(
        out["synonymous"] > 0,
        out["nonsynonymous"] / out["synonymous"],
        np.nan,
    )
    return out[columns]


outliers = pd.read_csv(outlier_csv)
required_outlier = {"DCC", "node"}
missing = required_outlier - set(outliers.columns)
if missing:
    raise SystemExit("outlier.csv missing columns: " + ", ".join(sorted(missing)))

outliers["DCC"] = outliers["DCC"].map(clean_text)
outliers["node"] = outliers["node"].map(clean_text)
outlier_pairs = set(zip(outliers["DCC"], outliers["node"]))
outlier_nodes_by_dcc = (
    outliers.groupby("DCC")["node"].apply(lambda x: set(map(clean_text, x))).to_dict()
)

raw_files = sorted(input_dir.glob("DCC*.raw_stem_post_events.csv"))
raw_files = [p for p in raw_files if not p.name.startswith("DCC1-7")]
if not raw_files:
    raise SystemExit(f"No per-DCC raw event files found in {input_dir}")

clean_raw_frames = []
stage_rows = []
report_rows = []

for raw_file in raw_files:
    raw = pd.read_csv(raw_file)
    required_raw = {"DCC", "stage", "node", "gene", "region", "mutation_type"}
    missing_raw = required_raw - set(raw.columns)
    if missing_raw:
        raise SystemExit(f"{raw_file.name} missing columns: {sorted(missing_raw)}")

    raw["DCC"] = raw["DCC"].map(clean_text)
    raw["stage"] = raw["stage"].map(clean_text)
    raw["node"] = raw["node"].map(clean_text)
    raw["gene"] = raw["gene"].map(clean_text)
    raw["region"] = raw["region"].map(clean_text).str.lower()
    raw["mutation_type"] = raw["mutation_type"].map(clean_text).str.lower()

    dcc = raw["DCC"].iloc[0] if raw.shape[0] else raw_file.name.split(".")[0]
    dcc_outlier_nodes = outlier_nodes_by_dcc.get(dcc, set())
    remove_mask = (
        (raw["stage"] == "post")
        & raw["node"].isin(dcc_outlier_nodes)
    )
    removed = raw[remove_mask].copy()
    clean = raw[~remove_mask].copy()

    clean.to_csv(output_dir / raw_file.name, index=False, na_rep="NA")
    clean_raw_frames.append(clean)

    for stage in ["pre", "post"]:
        stage_rows.append(count_stage(clean[clean["stage"] == stage], dcc, stage))

    report_rows.append(
        {
            "DCC": dcc,
            "outlier_node_n": len(dcc_outlier_nodes),
            "input_rows": raw.shape[0],
            "removed_rows": removed.shape[0],
            "removed_mutation_rows": int((removed["region"] == "mutation").sum()),
            "removed_recombination_rows": int((removed["region"] == "recombination").sum()),
            "output_rows": clean.shape[0],
            "removed_unique_nodes": removed["node"].nunique() if removed.shape[0] else 0,
            "removed_unique_genes": removed["gene"].nunique() if removed.shape[0] else 0,
        }
    )

report = pd.DataFrame(report_rows).sort_values("DCC")
report.to_csv(output_dir / "outlier_cleaning_report.csv", index=False)

gene_stage = pd.concat(stage_rows, ignore_index=True)
if gene_stage.shape[0] == 0:
    raise SystemExit("No events remained after cleaning.")

wide = gene_stage.pivot_table(
    index=["DCC", "gene"],
    columns="stage",
    values=[
        "rSNP", "mSNP", "r_syn", "r_nonsyn", "m_syn", "m_nonsyn",
        "synonymous", "nonsynonymous", "total_events",
        "rm_raw", "NS_S_raw",
    ],
    fill_value=0,
    aggfunc="sum",
)
wide.columns = [f"{metric}_{stage}" for metric, stage in wide.columns]
wide = wide.reset_index()

count_cols = [
    "rSNP_pre", "mSNP_pre", "r_syn_pre", "r_nonsyn_pre",
    "m_syn_pre", "m_nonsyn_pre", "synonymous_pre", "nonsynonymous_pre",
    "total_events_pre",
    "rSNP_post", "mSNP_post", "r_syn_post", "r_nonsyn_post",
    "m_syn_post", "m_nonsyn_post", "synonymous_post", "nonsynonymous_post",
    "total_events_post",
]
ratio_cols = ["rm_raw_pre", "rm_raw_post", "NS_S_raw_pre", "NS_S_raw_post"]
for col in count_cols:
    if col not in wide.columns:
        wide[col] = 0
for col in ratio_cols:
    if col not in wide.columns:
        wide[col] = np.nan

wide["rm_raw_pre"] = np.where(wide["mSNP_pre"] > 0, wide["rSNP_pre"] / wide["mSNP_pre"], np.nan)
wide["rm_raw_post"] = np.where(wide["mSNP_post"] > 0, wide["rSNP_post"] / wide["mSNP_post"], np.nan)
wide["NS_S_raw_pre"] = np.where(
    wide["synonymous_pre"] > 0,
    wide["nonsynonymous_pre"] / wide["synonymous_pre"],
    np.nan,
)
wide["NS_S_raw_post"] = np.where(
    wide["synonymous_post"] > 0,
    wide["nonsynonymous_post"] / wide["synonymous_post"],
    np.nan,
)
wide["all_total_events"] = wide["total_events_pre"] + wide["total_events_post"]

output_cols = [
    "DCC", "gene",
    "rSNP_pre", "mSNP_pre", "rm_raw_pre",
    "r_syn_pre", "r_nonsyn_pre", "m_syn_pre", "m_nonsyn_pre",
    "synonymous_pre", "nonsynonymous_pre", "NS_S_raw_pre",
    "total_events_pre",
    "rSNP_post", "mSNP_post", "rm_raw_post",
    "r_syn_post", "r_nonsyn_post", "m_syn_post", "m_nonsyn_post",
    "synonymous_post", "nonsynonymous_post", "NS_S_raw_post",
    "total_events_post",
    "all_total_events",
]
wide = wide[[c for c in output_cols if c in wide.columns]]
wide = wide.sort_values(
    ["DCC", "all_total_events", "total_events_pre", "total_events_post"],
    ascending=[True, False, False, False],
)

for dcc, sub in wide.groupby("DCC"):
    sub.to_csv(output_dir / f"{dcc}.gene_stem_post_summary.csv", index=False, na_rep="NA")

summary_cols = [c for c in output_cols if c not in ["DCC", "gene"]]
count_like_cols = [
    c for c in summary_cols
    if c.startswith(("rSNP_", "mSNP_", "r_syn_", "r_nonsyn_", "m_syn_", "m_nonsyn_"))
    or c.startswith(("synonymous_", "nonsynonymous_", "total_events_"))
    or c == "all_total_events"
]
genes = sorted(wide["gene"].unique())
wide_merged_parts = [pd.DataFrame({"gene": genes})]
for dcc in sorted(wide["DCC"].unique()):
    sub = wide[wide["DCC"] == dcc].set_index("gene")
    dcc_part = pd.DataFrame(index=genes)
    for col in summary_cols:
        if col in count_like_cols:
            dcc_part[f"{dcc}_{col}"] = sub[col].reindex(genes).fillna(0).astype(int)
        else:
            dcc_part[f"{dcc}_{col}"] = sub[col].reindex(genes)
    wide_merged_parts.append(dcc_part.reset_index(drop=True))
wide_merged = pd.concat(wide_merged_parts, axis=1)
wide_merged.to_csv(output_dir / "all_DCC_gene_stem_post_summary_wide.csv", index=False, na_rep="NA")

sum_cols = [
    "rSNP_pre", "mSNP_pre", "r_syn_pre", "r_nonsyn_pre",
    "m_syn_pre", "m_nonsyn_pre", "synonymous_pre", "nonsynonymous_pre",
    "total_events_pre",
    "rSNP_post", "mSNP_post", "r_syn_post", "r_nonsyn_post",
    "m_syn_post", "m_nonsyn_post", "synonymous_post", "nonsynonymous_post",
    "total_events_post", "all_total_events",
]
total = wide.groupby("gene", as_index=False)[sum_cols].sum()
total.insert(0, "DCC", "DCC1-7_total")
total.insert(
    1,
    "n_DCC_with_gene_events",
    wide.groupby("gene")["DCC"].nunique().reindex(total["gene"]).to_numpy(),
)
total["rm_raw_pre"] = np.where(total["mSNP_pre"] > 0, total["rSNP_pre"] / total["mSNP_pre"], np.nan)
total["rm_raw_post"] = np.where(total["mSNP_post"] > 0, total["rSNP_post"] / total["mSNP_post"], np.nan)
total["NS_S_raw_pre"] = np.where(
    total["synonymous_pre"] > 0,
    total["nonsynonymous_pre"] / total["synonymous_pre"],
    np.nan,
)
total["NS_S_raw_post"] = np.where(
    total["synonymous_post"] > 0,
    total["nonsynonymous_post"] / total["synonymous_post"],
    np.nan,
)
total_cols = [
    "DCC", "n_DCC_with_gene_events", "gene",
    "rSNP_pre", "mSNP_pre", "rm_raw_pre",
    "r_syn_pre", "r_nonsyn_pre", "m_syn_pre", "m_nonsyn_pre",
    "synonymous_pre", "nonsynonymous_pre", "NS_S_raw_pre",
    "total_events_pre",
    "rSNP_post", "mSNP_post", "rm_raw_post",
    "r_syn_post", "r_nonsyn_post", "m_syn_post", "m_nonsyn_post",
    "synonymous_post", "nonsynonymous_post", "NS_S_raw_post",
    "total_events_post",
    "all_total_events",
]
total = total[total_cols].sort_values(
    ["all_total_events", "total_events_pre", "total_events_post"],
    ascending=False,
)
total_file = output_dir / "DCC1-7_total.gene_stem_post_summary.csv"
total.to_csv(total_file, index=False, na_rep="NA")

shutil.copy2(rm_stats_script, output_dir / rm_stats_script.name)
subprocess.run(
    [
        sys.executable,
        str(rm_stats_script),
        str(total_file),
        str(output_dir / "DCC1-7_total.gene_stem_post_summary_with_rm_stats.csv"),
    ],
    check=True,
)

shutil.copy2(pns_script, pns_dir / pns_script.name)
shutil.copy2(pns_wide_script, pns_dir / pns_wide_script.name)
subprocess.run([sys.executable, str(pns_script), str(output_dir), str(pns_dir)], check=True)
subprocess.run(
    [
        sys.executable,
        str(pns_wide_script),
        str(pns_dir / "DCC1-7_total.pNS_by_gene.csv"),
        str(pns_dir / "DCC1-7_total.pNS_by_gene_pre_post_wide.csv"),
    ],
    check=True,
)
subprocess.run(
    [
        sys.executable,
        str(pns_wide_script),
        str(pns_dir / "all_DCC_pNS_by_gene_merged.csv"),
        str(pns_dir / "all_DCC_pNS_by_gene_merged_pre_post_wide.csv"),
    ],
    check=True,
)

shutil.copy2(Path(__file__), output_dir / Path(__file__).name)

print("Finished outlier cleaning")
print(output_dir)
print(output_dir / "outlier_cleaning_report.csv")
print(output_dir / "DCC1-7_total.gene_stem_post_summary_with_rm_stats.csv")
print(pns_dir / "DCC1-7_total.pNS_by_gene.csv")
