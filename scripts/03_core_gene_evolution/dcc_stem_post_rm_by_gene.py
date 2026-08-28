#!/usr/bin/env python3

import sys
from pathlib import Path
import pandas as pd
import numpy as np
from ete3 import Tree


# =========================================================
# Usage
# =========================================================
if len(sys.argv) != 5:
    print("""
Usage:
python dcc_stem_post_rm_by_gene.py tree.nwk mutation_detail.tsv target_nodes.txt output_dir

Input:
1. tree.nwk
   - node-labelled tree, readable by ete3 with format=1

2. mutation_detail.tsv
   - required columns: node, gene, region
   - region should contain: recombination or mutation
   - each row is one SNP/event assigned to a branch/node

3. target_nodes.txt
   - one target node per line, two columns, or three columns:
     node_name
     DCC_name<TAB>post_root_node
     DCC_name<TAB>post_root_node<TAB>stem_node1,stem_node2
   - if one column is supplied, DCC_name and post_root_node will both
     be the supplied node name
   - if the third column is absent, stem branch is post_root_node itself

Definition:
stem branch:
   events assigned to the stem node(s); by default this is the
   post_root_node itself

post-DCC branches:
   events assigned to all named descendant branches below the post_root_node,
   excluding the post_root_node itself

Output:
<DCC>.raw_stem_post_events.csv
<DCC>.gene_stem_post_summary.csv
all_DCC_gene_stem_post_summary_wide.csv
DCC1-7_total.gene_stem_post_summary.csv
""")
    sys.exit(1)

tree_file = sys.argv[1]
mut_file = sys.argv[2]
target_file = sys.argv[3]
output_dir = Path(sys.argv[4])
output_dir.mkdir(parents=True, exist_ok=True)


# =========================================================
# Helpers
# =========================================================
def clean_node_name(x):
    return str(x).strip()


def read_target_nodes(path):
    targets = []

    with open(path) as f:
        for line in f:
            line = line.strip()

            if line == "" or line.startswith("#"):
                continue

            parts = line.split()

            if len(parts) == 1:
                dcc = parts[0]
                node = parts[0]
                stem_nodes = [node]
            elif len(parts) == 2:
                dcc = parts[0]
                node = parts[1]
                stem_nodes = [node]
            else:
                dcc = parts[0]
                node = parts[1]
                stem_nodes = [
                    x.strip() for x in parts[2].split(",")
                    if x.strip() != ""
                ]

            targets.append((dcc, node, stem_nodes))

    return targets


def count_stage(sub, dcc, target_node, stage):
    if sub.shape[0] == 0:
        return pd.DataFrame(columns=[
            "DCC", "target_node", "gene", "stage",
            "rSNP", "mSNP", "r_syn", "r_nonsyn", "m_syn", "m_nonsyn",
            "synonymous", "nonsynonymous", "total_events", "rm_raw",
            "NS_S_raw"
        ])

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
            total_events=("region", "size")
        )
        .reset_index()
    )

    out["DCC"] = dcc
    out["target_node"] = target_node
    out["stage"] = stage
    out["rm_raw"] = np.where(
        out["mSNP"] > 0,
        out["rSNP"] / out["mSNP"],
        np.nan
    )
    out["NS_S_raw"] = np.where(
        out["synonymous"] > 0,
        out["nonsynonymous"] / out["synonymous"],
        np.nan
    )

    out = out[[
        "DCC", "target_node", "gene", "stage",
        "rSNP", "mSNP", "r_syn", "r_nonsyn", "m_syn", "m_nonsyn",
        "synonymous", "nonsynonymous", "total_events", "rm_raw",
        "NS_S_raw"
    ]]

    return out


# =========================================================
# Load tree
# =========================================================
print("[INFO] loading tree...")

T = Tree(tree_file, format=1)

tree_node_names = set()
for n in T.traverse():
    if n.name.strip() != "":
        tree_node_names.add(n.name.strip())

print(f"[INFO] named tree nodes: {len(tree_node_names)}")


# =========================================================
# Load mutation table
# =========================================================
print("[INFO] loading mutation table...")

df = pd.read_csv(mut_file, sep="\t")

required_cols = {"node", "gene", "region"}
missing_cols = required_cols - set(df.columns)

if len(missing_cols) > 0:
    sys.exit(
        "Error: mutation_detail.tsv missing required columns: "
        + ", ".join(sorted(missing_cols))
    )

df["node"] = df["node"].astype(str).str.strip()
df["gene"] = df["gene"].astype(str).str.strip()
df["region"] = df["region"].astype(str).str.strip().str.lower()

if "mutation_type" not in df.columns:
    df["mutation_type"] = "unknown"
else:
    df["mutation_type"] = df["mutation_type"].astype(str).str.strip().str.lower()

df = df[df["region"].isin(["recombination", "mutation"])].copy()

print(f"[INFO] usable mutation/recombination rows: {df.shape[0]}")
print(f"[INFO] genes in mutation table: {df['gene'].nunique()}")


# =========================================================
# Load targets
# =========================================================
targets = read_target_nodes(target_file)

if len(targets) == 0:
    sys.exit("Error: no target nodes found")

print(f"[INFO] target nodes: {len(targets)}")


# =========================================================
# Extract stem and post-DCC branches
# =========================================================
stage_rows = []
raw_files = []

for dcc, target_name, target_stem_nodes in targets:
    target_name = clean_node_name(target_name)

    try:
        target = T & target_name
    except Exception:
        print(f"[WARNING] target node not found in tree: {target_name}")
        continue

    stem_nodes = set()
    for stem_name in target_stem_nodes:
        stem_name = clean_node_name(stem_name)
        if stem_name == "":
            continue

        try:
            stem_node = T & stem_name
        except Exception:
            print(f"[WARNING] stem node not found in tree: {stem_name}")
            continue

        if stem_node.name.strip() != "":
            stem_nodes.add(stem_node.name.strip())

    post_nodes = set()
    for n in target.traverse():
        name = n.name.strip()
        if name == "" or name == target_name:
            continue
        post_nodes.add(name)

    stem_df = df[df["node"].isin(stem_nodes)].copy()
    post_df = df[df["node"].isin(post_nodes)].copy()

    raw_stem = stem_df.copy()
    raw_post = post_df.copy()
    raw_stem.insert(0, "stage", "pre")
    raw_post.insert(0, "stage", "post")
    raw_stem.insert(0, "target_node", target_name)
    raw_post.insert(0, "target_node", target_name)
    raw_stem.insert(0, "DCC", dcc)
    raw_post.insert(0, "DCC", dcc)
    raw_out = pd.concat([raw_stem, raw_post], ignore_index=True)

    safe_dcc = str(dcc).replace("/", "_").replace(" ", "_")
    raw_file = output_dir / f"{safe_dcc}.raw_stem_post_events.csv"
    raw_out.to_csv(raw_file, index=False, na_rep="NA")
    raw_files.append(raw_file)

    print(
        f"[INFO] {dcc} ({target_name}): "
        f"stem_nodes={len(stem_nodes)}, post_nodes={len(post_nodes)}, "
        f"stem_events={stem_df.shape[0]}, post_events={post_df.shape[0]}"
    )

    stage_rows.append(count_stage(stem_df, dcc, target_name, "pre"))
    stage_rows.append(count_stage(post_df, dcc, target_name, "post"))


if len(stage_rows) == 0:
    sys.exit("Error: no valid target node was processed")

gene_stage = pd.concat(stage_rows, ignore_index=True)


# =========================================================
# Make wide stem/post summary per DCC and gene
# =========================================================
print("[INFO] making stem/post summary...")

if gene_stage.shape[0] == 0:
    sys.exit("Error: no events found for target nodes")

wide = gene_stage.pivot_table(
    index=["DCC", "target_node", "gene"],
    columns="stage",
    values=[
        "rSNP", "mSNP", "r_syn", "r_nonsyn", "m_syn", "m_nonsyn",
        "synonymous", "nonsynonymous", "total_events",
        "rm_raw", "NS_S_raw"
    ],
    fill_value=0,
    aggfunc="sum"
)

wide.columns = [f"{metric}_{stage}" for metric, stage in wide.columns]
wide = wide.reset_index()

for col in [
    "rSNP_pre", "mSNP_pre", "r_syn_pre", "r_nonsyn_pre",
    "m_syn_pre", "m_nonsyn_pre", "synonymous_pre",
    "nonsynonymous_pre", "total_events_pre",
    "rSNP_post", "mSNP_post", "r_syn_post", "r_nonsyn_post",
    "m_syn_post", "m_nonsyn_post", "synonymous_post",
    "nonsynonymous_post", "total_events_post"
]:
    if col not in wide.columns:
        wide[col] = 0

for col in ["rm_raw_pre", "rm_raw_post", "NS_S_raw_pre", "NS_S_raw_post"]:
    if col not in wide.columns:
        wide[col] = np.nan

# Recalculate ratios after pivot so missing stages are handled consistently.
wide["rm_raw_pre"] = np.where(
    wide["mSNP_pre"] > 0,
    wide["rSNP_pre"] / wide["mSNP_pre"],
    np.nan
)
wide["rm_raw_post"] = np.where(
    wide["mSNP_post"] > 0,
    wide["rSNP_post"] / wide["mSNP_post"],
    np.nan
)
wide["NS_S_raw_pre"] = np.where(
    wide["synonymous_pre"] > 0,
    wide["nonsynonymous_pre"] / wide["synonymous_pre"],
    np.nan
)
wide["NS_S_raw_post"] = np.where(
    wide["synonymous_post"] > 0,
    wide["nonsynonymous_post"] / wide["synonymous_post"],
    np.nan
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
    by=["DCC", "all_total_events", "total_events_pre", "total_events_post"],
    ascending=[True, False, False, False]
)


# =========================================================
# Output
# =========================================================
dcc_files = []

# Per-DCC clean gene summary files.
for dcc, sub in wide.groupby("DCC"):
    safe_dcc = str(dcc).replace("/", "_").replace(" ", "_")
    out_file = output_dir / f"{safe_dcc}.gene_stem_post_summary.csv"
    sub.to_csv(out_file, index=False, na_rep="NA")
    dcc_files.append(out_file)

# Wide merged table: one gene per row, preserving each DCC's summary columns.
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
    dcc_part = dcc_part.reset_index(drop=True)
    wide_merged_parts.append(dcc_part)

wide_merged = pd.concat(wide_merged_parts, axis=1)
wide_merged_file = output_dir / "all_DCC_gene_stem_post_summary_wide.csv"
wide_merged.to_csv(wide_merged_file, index=False, na_rep="NA")

# Total DCC1-7 counts per gene; ratios are recalculated after summing counts.
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
    wide.groupby("gene")["DCC"].nunique().reindex(total["gene"]).to_numpy()
)
total["rm_raw_pre"] = np.where(
    total["mSNP_pre"] > 0,
    total["rSNP_pre"] / total["mSNP_pre"],
    np.nan
)
total["rm_raw_post"] = np.where(
    total["mSNP_post"] > 0,
    total["rSNP_post"] / total["mSNP_post"],
    np.nan
)
total["NS_S_raw_pre"] = np.where(
    total["synonymous_pre"] > 0,
    total["nonsynonymous_pre"] / total["synonymous_pre"],
    np.nan
)
total["NS_S_raw_post"] = np.where(
    total["synonymous_post"] > 0,
    total["nonsynonymous_post"] / total["synonymous_post"],
    np.nan
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
    ascending=False
)
total_file = output_dir / "DCC1-7_total.gene_stem_post_summary.csv"
total.to_csv(total_file, index=False, na_rep="NA")

print("\n==============================")
print("Finished")
print("==============================")
print("Output:")
for x in raw_files:
    print(x)
for x in dcc_files:
    print(x)
print(wide_merged_file)
print(total_file)
