#!/usr/bin/env python3
import argparse
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd


DCC_ORDER = ["DCC1", "DCC2", "DCC4", "DCC5", "DCC3", "DCC6", "DCC7"]
ABS_DCC = {"DCC1", "DCC2", "DCC4", "DCC5"}


def clean(x):
    if pd.isna(x):
        return ""
    x = str(x).strip()
    if x.lower() in {"nan", "none", ""}:
        return ""
    return x


def nseries(x):
    return pd.to_numeric(x, errors="coerce")


def pns_from_counts(pns, obs_syn, obs_nsy, exp_syn, exp_nsy, total):
    pns = pd.to_numeric(pd.Series([pns]), errors="coerce").iloc[0]
    obs_syn = pd.to_numeric(pd.Series([obs_syn]), errors="coerce").iloc[0]
    obs_nsy = pd.to_numeric(pd.Series([obs_nsy]), errors="coerce").iloc[0]
    exp_syn = pd.to_numeric(pd.Series([exp_syn]), errors="coerce").iloc[0]
    exp_nsy = pd.to_numeric(pd.Series([exp_nsy]), errors="coerce").iloc[0]
    total = pd.to_numeric(pd.Series([total]), errors="coerce").iloc[0]

    if pd.isna(total) or total <= 0:
        return 0.01
    if not pd.isna(pns) and np.isfinite(pns):
        return float(pns)
    if obs_nsy > 0 and exp_syn > 0 and exp_nsy > 0:
        syn_for_pns = obs_syn if obs_syn > 0 else 1
        return float((obs_nsy / exp_nsy) / (syn_for_pns / exp_syn))
    return 0.01


def call_switch(pre, post):
    if pre < 1 and post > 1:
        return "purifying_to_positive"
    if pre > 1 and post < 1:
        return "positive_to_purifying"
    return "unchanged"


def read_core_table(path):
    df = pd.read_csv(path, sep="\t")
    df.columns = ["gene", "HMM_call"]
    df["gene"] = df["gene"].map(clean)
    df["HMM_call"] = df["HMM_call"].map(clean)
    return df[df["gene"] != ""].drop_duplicates("gene")


def build_mapping(path, from_col, to_cols):
    df = pd.read_csv(path)
    sort_cols = [from_col]
    ascending = [True]
    for col in ["bitscore", "pident", "qcov", "scov"]:
        if col in df.columns:
            sort_cols.append(col)
            ascending.append(False)
    df = df.sort_values(sort_cols, ascending=ascending)
    out = defaultdict(list)
    for _, row in df.iterrows():
        src = clean(row.get(from_col))
        for to_col in to_cols:
            dst = clean(row.get(to_col))
            if src and dst and dst not in out[src]:
                out[src].append(dst)
    return out


def load_pre_by_gene(path):
    df = pd.read_csv(path)
    df = df[df["stage"].astype(str).str.lower() == "pre"].copy()
    df["gene"] = df["gene"].map(clean)
    rows = {}
    for gene, sub in df.groupby("gene"):
        obs_syn = nseries(sub["observed_syn_raw"]).fillna(0).sum()
        obs_nsy = nseries(sub["observed_nonsyn"]).fillna(0).sum()
        exp_syn = nseries(sub["expected_syn"]).fillna(0).sum()
        exp_nsy = nseries(sub["expected_nonsyn"]).fillna(0).sum()
        total = nseries(sub["total_mut_raw"]).fillna(0).sum()
        rows[gene] = {
            "pNS": pns_from_counts(np.nan, obs_syn, obs_nsy, exp_syn, exp_nsy, total),
            "total": int(total),
            "obs_syn": int(obs_syn),
            "obs_nsy": int(obs_nsy),
        }
    return rows


def load_post_by_token(path):
    df = pd.read_csv(path)
    df["GENE"] = df["GENE"].map(clean)
    expanded = []
    for _, row in df.iterrows():
        tokens = re.findall(r"[A-Za-z0-9]+_RS\d+|[A-Za-z][A-Za-z0-9_]+", row["GENE"])
        if not tokens:
            tokens = [row["GENE"]]
        for token in tokens:
            expanded.append({
                "token": token,
                "obs_syn": row.get("OBSERVED_SYN"),
                "obs_nsy": row.get("OBSERVED_NSY"),
                "exp_syn": row.get("EXPECTED_SYN"),
                "exp_nsy": row.get("EXPECTED_NSY"),
                "total": row.get("TOTAL"),
            })
    exp = pd.DataFrame(expanded)
    rows = {}
    for token, sub in exp.groupby("token"):
        obs_syn = nseries(sub["obs_syn"]).fillna(0).sum()
        obs_nsy = nseries(sub["obs_nsy"]).fillna(0).sum()
        exp_syn = nseries(sub["exp_syn"]).fillna(0).sum()
        exp_nsy = nseries(sub["exp_nsy"]).fillna(0).sum()
        total = nseries(sub["total"]).fillna(0).sum()
        rows[token] = {
            "obs_syn": obs_syn,
            "obs_nsy": obs_nsy,
            "exp_syn": exp_syn,
            "exp_nsy": exp_nsy,
            "total": total,
        }
    return rows


def aggregate_post(candidates, post_by_token):
    obs_syn = obs_nsy = exp_syn = exp_nsy = total = 0.0
    mapped = False
    for token in candidates:
        hit = post_by_token.get(token)
        if hit is None:
            continue
        mapped = True
        obs_syn += hit["obs_syn"]
        obs_nsy += hit["obs_nsy"]
        exp_syn += hit["exp_syn"]
        exp_nsy += hit["exp_nsy"]
        total += hit["total"]
    return {
        "pNS": pns_from_counts(np.nan, obs_syn, obs_nsy, exp_syn, exp_nsy, total),
        "total": int(total),
        "obs_syn": int(obs_syn),
        "obs_nsy": int(obs_nsy),
        "mapping_available": len(candidates) > 0,
        "event_observed": total > 0,
    }


def classify_gene(states):
    changed = [s for s in states if s != "unchanged"]
    if len(changed) == 0:
        return "all_unchanged", 0, ""
    unique_changed = sorted(set(changed))
    if len(changed) == 7 and len(unique_changed) == 1:
        return f"all_7_{unique_changed[0]}", 7, unique_changed[0]
    if len(unique_changed) == 1:
        return f"changed_DCCs_same_direction_n{len(changed)}", len(changed), unique_changed[0]
    return "mixed_direction", len(changed), "mixed"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pns-dir", required=True)
    ap.add_argument("--core-table", required=True)
    ap.add_argument("--abs-map", required=True)
    ap.add_argument("--mas-map", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    core = read_core_table(args.core_table)
    abs_map = build_mapping(args.abs_map, "ref2_locus", ["ref1_locus", "ref1_gene", "ref1_old_locus"])
    mas_map = build_mapping(args.mas_map, "ref1_locus", ["ref2_locus", "ref2_gene", "ref2_old_locus"])

    pre_tables = {dcc: load_pre_by_gene(os.path.join(args.pns_dir, f"{dcc}.pNS_by_gene.csv")) for dcc in DCC_ORDER}
    post_file = {
        "DCC1": "Mab_DCC1_gene.csv",
        "DCC2": "Mab_DCC2_gene.csv",
        "DCC4": "Mab_DCC4_gene.csv",
        "DCC5": "Mab_DCC5_gene.csv",
        "DCC3": "Mas_DCC3_gene.csv",
        "DCC6": "Mas_DCC6_gene.csv",
        "DCC7": "Mas_DCC7_gene.csv",
    }
    post_tables = {dcc: load_post_by_token(os.path.join(args.pns_dir, fn)) for dcc, fn in post_file.items()}

    long_rows = []
    summary_rows = []
    for _, row in core.iterrows():
        gene = row["gene"]
        states = []
        compact = {}
        for dcc in DCC_ORDER:
            pre = pre_tables[dcc].get(gene, {"pNS": 0.01, "total": 0, "obs_syn": 0, "obs_nsy": 0})
            candidates = abs_map.get(gene, []) if dcc in ABS_DCC else mas_map.get(gene, [])
            post = aggregate_post(candidates, post_tables[dcc])
            state = call_switch(pre["pNS"], post["pNS"])
            states.append(state)
            compact[dcc] = state
            long_rows.append({
                "gene": gene,
                "HMM_call": row["HMM_call"],
                "DCC": dcc,
                "pNS_pre": pre["pNS"],
                "pNS_post_full": post["pNS"],
                "pre_total_mut": pre["total"],
                "post_total_mut": post["total"],
                "pre_observed_syn": pre["obs_syn"],
                "pre_observed_nonsyn": pre["obs_nsy"],
                "post_observed_syn": post["obs_syn"],
                "post_observed_nonsyn": post["obs_nsy"],
                "post_reference_mapping_available": post["mapping_available"],
                "post_event_observed": post["event_observed"],
                "switch_call": state,
            })
        consistency_class, n_changed, shared_direction = classify_gene(states)
        summary_rows.append({
            "gene": gene,
            "HMM_call": row["HMM_call"],
            "n_changed_DCC": n_changed,
            "shared_changed_direction": shared_direction,
            "consistency_class": consistency_class,
            **compact,
        })

    long = pd.DataFrame(long_rows)
    summary = pd.DataFrame(summary_rows)

    class_counts = summary["consistency_class"].value_counts().rename_axis("consistency_class").reset_index(name="gene_count")
    n_direction = (
        summary.groupby(["n_changed_DCC", "shared_changed_direction"], dropna=False)
        .size()
        .reset_index(name="gene_count")
        .sort_values(["n_changed_DCC", "shared_changed_direction"], ascending=[False, True])
    )

    long.to_csv(os.path.join(args.outdir, "all_core_DCC_pNS_switch_long_fullpost.csv"), index=False)
    summary.to_csv(os.path.join(args.outdir, "all_core_DCC_pNS_switch_direction_consistency.csv"), index=False)
    class_counts.to_csv(os.path.join(args.outdir, "all_core_DCC_pNS_switch_consistency_class_counts.csv"), index=False)
    n_direction.to_csv(os.path.join(args.outdir, "all_core_DCC_pNS_switch_n_changed_direction_counts.csv"), index=False)

    print("core_genes", len(core))
    print(class_counts.to_string(index=False))
    print("\nBy n_changed and direction:")
    print(n_direction.to_string(index=False))


if __name__ == "__main__":
    main()
