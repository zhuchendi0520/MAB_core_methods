#!/usr/bin/env python3
import argparse
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd


DCC_ORDER = ["DCC1", "DCC2", "DCC4", "DCC5", "DCC3", "DCC6", "DCC7"]
ABS_DCC = {"DCC1", "DCC2", "DCC4", "DCC5"}
MAS_DCC = {"DCC3", "DCC6", "DCC7"}


def read_csv(path):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def num(x):
    return pd.to_numeric(x, errors="coerce")


def clean_gene(x):
    if pd.isna(x):
        return ""
    x = str(x).strip()
    if x.lower() in {"nan", "none", ""}:
        return ""
    return x


def pns_from_counts(pns, obs_syn, obs_nsy, exp_syn, exp_nsy, total):
    pns = pd.to_numeric(pd.Series([pns]), errors="coerce").iloc[0]
    obs_syn = pd.to_numeric(pd.Series([obs_syn]), errors="coerce").iloc[0]
    obs_nsy = pd.to_numeric(pd.Series([obs_nsy]), errors="coerce").iloc[0]
    exp_syn = pd.to_numeric(pd.Series([exp_syn]), errors="coerce").iloc[0]
    exp_nsy = pd.to_numeric(pd.Series([exp_nsy]), errors="coerce").iloc[0]
    total = pd.to_numeric(pd.Series([total]), errors="coerce").iloc[0]

    if pd.isna(total):
        total = 0
    if total <= 0:
        return 0.01
    if not pd.isna(pns) and np.isfinite(pns):
        return float(pns)
    if obs_nsy > 0 and exp_syn > 0 and exp_nsy > 0:
        syn_for_pns = obs_syn if obs_syn > 0 else 1
        return float((obs_nsy / exp_nsy) / (syn_for_pns / exp_syn))
    return 0.01


def call_switch(pre, post):
    if pre < 1 and post > 1:
        return "negative->positive"
    if pre > 1 and post < 1:
        return "positive->negative"
    return "unchanged"


def build_abs_map(path):
    df = read_csv(path)
    df = df.sort_values(["ref2_locus", "bitscore", "pident", "qcov", "scov"], ascending=[True, False, False, False, False])
    out = defaultdict(list)
    for _, row in df.iterrows():
        mab = clean_gene(row.get("ref2_locus"))
        for key in ["ref1_locus", "ref1_gene", "ref1_old_locus"]:
            val = clean_gene(row.get(key))
            if mab and val and val not in out[mab]:
                out[mab].append(val)
    return out


def build_mas_map(path):
    df = read_csv(path)
    df = df.sort_values(["ref1_locus", "bitscore", "pident", "qcov", "scov"], ascending=[True, False, False, False, False])
    out = defaultdict(list)
    for _, row in df.iterrows():
        mab = clean_gene(row.get("ref1_locus"))
        for key in ["ref2_locus", "ref2_gene", "ref2_old_locus"]:
            val = clean_gene(row.get(key))
            if mab and val and val not in out[mab]:
                out[mab].append(val)
    return out


def load_post_file(path):
    df = read_csv(path)
    df["GENE"] = df["GENE"].map(clean_gene)
    df["_tokens"] = df["GENE"].map(lambda x: re.findall(r"[A-Za-z0-9]+_RS\d+|[A-Za-z][A-Za-z0-9_]+", x))
    rows = []
    for _, row in df.iterrows():
        tokens = row["_tokens"] or [row["GENE"]]
        for token in tokens:
            rows.append({
                "token": token,
                "pNS": row.get("pNS"),
                "OBSERVED_SYN": row.get("OBSERVED_SYN"),
                "OBSERVED_NSY": row.get("OBSERVED_NSY"),
                "EXPECTED_SYN": row.get("EXPECTED_SYN"),
                "EXPECTED_NSY": row.get("EXPECTED_NSY"),
                "TOTAL": row.get("TOTAL"),
            })
    return pd.DataFrame(rows)


def aggregate_post_pns(post_df, candidate_tokens):
    hit = post_df[post_df["token"].isin(candidate_tokens)].copy()
    if hit.empty:
        return 0.01, 0, 0, 0, "unmapped_or_absent"

    obs_syn = num(hit["OBSERVED_SYN"]).fillna(0).sum()
    obs_nsy = num(hit["OBSERVED_NSY"]).fillna(0).sum()
    exp_syn = num(hit["EXPECTED_SYN"]).fillna(0).sum()
    exp_nsy = num(hit["EXPECTED_NSY"]).fillna(0).sum()
    total = num(hit["TOTAL"]).fillna(0).sum()
    pns = pns_from_counts(np.nan, obs_syn, obs_nsy, exp_syn, exp_nsy, total)
    return pns, int(total), int(obs_syn), int(obs_nsy), "mapped"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pns-dir", required=True)
    ap.add_argument("--abs-map", required=True, help="EFV83/ref1 to MAB/ref2 mapping CSV")
    ap.add_argument("--mas-map", required=True, help="MAB/ref1 to MMASJCM/ref2 mapping CSV")
    ap.add_argument("--total-switch", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threshold", type=float, default=0.9)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    abs_map = build_abs_map(args.abs_map)
    mas_map = build_mas_map(args.mas_map)

    total = read_csv(args.total_switch)
    total["posterior_switch_probability"] = num(total["posterior_switch_probability"])
    target = total[total["posterior_switch_probability"] >= args.threshold].copy()
    target["locus"] = target["locus"].map(clean_gene)
    target = target[target["locus"] != ""].drop_duplicates("locus")

    post_files = {
        "DCC1": "Mab_DCC1_gene.csv",
        "DCC2": "Mab_DCC2_gene.csv",
        "DCC4": "Mab_DCC4_gene.csv",
        "DCC5": "Mab_DCC5_gene.csv",
        "DCC3": "Mas_DCC3_gene.csv",
        "DCC6": "Mas_DCC6_gene.csv",
        "DCC7": "Mas_DCC7_gene.csv",
    }
    post_tables = {dcc: load_post_file(os.path.join(args.pns_dir, fn)) for dcc, fn in post_files.items()}

    records = []
    for dcc in DCC_ORDER:
        pre = read_csv(os.path.join(args.pns_dir, f"{dcc}.pNS_by_gene.csv"))
        pre = pre[pre["stage"].astype(str).str.lower() == "pre"].copy()
        pre["gene"] = pre["gene"].map(clean_gene)
        pre_by_gene = pre.set_index("gene", drop=False)

        for _, gene_row in target.iterrows():
            gene = gene_row["locus"]
            if gene in pre_by_gene.index:
                pr = pre_by_gene.loc[gene]
                if isinstance(pr, pd.DataFrame):
                    obs_syn_pre = num(pr["observed_syn_raw"]).fillna(0).sum()
                    obs_nsy_pre = num(pr["observed_nonsyn"]).fillna(0).sum()
                    exp_syn_pre = num(pr["expected_syn"]).fillna(0).sum()
                    exp_nsy_pre = num(pr["expected_nonsyn"]).fillna(0).sum()
                    total_pre = num(pr["total_mut_raw"]).fillna(0).sum()
                    pns_pre = pns_from_counts(np.nan, obs_syn_pre, obs_nsy_pre, exp_syn_pre, exp_nsy_pre, total_pre)
                else:
                    total_pre = pd.to_numeric(pr.get("total_mut_raw"), errors="coerce")
                    pns_pre = pns_from_counts(pr.get("pNS"), pr.get("observed_syn_raw"), pr.get("observed_nonsyn"), pr.get("expected_syn"), pr.get("expected_nonsyn"), total_pre)
                    obs_syn_pre = pd.to_numeric(pr.get("observed_syn_raw"), errors="coerce")
                    obs_nsy_pre = pd.to_numeric(pr.get("observed_nonsyn"), errors="coerce")
            else:
                pns_pre, total_pre, obs_syn_pre, obs_nsy_pre = 0.01, 0, 0, 0

            if dcc in ABS_DCC:
                candidates = abs_map.get(gene, [])
            else:
                candidates = mas_map.get(gene, [])
            pns_post, total_post, obs_syn_post, obs_nsy_post, map_status = aggregate_post_pns(post_tables[dcc], candidates)

            records.append({
                "gene": gene,
                "Preferred_name": gene_row.get("gene", ""),
                "DCC": dcc,
                "total_switch_direction": gene_row.get("posterior_switch_direction", ""),
                "pNS_pre": pns_pre,
                "pNS_post_full": pns_post,
                "pre_total_mut": int(0 if pd.isna(total_pre) else total_pre),
                "post_total_mut": total_post,
                "pre_observed_syn": int(0 if pd.isna(obs_syn_pre) else obs_syn_pre),
                "pre_observed_nonsyn": int(0 if pd.isna(obs_nsy_pre) else obs_nsy_pre),
                "post_observed_syn": obs_syn_post,
                "post_observed_nonsyn": obs_nsy_post,
                "post_mapping_status": map_status,
                "switch_call": call_switch(float(pns_pre), float(pns_post)),
            })

    long = pd.DataFrame(records)
    long["switch_numeric"] = long["switch_call"].map({
        "negative->positive": 1,
        "positive->negative": -1,
        "unchanged": 0,
    })
    long["dcc_order"] = long["DCC"].map({d: i for i, d in enumerate(DCC_ORDER)})
    gene_order = (
        long.groupby("gene")["switch_numeric"]
        .apply(lambda s: int((s != 0).sum()))
        .sort_values(ascending=False)
        .index.tolist()
    )
    long["gene"] = pd.Categorical(long["gene"], categories=gene_order, ordered=True)
    long = long.sort_values(["gene", "dcc_order"]).drop(columns=["dcc_order"])

    matrix = long.pivot(index="gene", columns="DCC", values="switch_call").reset_index()
    for dcc in DCC_ORDER:
        if dcc not in matrix.columns:
            matrix[dcc] = "unchanged"
    matrix = matrix[["gene"] + DCC_ORDER]

    summary = (
        long.groupby(["DCC", "switch_call"], as_index=False)
        .size()
        .pivot(index="DCC", columns="switch_call", values="size")
        .fillna(0)
        .astype(int)
        .reset_index()
    )
    summary["mapped_post_genes"] = long.groupby("DCC")["post_mapping_status"].apply(lambda s: int((s == "mapped").sum())).values
    summary["target_genes"] = len(target)

    long.to_csv(os.path.join(args.outdir, "DCC_pNS_switch_convergence_long_fullpost_simple.csv"), index=False)
    matrix.to_csv(os.path.join(args.outdir, "DCC_pNS_switch_convergence_matrix_fullpost_simple.csv"), index=False)
    summary.to_csv(os.path.join(args.outdir, "DCC_pNS_switch_convergence_summary_fullpost_simple.csv"), index=False)


if __name__ == "__main__":
    main()
