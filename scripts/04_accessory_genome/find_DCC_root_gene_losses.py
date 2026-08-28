#!/usr/bin/env python3
"""Find accessory-gene losses on DCC stem branches.

A candidate loss must be rare inside a monophyletic DCC and common both in its
immediate sister lineage and across all samples outside that DCC. The script
reports strict (<=5% vs >=95%) and relaxed (<=30% vs >=70%) calls. DCC7's
known off-clade label is automatically removed by selecting the largest pure
monophyletic clade for each DCC.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo


DCCS = [f"DCC{i}" for i in range(1, 8)]


def parent_map(tree):
    return {child: parent for parent in tree.find_clades(order="level") for child in parent.clades}


def largest_pure_dcc_clade(tree, group_map, dcc):
    """Return largest clade whose every descendant is explicitly labelled dcc."""
    best = None
    best_key = (-1, -1)
    for clade in tree.find_clades(order="postorder"):
        tips = [x.name for x in clade.get_terminals()]
        if not tips or any(group_map.get(x) != dcc for x in tips):
            continue
        n_dcc = len(tips)
        key = (n_dcc, len(tips))
        if key > best_key:
            best, best_key = clade, key
    if best is None:
        raise ValueError(f"No pure clade found for {dcc}")
    return best


def pct_present(values):
    if len(values) == 0:
        return np.nan
    return 100.0 * np.count_nonzero(values) / len(values)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True)
    ap.add_argument("--groups", required=True)
    ap.add_argument("--rtab", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    groups = pd.read_csv(args.groups, sep="\t", dtype=str)
    if not {"strain", "DCC"}.issubset(groups.columns):
        raise ValueError("Group file must contain strain and DCC columns")
    group_map = dict(zip(groups["strain"], groups["DCC"]))

    tree = Phylo.read(args.tree, "newick")
    pmap = parent_map(tree)
    tree_tips = {x.name for x in tree.get_terminals()}

    with open(args.rtab, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
    matrix_samples = header[1:]
    sample_index = {s: i for i, s in enumerate(matrix_samples)}

    clade_info = {}
    qc_rows = []
    for dcc in DCCS:
        assigned = {s for s, g in group_map.items() if g == dcc and s in tree_tips}
        clade = largest_pure_dcc_clade(tree, group_map, dcc)
        inside_tree = {x.name for x in clade.get_terminals()}
        inside = sorted(inside_tree & set(matrix_samples))
        parent = pmap.get(clade)
        if parent is None:
            raise ValueError(f"Selected clade for {dcc} is the tree root")
        sister_tree = {x.name for x in parent.get_terminals()} - inside_tree
        sister = sorted(sister_tree & set(matrix_samples))
        outside = sorted(set(matrix_samples) - set(inside))
        excluded = sorted(assigned - inside_tree)
        clade_info[dcc] = {
            "inside_idx": np.array([sample_index[x] for x in inside]),
            "sister_idx": np.array([sample_index[x] for x in sister]),
            "outside_idx": np.array([sample_index[x] for x in outside]),
            "n_inside": len(inside),
            "n_sister": len(sister),
        }
        qc_rows.append({
            "DCC": dcc,
            "n_group_assigned_on_tree": len(assigned),
            "n_used_in_monophyletic_clade": len(inside),
            "n_immediate_sister_samples": len(sister),
            "excluded_off_clade_strains": "|".join(excluded),
        })

    pd.DataFrame(qc_rows).to_csv(outdir / "DCC_root_loss_clade_QC.csv", index=False)

    calls = []
    chunks = pd.read_csv(args.rtab, sep="\t", chunksize=500, dtype={"Gene": str})
    for chunk in chunks:
        genes = chunk.iloc[:, 0].astype(str).to_numpy()
        values = chunk.iloc[:, 1:].to_numpy(dtype=np.uint8, copy=False)
        for dcc, info in clade_info.items():
            inside = values[:, info["inside_idx"]]
            sister = values[:, info["sister_idx"]]
            outside = values[:, info["outside_idx"]]
            in_pct = 100 * inside.mean(axis=1)
            sis_pct = 100 * sister.mean(axis=1)
            out_pct = 100 * outside.mean(axis=1)

            strict = (in_pct <= 5) & (sis_pct >= 95) & (out_pct >= 95)
            relaxed = (in_pct <= 30) & (sis_pct >= 70) & (out_pct >= 70)
            selected = relaxed
            for i in np.flatnonzero(selected):
                calls.append({
                    "DCC": dcc,
                    "Gene": genes[i],
                    "threshold_class": "strict_95_5" if strict[i] else "relaxed_only_70_30",
                    "DCC_present_n": int(inside[i].sum()),
                    "DCC_total_n": inside.shape[1],
                    "DCC_carriage_pct": round(float(in_pct[i]), 6),
                    "sister_present_n": int(sister[i].sum()),
                    "sister_total_n": sister.shape[1],
                    "sister_carriage_pct": round(float(sis_pct[i]), 6),
                    "outside_present_n": int(outside[i].sum()),
                    "outside_total_n": outside.shape[1],
                    "outside_carriage_pct": round(float(out_pct[i]), 6),
                })

    result = pd.DataFrame(calls)
    if result.empty:
        result = pd.DataFrame(columns=[
            "DCC", "Gene", "threshold_class", "DCC_present_n", "DCC_total_n",
            "DCC_carriage_pct", "sister_present_n", "sister_total_n",
            "sister_carriage_pct", "outside_present_n", "outside_total_n",
            "outside_carriage_pct",
        ])
    result = result.sort_values(["DCC", "threshold_class", "DCC_carriage_pct", "Gene"])
    result.to_csv(outdir / "all_DCC_root_gene_losses_strict_and_relaxed.csv", index=False)

    for dcc in DCCS:
        sub = result[result["DCC"] == dcc]
        sub[sub["threshold_class"] == "strict_95_5"].to_csv(
            outdir / f"{dcc}_root_gene_losses_strict_95_5.csv", index=False
        )
        sub.to_csv(outdir / f"{dcc}_root_gene_losses_relaxed_70_30_including_strict.csv", index=False)

    summary = (
        result.assign(strict=result["threshold_class"].eq("strict_95_5"))
        .groupby("DCC", as_index=False)
        .agg(
            strict_loss_genes=("strict", "sum"),
            relaxed_total_loss_genes=("Gene", "size"),
        )
    )
    summary["relaxed_only_loss_genes"] = summary["relaxed_total_loss_genes"] - summary["strict_loss_genes"]
    summary.to_csv(outdir / "DCC_root_gene_loss_counts.csv", index=False)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
