#!/usr/bin/env python3
from __future__ import annotations

import os
import argparse
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Rectangle
import textwrap
import numpy as np
import pandas as pd


def clean_num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def classify_rm(row: pd.Series) -> str:
    p = row.get("fisher_p_adj_BH")
    pre = row.get("rm_raw_pre")
    post = row.get("rm_raw_post")
    if pd.isna(p) or pd.isna(pre) or pd.isna(post) or p >= 0.05 or pre == post:
        return "No significant r/m shift"
    if post > pre:
        return "r/m increased"
    return "r/m decreased"


def classify_pns(row: pd.Series) -> str:
    prob = row.get("posterior_switch_probability")
    direction = str(row.get("posterior_switch_direction", ""))
    if pd.isna(prob) or prob < 0.9:
        return "No pN/pS switch"
    if direction == "purifying_to_positive":
        return "Purifying -> Positive selection"
    if direction == "positive_to_purifying":
        return "Positive -> Purifying selection"
    return "No pN/pS switch"


def add_band(ax, x0, x1, y0a, y0b, y1a, y1b, color, alpha=0.72):
    dx = x1 - x0
    verts = [
        (x0, y0a),
        (x0 + dx * 0.45, y0a),
        (x1 - dx * 0.45, y1a),
        (x1, y1a),
        (x1, y1b),
        (x1 - dx * 0.45, y1b),
        (x0 + dx * 0.45, y0b),
        (x0, y0b),
        (x0, y0a),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.LINETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CLOSEPOLY,
    ]
    patch = PathPatch(
        MplPath(verts, codes),
        facecolor=color,
        edgecolor="white",
        linewidth=0.8,
        alpha=alpha,
    )
    ax.add_patch(patch)


def layout_segments(counts: pd.DataFrame, cats: list[str], side: str, gap=0.035, gap_after=None):
    totals = counts.groupby(side)["n"].sum().reindex(cats).fillna(0)
    total_n = float(totals.sum())
    if gap_after is None:
        gap_after = [gap] * (len(cats) - 1)
    available = 1.0 - sum(gap_after)
    heights = totals / total_n * available if total_n else totals
    y_top = 1.0
    bounds = {}
    for i, (cat, h) in enumerate(zip(cats, heights)):
        y_bottom = y_top - float(h)
        bounds[cat] = [y_bottom, y_top]
        if i < len(cats) - 1:
            y_top = y_bottom - gap_after[i]
    return bounds, totals


def wrapped_gene_label(row: pd.Series) -> str:
    gene = str(row["gene"])
    symbol = row.get("gene_pns")
    if isinstance(symbol, str) and symbol.strip() and symbol.strip().lower() != "nan":
        return f"{symbol} ({gene})"
    return gene


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot coupled r/m and pN/pS shifts")
    parser.add_argument("rm_csv", type=Path)
    parser.add_argument("pns_csv", type=Path)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    rm_file = args.rm_csv
    pns_file = args.pns_csv
    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    rm = pd.read_csv(rm_file)
    pns = pd.read_csv(pns_file)

    rm = rm[rm["DCC"].astype(str).eq("DCC1-7_total")].copy()
    pns = pns[pns["DCC"].astype(str).eq("DCC1-7_total")].copy()

    for col in ["rm_raw_pre", "rm_raw_post", "fisher_p_adj_BH"]:
        rm[col] = clean_num(rm[col])
    for col in ["posterior_switch_probability", "pNS_pre", "pNS_post"]:
        pns[col] = clean_num(pns[col])

    rm["rm_state"] = rm.apply(classify_rm, axis=1)
    pns["pns_state"] = pns.apply(classify_pns, axis=1)

    merged = rm.merge(
        pns[[
            "locus",
            "gene",
            "product",
            "pNS_pre",
            "pNS_post",
            "posterior_switch_probability",
            "posterior_switch_direction",
            "pns_state",
        ]],
        left_on="gene",
        right_on="locus",
        how="inner",
        suffixes=("", "_pns"),
    )

    merged_out = out_dir / "coupled_rm_pns_gene_categories.csv"
    merged.to_csv(merged_out, index=False)

    pns_order = [
        "Purifying \u2192 Positive selection",
        "Positive \u2192 Purifying selection",
        "No pN/pS switch",
    ]
    rm_order = ["r/m increased", "r/m decreased", "No significant r/m shift"]

    changed = merged[
        (merged["rm_state"] != "No significant r/m shift")
        | (merged["pns_state"] != "No pN/pS switch")
    ].copy()

    counts = (
        changed.groupby(["rm_state", "pns_state"], dropna=False)
        .size()
        .reset_index(name="n")
    )
    counts_all = (
        merged.groupby(["rm_state", "pns_state"], dropna=False)
        .size()
        .reset_index(name="n")
    )
    counts = counts[counts["rm_state"].isin(rm_order) & counts["pns_state"].isin(pns_order)].copy()
    counts_all = counts_all[counts_all["rm_state"].isin(rm_order) & counts_all["pns_state"].isin(pns_order)].copy()
    counts_out = out_dir / "coupled_rm_pns_transition_counts.csv"
    counts_all.to_csv(counts_out, index=False)

    left_bounds, left_totals = layout_segments(counts, pns_order, "pns_state", gap_after=[0.12, 0.055])
    mid_bounds, mid_totals = layout_segments(counts, rm_order, "rm_state", gap_after=[0.05, 0.05])

    left_cursor = {k: v[0] for k, v in left_bounds.items()}
    mid_cursor = {k: v[0] for k, v in mid_bounds.items()}

    palette = {
        "Purifying \u2192 Positive selection": "#E59673",
        "Positive \u2192 Purifying selection": "#8FA7DB",
        "No pN/pS switch": "#D6D6D6",
    }
    coupled_palette = {
        ("r/m increased", "Purifying \u2192 Positive selection"): "#D77854",
        ("r/m increased", "Positive \u2192 Purifying selection"): "#7F97D2",
        ("r/m decreased", "Purifying \u2192 Positive selection"): "#F0B391",
        ("r/m decreased", "Positive \u2192 Purifying selection"): "#AFC0E8",
    }

    total_n = counts["n"].sum()
    scale = (1.0 - 0.035 * (len(pns_order) - 1)) / total_n

    plt.rcParams.update({
        "font.family": "Arial",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, ax = plt.subplots(figsize=(10.2, 5.6))
    x_left, x_mid, x_cat = 0.13, 0.48, 0.77
    bar_w = 0.028

    for _, row in counts.iterrows():
        left = row["pns_state"]
        mid = row["rm_state"]
        h = row["n"] * scale
        y0a, y0b = left_cursor[left], left_cursor[left] + h
        y1a, y1b = mid_cursor[mid], mid_cursor[mid] + h
        color = coupled_palette.get((mid, left), "#BFBFBF")
        alpha = 0.86 if (mid, left) in coupled_palette else 0.20
        add_band(ax, x_left + bar_w / 2, x_mid - bar_w / 2, y0a, y0b, y1a, y1b, color, alpha=alpha)
        left_cursor[left] += h
        mid_cursor[mid] += h

    for cat in pns_order:
        y0, y1 = left_bounds[cat]
        ax.add_patch(Rectangle((x_left - bar_w / 2, y0), bar_w, y1 - y0,
                               facecolor=palette[cat], edgecolor="black", linewidth=0.7, alpha=0.95))
        label = cat.replace("Purifying \u2192 Positive selection", "Purifying \u2192\nPositive selection")
        label = label.replace("Positive \u2192 Purifying selection", "Positive \u2192\nPurifying selection")
        ax.text(x_left - 0.035, (y0 + y1) / 2, f"{label}\n(n={int(left_totals[cat])})",
                ha="right", va="center", fontsize=9.1)

    for cat in rm_order:
        y0, y1 = mid_bounds[cat]
        ax.add_patch(Rectangle((x_mid - bar_w / 2, y0), bar_w, y1 - y0,
                               facecolor="#F2F2F2", edgecolor="black", linewidth=0.7))
        label = cat.replace("No significant r/m shift", "No significant\nr/m shift")
        ax.text(
            x_mid + 0.055,
            (y0 + y1) / 2,
            f"{label}\n(n={int(mid_totals[cat])})",
            ha="left",
            va="center",
            fontsize=9.1,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.86, pad=1.8),
        )

    coupled = changed[
        (changed["rm_state"] != "No significant r/m shift")
        & (changed["pns_state"] != "No pN/pS switch")
    ].copy()
    coupled["coupled_class"] = list(zip(coupled["rm_state"], coupled["pns_state"]))
    coupled_out = out_dir / "coupled_rm_pns_double_shift_genes.csv"
    coupled.to_csv(coupled_out, index=False)

    cat_order = [
        ("r/m decreased", "Purifying \u2192 Positive selection"),
        ("r/m increased", "Purifying \u2192 Positive selection"),
        ("r/m decreased", "Positive \u2192 Purifying selection"),
        ("r/m increased", "Positive \u2192 Purifying selection"),
    ]
    cat_labels = {
        ("r/m decreased", "Purifying \u2192 Positive selection"): "Purifying \u2192 Positive\nr/m decreased",
        ("r/m increased", "Purifying \u2192 Positive selection"): "Purifying \u2192 Positive\nr/m increased",
        ("r/m decreased", "Positive \u2192 Purifying selection"): "Positive \u2192 Purifying\nr/m decreased",
        ("r/m increased", "Positive \u2192 Purifying selection"): "Positive \u2192 Purifying\nr/m increased",
    }
    cat_gap = 0.022
    cat_h = (0.94 - cat_gap * 3) / 4
    y_top = 0.98
    cat_bounds = {}
    for cat in cat_order:
        y0 = y_top - cat_h
        cat_bounds[cat] = (y0, y_top)
        y_top = y0 - cat_gap

    for cat in cat_order:
        y0, y1 = cat_bounds[cat]
        color = coupled_palette[cat]
        group = coupled[coupled["coupled_class"].apply(lambda x: x == cat)].copy()
        group = group.sort_values(["posterior_switch_probability", "gene"], ascending=[False, True])
        ax.add_patch(Rectangle((x_cat - 0.012, y0), 0.21, y1 - y0,
                               facecolor=color, edgecolor="#BFBFBF", linewidth=0.8, alpha=0.22))
        ax.add_patch(Rectangle((x_cat - 0.012, y0), 0.012, y1 - y0,
                               facecolor=color, edgecolor=color, linewidth=0.5, alpha=0.95))
        ax.text(x_cat + 0.006, y1 - 0.016, f"{cat_labels[cat]} (n={len(group)})",
                ha="left", va="top", fontsize=8.0, fontweight="bold", color="black")
        genes = [wrapped_gene_label(r) for _, r in group.iterrows()]
        gene_text = "\n".join(textwrap.shorten(g, width=32, placeholder="...") for g in genes) if genes else "None"
        ax.text(x_cat + 0.006, y1 - 0.078, gene_text,
                ha="left", va="top", fontsize=7.5, color="black", linespacing=1.25)

    # Draw emphasized connectors from r/m bars into the corresponding candidate boxes.
    for cat in cat_order:
        rm_state, pns_state = cat
        n = int(((coupled["rm_state"] == rm_state) & (coupled["pns_state"] == pns_state)).sum())
        if n == 0:
            continue
        y_mid_source = np.mean(mid_bounds[rm_state])
        y0, y1 = cat_bounds[cat]
        y_mid_target = (y0 + y1) / 2
        h = max(0.008, n * scale * 2.3)
        add_band(
            ax,
            x_mid + bar_w / 2,
            x_cat - 0.012,
            y_mid_source - h / 2,
            y_mid_source + h / 2,
            y_mid_target - h / 2,
            y_mid_target + h / 2,
            coupled_palette[cat],
            alpha=0.72,
        )

    ax.text(x_left, 1.08, "pN/pS switch", ha="center", va="bottom", fontsize=12.5, fontweight="bold")
    ax.text(x_mid, 1.08, "r/m shift", ha="center", va="bottom", fontsize=12.5, fontweight="bold")
    ax.text(x_cat + 0.075, 1.08, "Coupled switch genes", ha="center", va="bottom", fontsize=12.5, fontweight="bold")

    legend_handles = [
        Rectangle((0, 0), 1, 1, facecolor="#E59673", edgecolor="black", linewidth=0.6, label="Purifying \u2192 Positive selection"),
        Rectangle((0, 0), 1, 1, facecolor="#8FA7DB", edgecolor="black", linewidth=0.6, label="Positive \u2192 Purifying selection"),
        Rectangle((0, 0), 1, 1, facecolor="#D6D6D6", edgecolor="black", linewidth=0.6, label="Single-axis shift"),
    ]
    ax.legend(
        handles=legend_handles,
        frameon=False,
        ncol=3,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.13),
        fontsize=9.5,
        handlelength=1.1,
        columnspacing=1.1,
    )

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-0.02, 1.12)
    ax.axis("off")
    fig.tight_layout()

    png = out_dir / "coupled_rm_pns_alluvial.png"
    pdf = out_dir / "coupled_rm_pns_alluvial.pdf"
    fig.savefig(png, dpi=350, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    print(merged_out)
    print(counts_out)
    print(coupled_out)
    print(png)
    print(pdf)


if __name__ == "__main__":
    main()
