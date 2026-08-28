#!/usr/bin/env python3
from __future__ import annotations

import math
import os
import argparse
from pathlib import Path
import textwrap

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, Wedge, Patch
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch
import numpy as np
import pandas as pd


parser = argparse.ArgumentParser(description="Plot circular coupling of r/m and pN/pS shifts")
parser.add_argument("category_csv", type=Path)
parser.add_argument("rm_csv", type=Path)
parser.add_argument("pns_csv", type=Path)
parser.add_argument("output_dir", type=Path)
cli = parser.parse_args()
IN_FILE = cli.category_csv
RM_FILE = cli.rm_csv
PNS_FILE = cli.pns_csv
OUT_DIR = cli.output_dir


PNS_ORDER = [
    "Purifying -> Positive selection",
    "Positive -> Purifying selection",
    "No pN/pS switch",
]
RM_ORDER = ["r/m increased", "r/m decreased", "No significant r/m shift"]

ANGLES = {
    "Purifying -> Positive selection": 140,
    "Positive -> Purifying selection": 220,
    "No pN/pS switch": 180,
    "r/m increased": 40,
    "r/m decreased": 320,
    "No significant r/m shift": 0,
}

COLORS = {
    ("r/m increased", "Purifying -> Positive selection"): "#D77854",
    ("r/m decreased", "Purifying -> Positive selection"): "#F0B391",
    ("r/m increased", "Positive -> Purifying selection"): "#7F97D2",
    ("r/m decreased", "Positive -> Purifying selection"): "#AFC0E8",
    "single": "#CFCFCF",
    "pns_pos": "#E59673",
    "pns_neg": "#8FA7DB",
    "rm": "#F2F2F2",
}


def polar_xy(angle_deg: float, radius: float = 1.0) -> tuple[float, float]:
    a = math.radians(angle_deg)
    return radius * math.cos(a), radius * math.sin(a)


def draw_curve(ax, a0, a1, lw, color, alpha, zorder=1):
    x0, y0 = polar_xy(a0, 0.93)
    x1, y1 = polar_xy(a1, 0.93)
    c0 = polar_xy(a0, 0.18)
    c1 = polar_xy(a1, 0.18)
    path = MplPath(
        [(x0, y0), c0, c1, (x1, y1)],
        [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4],
    )
    ax.add_patch(
        PathPatch(
            path,
            facecolor="none",
            edgecolor=color,
            lw=lw,
            alpha=alpha,
            capstyle="round",
            zorder=zorder,
        )
    )


def arc_points(a0: float, a1: float, radius: float, n: int = 12) -> list[tuple[float, float]]:
    return [polar_xy(a, radius) for a in np.linspace(a0, a1, n)]


def draw_ribbon(ax, source_arc, target_arc, color, alpha, zorder=1):
    s0, s1 = source_arc
    t0, t1 = target_arc
    r = 0.925
    source_pts = arc_points(s0, s1, r, 10)
    target_pts = arc_points(t1, t0, r, 10)
    c_s1 = polar_xy(s1, 0.18)
    c_t1 = polar_xy(t1, 0.18)
    c_t0 = polar_xy(t0, 0.18)
    c_s0 = polar_xy(s0, 0.18)

    verts = [source_pts[0]]
    codes = [MplPath.MOVETO]
    for pt in source_pts[1:]:
        verts.append(pt)
        codes.append(MplPath.LINETO)
    verts.extend([c_s1, c_t1, target_pts[0]])
    codes.extend([MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])
    for pt in target_pts[1:]:
        verts.append(pt)
        codes.append(MplPath.LINETO)
    verts.extend([c_t0, c_s0, source_pts[0]])
    codes.extend([MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])
    verts.append(source_pts[0])
    codes.append(MplPath.CLOSEPOLY)

    ax.add_patch(
        PathPatch(
            MplPath(verts, codes),
            facecolor=color,
            edgecolor="white",
            linewidth=0.55,
            alpha=alpha,
            zorder=zorder,
        )
    )


def node_arc_bounds(node: str, span: float = 34.0) -> tuple[float, float]:
    center = ANGLES[node]
    return center - span / 2, center + span / 2


def allocate_subarcs(counts: pd.DataFrame, nodes: list[str], node_col: str, other_order: list[str], other_col: str):
    allocations = {}
    for node in nodes:
        a0, a1 = node_arc_bounds(node)
        total = int(counts.loc[counts[node_col].eq(node), "n"].sum())
        cursor = a0
        if total <= 0:
            continue
        for other in other_order:
            row = counts[(counts[node_col] == node) & (counts[other_col] == other)]
            if row.empty:
                continue
            n = int(row["n"].iloc[0])
            width = (a1 - a0) * n / total
            allocations[(node, other)] = (cursor, cursor + width)
            cursor += width
    return allocations


def clean_label(label: str) -> str:
    return (
        label.replace("Purifying -> Positive selection", "Purifying ->\nPositive selection")
        .replace("Positive -> Purifying selection", "Positive ->\nPurifying selection")
        .replace("No significant r/m shift", "No significant\nr/m shift")
        .replace("No pN/pS switch", "No pN/pS\nswitch")
    )


def gene_label(row: pd.Series) -> str:
    symbol = row.get("gene_pns")
    gene = str(row["gene"])
    if isinstance(symbol, str) and symbol.strip() and symbol.strip().lower() != "nan":
        return f"{symbol} ({gene})"
    return gene


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


def load_or_build_gene_categories() -> pd.DataFrame:
    if IN_FILE.exists():
        return pd.read_csv(IN_FILE)

    rm = pd.read_csv(RM_FILE)
    pns = pd.read_csv(PNS_FILE)
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
    merged.to_csv(IN_FILE, index=False)
    return merged


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = load_or_build_gene_categories()

    all_genes = df.copy()
    changed = df[
        (df["rm_state"] != "No significant r/m shift")
        | (df["pns_state"] != "No pN/pS switch")
    ].copy()

    counts = (
        all_genes.groupby(["pns_state", "rm_state"])
        .size()
        .reset_index(name="n")
        .sort_values("n", ascending=False)
    )
    counts.to_csv(OUT_DIR / "coupled_rm_pns_circular_counts.csv", index=False)

    dual = changed[
        (changed["rm_state"] != "No significant r/m shift")
        & (changed["pns_state"] != "No pN/pS switch")
    ].copy()
    dual.to_csv(OUT_DIR / "coupled_rm_pns_circular_double_shift_genes.csv", index=False)

    node_totals = {}
    for node in PNS_ORDER:
        node_totals[node] = int(counts.loc[counts["pns_state"].eq(node), "n"].sum())
    for node in RM_ORDER:
        node_totals[node] = int(counts.loc[counts["rm_state"].eq(node), "n"].sum())

    plt.rcParams.update(
        {
            "font.family": "Arial",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )

    fig = plt.figure(figsize=(8.8, 8.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.05, 0.58], hspace=0.03)
    ax = fig.add_subplot(gs[0])
    ax.set_aspect("equal")
    ax.axis("off")

    pns_arcs = allocate_subarcs(
        counts=counts,
        nodes=PNS_ORDER,
        node_col="pns_state",
        other_order=RM_ORDER,
        other_col="rm_state",
    )
    rm_arcs = allocate_subarcs(
        counts=counts,
        nodes=RM_ORDER,
        node_col="rm_state",
        other_order=PNS_ORDER,
        other_col="pns_state",
    )

    # Draw single-axis ribbons first, then dual-change ribbons on top.
    flow_rows = []
    for _, row in counts.iterrows():
        pns = row["pns_state"]
        rm = row["rm_state"]
        is_dual = rm != "No significant r/m shift" and pns != "No pN/pS switch"
        flow_rows.append((is_dual, pns, rm))

    for is_dual, pns, rm in sorted(flow_rows, key=lambda x: x[0]):
        if pns == "Purifying -> Positive selection":
            color = COLORS["pns_pos"]
        elif pns == "Positive -> Purifying selection":
            color = COLORS["pns_neg"]
        else:
            color = COLORS["single"]
        alpha = 0.86 if is_dual else 0.22
        draw_ribbon(
            ax,
            pns_arcs[(pns, rm)],
            rm_arcs[(rm, pns)],
            color=color,
            alpha=alpha,
            zorder=3 if is_dual else 1,
        )

    # Outer ring nodes.
    for node in PNS_ORDER + RM_ORDER:
        angle = ANGLES[node]
        x, y = polar_xy(angle, 1.0)
        if node.startswith("Purifying"):
            face = COLORS["pns_pos"]
        elif node.startswith("Positive"):
            face = COLORS["pns_neg"]
        elif node.startswith("r/m"):
            face = COLORS["rm"]
        else:
            face = "#EEEEEE"
        ax.add_patch(
            Wedge(
                (0, 0),
                1.07,
                *node_arc_bounds(node),
                width=0.12,
                facecolor=face,
                edgecolor="black",
                linewidth=0.8,
                alpha=0.96,
                zorder=4,
            )
        )
        lx, ly = polar_xy(angle, 1.30)
        ha = "left" if lx > 0.15 else "right" if lx < -0.15 else "center"
        ax.text(
            lx,
            ly,
            f"{clean_label(node)}\n(n={node_totals[node]})",
            ha=ha,
            va="center",
            fontsize=10.0,
            fontweight="bold" if "No " not in node else "normal",
            linespacing=1.05,
        )

    ax.text(-0.62, 1.22, "pN/pS switch", ha="center", va="center", fontsize=12.5, fontweight="bold")
    ax.text(0.68, 1.22, "r/m shift", ha="center", va="center", fontsize=12.5, fontweight="bold")
    ax.set_xlim(-1.58, 1.58)
    ax.set_ylim(-1.42, 1.46)

    ax2 = fig.add_subplot(gs[1])
    ax2.axis("off")
    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.text(
        0.5,
        0.96,
        "Genes with coupled pN/pS and r/m changes",
        ha="center",
        va="top",
        fontsize=12.0,
        fontweight="bold",
    )

    boxes = [
        ("Purifying -> Positive\n+ r/m decreased", ("r/m decreased", "Purifying -> Positive selection")),
        ("Purifying -> Positive\n+ r/m increased", ("r/m increased", "Purifying -> Positive selection")),
        ("Positive -> Purifying\n+ r/m decreased", ("r/m decreased", "Positive -> Purifying selection")),
        ("Positive -> Purifying\n+ r/m increased", ("r/m increased", "Positive -> Purifying selection")),
    ]
    x_positions = [0.02, 0.265, 0.51, 0.755]
    box_w = 0.225
    for x, (title, key) in zip(x_positions, boxes):
        rm, pns = key
        color = COLORS[key]
        sub = dual[(dual["rm_state"] == rm) & (dual["pns_state"] == pns)].copy()
        sub = sub.sort_values(["posterior_switch_probability", "gene"], ascending=[False, True])
        genes = [gene_label(r) for _, r in sub.iterrows()]
        gene_text = "\n".join(textwrap.shorten(g, width=28, placeholder="...") for g in genes) if genes else "None"
        ax2.add_patch(
            FancyBboxPatch(
                (x, 0.08),
                box_w,
                0.74,
                boxstyle="round,pad=0.008,rounding_size=0.008",
                facecolor=color,
                edgecolor="#BDBDBD",
                linewidth=0.8,
                alpha=0.24,
            )
        )
        ax2.add_patch(Rectangle := FancyBboxPatch(
            (x, 0.78),
            box_w,
            0.04,
            boxstyle="round,pad=0,rounding_size=0.004",
            facecolor=color,
            edgecolor=color,
            linewidth=0.0,
            alpha=0.95,
        ))
        ax2.text(x + 0.012, 0.72, f"{title} (n={len(sub)})", ha="left", va="top", fontsize=8.8, fontweight="bold")
        ax2.text(x + 0.012, 0.50, gene_text, ha="left", va="top", fontsize=8.4, linespacing=1.22)

    legend = [
        Patch(facecolor="#E59673", edgecolor="black", label="Purifying -> Positive selection"),
        Patch(facecolor="#8FA7DB", edgecolor="black", label="Positive -> Purifying selection"),
        Patch(facecolor="#CFCFCF", edgecolor="black", label="Single-axis shift"),
    ]
    fig.legend(
        handles=legend,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.006),
        ncol=3,
        frameon=False,
        fontsize=9.6,
        handlelength=1.0,
        columnspacing=1.3,
    )

    png = OUT_DIR / "coupled_rm_pns_circular.png"
    pdf = OUT_DIR / "coupled_rm_pns_circular.pdf"
    fig.savefig(png, dpi=350, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    print(png)
    print(pdf)
    print(OUT_DIR / "coupled_rm_pns_circular_double_shift_genes.csv")
    print(OUT_DIR / "coupled_rm_pns_circular_counts.csv")


if __name__ == "__main__":
    main()
