#!/usr/bin/env python3
"""Plot absolute DCC pangenome composition with a compressed angular scale."""

from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cache")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


if len(sys.argv) != 3:
    raise SystemExit(
        "Usage: python 06_plot_DCC_pangenome_composition_absolute.py "
        "INPUT_CSV OUTPUT_DIR"
    )
INPUT = Path(sys.argv[1])
OUT_DIR = Path(sys.argv[2])
OUT_DIR.mkdir(parents=True, exist_ok=True)

DCC_ORDER = [f"DCC{i}" for i in range(1, 8)]
COMPONENTS = ["Core", "Soft Core", "Shell", "Cloud"]
COLORS = {
    "Core": "#717CB5",
    "Soft Core": "#8CBCD3",
    "Shell": "#D5E2B5",
    "Cloud": "#F1DDA2",
}

BREAK_VALUE = 6000.0
AXIS_MAX = 20000.0
ARC_RADIANS = np.deg2rad(300.0)
BREAK_ANGLE = np.pi
BREAK_GAP_OUTER_DEGREES = 4.5
BREAK_GAP_INNER_DEGREES = 8.0


def value_to_angle(value: float | np.ndarray) -> float | np.ndarray:
    """Piecewise angular transform that compresses values above 6,000."""
    values = np.asarray(value, dtype=float)
    angles = np.where(
        values <= BREAK_VALUE,
        values / BREAK_VALUE * BREAK_ANGLE,
        BREAK_ANGLE
        + (values - BREAK_VALUE)
        / (AXIS_MAX - BREAK_VALUE)
        * (ARC_RADIANS - BREAK_ANGLE),
    )
    return float(angles) if angles.ndim == 0 else angles


wide = pd.read_csv(INPUT)
required = ["DCC", *COMPONENTS]
missing = [column for column in required if column not in wide.columns]
if missing:
    raise SystemExit(f"Missing columns: {missing}")

wide = wide.set_index("DCC").reindex(DCC_ORDER)
counts = wide[COMPONENTS].apply(pd.to_numeric, errors="coerce")
if counts.isna().any().any():
    raise SystemExit("Composition columns contain missing or non-numeric values")
if (counts < 0).any().any():
    raise SystemExit("Composition columns contain negative values")
if counts.sum(axis=1).max() > AXIS_MAX:
    raise SystemExit("AXIS_MAX is smaller than the largest DCC total")

records: list[dict[str, float | str]] = []
for dcc in DCC_ORDER:
    cumulative = 0.0
    total = float(counts.loc[dcc].sum())
    for component in COMPONENTS:
        count = float(counts.loc[dcc, component])
        start = cumulative
        end = cumulative + count
        records.append(
            {
                "DCC": dcc,
                "Component": component,
                "Count": int(count),
                "Total": int(total),
                "Cumulative_start": int(start),
                "Cumulative_end": int(end),
                "Angle_start_degrees": np.degrees(value_to_angle(start)),
                "Angle_end_degrees": np.degrees(value_to_angle(end)),
            }
        )
        cumulative = end

plot_table = pd.DataFrame(records)
plot_table.to_csv(
    OUT_DIR / "DCC_pangenome_composition_absolute_broken_scale.csv",
    index=False,
)

fig = plt.figure(figsize=(8.5, 8.5), facecolor="white")
ax = fig.add_axes([0.08, 0.07, 0.86, 0.84], projection="polar")
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.set_axis_off()

inner_radius = 3.15
ring_height = 0.62
ring_gap = 0.21
ring_centers: dict[str, float] = {}
ring_bottoms: dict[str, float] = {}

for idx, dcc in enumerate(reversed(DCC_ORDER)):
    bottom = inner_radius + idx * (ring_height + ring_gap)
    ring_centers[dcc] = bottom + ring_height / 2
    ring_bottoms[dcc] = bottom
    dcc_rows = plot_table[plot_table["DCC"] == dcc]

    for row in dcc_rows.itertuples(index=False):
        theta_start = value_to_angle(row.Cumulative_start)
        theta_end = value_to_angle(row.Cumulative_end)
        ax.bar(
            theta_start,
            ring_height,
            width=theta_end - theta_start,
            bottom=bottom,
            align="edge",
            color=COLORS[row.Component],
            edgecolor="#262626",
            linewidth=1.05,
        )

outer_radius = inner_radius + 6 * (ring_height + ring_gap) + ring_height
guide_radius = outer_radius + 0.54
theta_guide = np.linspace(0, ARC_RADIANS, 500)
ax.plot(theta_guide, np.full_like(theta_guide, guide_radius), color="#262626", lw=0.8)

tick_values = [0, 2000, 4000, 6000, 10000, 14000, 18000, 20000]
for value in tick_values:
    theta = value_to_angle(value)
    ax.plot(
        [theta, theta],
        [guide_radius - 0.08, guide_radius + 0.08],
        color="#262626",
        lw=0.8,
    )
    ax.text(
        theta,
        guide_radius + 0.36,
        f"{value:,}",
        fontsize=10,
        ha="center",
        va="center",
        family="Helvetica",
        rotation=(((-np.degrees(theta)) + 90) % 180) - 90,
        rotation_mode="anchor",
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.08},
    )

ax.set_ylim(0, guide_radius + 0.75)

# Cut a true white break into each ring. The inner/outer edges are white,
# while the two angular boundaries remain black.
for dcc in DCC_ORDER:
    bottom = ring_bottoms[dcc]
    outer_to_inner_index = DCC_ORDER.index(dcc)
    gap_degrees = np.interp(
        outer_to_inner_index,
        [0, len(DCC_ORDER) - 1],
        [BREAK_GAP_OUTER_DEGREES, BREAK_GAP_INNER_DEGREES],
    )
    gap_angle = np.deg2rad(gap_degrees)
    left_edge = BREAK_ANGLE - gap_angle / 2
    right_edge = BREAK_ANGLE + gap_angle / 2
    ax.bar(
        left_edge,
        ring_height,
        width=gap_angle,
        bottom=bottom,
        align="edge",
        color="white",
        edgecolor="white",
        linewidth=1.8,
        zorder=8,
    )
    for theta_edge in [left_edge, right_edge]:
        ax.plot(
            [theta_edge, theta_edge],
            [bottom, bottom + ring_height],
            color="#262626",
            linewidth=1.05,
            zorder=9,
        )
    ax.text(
        BREAK_ANGLE,
        ring_centers[dcc],
        "//",
        fontsize=9.5,
        fontweight="bold",
        family="Helvetica",
        ha="center",
        va="center",
        rotation=-18,
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.04},
        zorder=10,
    )

# Also mark the corresponding position on the outer guide without changing it.
ax.text(
    BREAK_ANGLE,
    guide_radius,
    "//",
    fontsize=10.5,
    fontweight="bold",
    family="Helvetica",
    ha="center",
    va="center",
    rotation=-18,
    bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.10},
    zorder=10,
)

fig.canvas.draw()
for dcc in DCC_ORDER:
    display_xy = ax.transData.transform((0, ring_centers[dcc]))
    figure_xy = fig.transFigure.inverted().transform(display_xy)
    fig.text(
        figure_xy[0] - 0.018,
        figure_xy[1],
        dcc,
        ha="right",
        va="center",
        fontsize=13,
        fontweight="bold",
        family="Helvetica",
        color="#202020",
    )

legend_handles = [
    Patch(facecolor=COLORS[component], edgecolor="#262626", label=component)
    for component in COMPONENTS
]
fig.legend(
    handles=legend_handles,
    loc="upper left",
    bbox_to_anchor=(0.225, 0.805),
    frameon=False,
    fontsize=12,
    handlelength=1.25,
    handleheight=1.25,
    borderaxespad=0,
    labelspacing=0.38,
)

fig.text(
    0.185,
    0.555,
    "//  Scale compressed\n     above 6,000 genes",
    fontsize=9.5,
    family="Helvetica",
    color="#303030",
)

fig.suptitle(
    "Absolute pangenome composition across DCCs",
    x=0.52,
    y=0.965,
    fontsize=18,
    fontweight="bold",
    family="Helvetica",
)

png_path = OUT_DIR / "DCC_pangenome_composition_absolute_broken_polar.png"
pdf_path = OUT_DIR / "DCC_pangenome_composition_absolute_broken_polar.pdf"
fig.savefig(png_path, dpi=400, facecolor="white")
fig.savefig(pdf_path, facecolor="white")
plt.close(fig)

print(png_path)
print(pdf_path)
