#!/usr/bin/env python3
"""Plot repeat-cleaned tree-based indel events along the MAB genome."""

from pathlib import Path
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from Bio import SeqIO


DCC_COLUMNS = [f"DCC{i}_event_count" for i in range(1, 8)]


def gene_positions(genbank: Path) -> pd.DataFrame:
    rows = []
    for record in SeqIO.parse(genbank, "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            tags = feature.qualifiers.get("locus_tag", [])
            if not tags:
                continue
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            rows.append({
                "locus_tag": tags[0],
                "gene_position": (start + end) / 2,
                "gene_start": start,
                "gene_end": end,
            })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--genbank", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--min-dccs", type=int, default=6)
    parser.add_argument("--label-min-events", type=int, default=70)
    args = parser.parse_args()

    counts = pd.read_csv(args.counts)
    counts["n_DCCs_with_event"] = (counts[DCC_COLUMNS] > 0).sum(axis=1)
    counts["highlight_ge6_DCCs"] = (
        (counts["n_DCCs_with_event"] >= args.min_dccs)
        & (counts["total_event_count"] >= args.label_min_events)
    )
    positions = gene_positions(args.genbank)
    df = counts.merge(positions, on="locus_tag", how="inner")
    highlighted = df[df.highlight_ge6_DCCs].copy()

    plt.rcParams.update({"font.family": "Arial", "font.size": 14})
    fig, ax = plt.subplots(figsize=(5.4, 5.85))

    background = df[~df.highlight_ge6_DCCs]
    background_sizes = 22 + background.n_DCCs_with_event * 9
    highlighted_sizes = 22 + highlighted.n_DCCs_with_event * 9
    ax.scatter(
        background.gene_position, background.total_event_count,
        s=background_sizes, facecolors="white", edgecolors="#79BFE8",
        linewidths=1.35, alpha=0.82, zorder=1,
        label="Other genes",
    )
    ax.scatter(
        highlighted.gene_position, highlighted.total_event_count,
        s=highlighted_sizes, facecolors="white", edgecolors="#E56F68",
        linewidths=1.65, alpha=0.95, zorder=3,
        label=f"Present in >={args.min_dccs} DCCs",
    )

    ax.set_xlabel("Gene position")
    ax.set_ylabel("Independent indel events")
    x_ticks = [0, 2_000_000, 4_000_000]
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([
        "0" if value == 0 else rf"${int(value / 1_000_000)}\times10^{{6}}$"
        for value in x_ticks
    ])
    ax.set_xlim(0, positions.gene_end.max() * 1.015)
    ax.set_ylim(-2, df.total_event_count.max() * 1.28)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.25)
    ax.spines["bottom"].set_linewidth(1.25)
    ax.tick_params(width=1.1, length=5)
    color_legend = ax.legend(
        frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0, fontsize=10.5, handletextpad=0.4
    )
    ax.add_artist(color_legend)
    size_handles = [
        ax.scatter([], [], s=22 + n * 9, facecolors="white", edgecolors="#777777",
                   linewidths=1.1, label=str(n))
        for n in (1, 4, 7)
    ]
    ax.legend(
        handles=size_handles, title="Number of DCCs", frameon=False,
        loc="upper left", bbox_to_anchor=(1.02, 0.68), borderaxespad=0,
        fontsize=9.5, title_fontsize=10.5,
        labelspacing=0.8, handletextpad=0.7,
    )
    fig.tight_layout()

    args.output_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(args.output_prefix) + ".png", dpi=600, bbox_inches="tight")
    fig.savefig(str(args.output_prefix) + ".pdf", bbox_inches="tight")
    df.to_csv(str(args.output_prefix) + ".plot_data.csv", index=False)
    highlighted.sort_values("total_event_count", ascending=False).to_csv(
        str(args.output_prefix) + ".highlighted_genes.csv", index=False
    )
    print(f"Plotted {len(df)} genes; highlighted {len(highlighted)} genes.")


if __name__ == "__main__":
    main()
