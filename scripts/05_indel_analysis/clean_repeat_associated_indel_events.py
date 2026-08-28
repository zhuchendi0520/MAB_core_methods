#!/usr/bin/env python3
"""Audit high-event genes and remove repeat-associated indel artifacts."""

from __future__ import annotations

import argparse
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from Bio import Phylo, SeqIO


DCCS = [f"DCC{i}" for i in range(1, 8)]


def primitive_motif(seq: str, max_period: int = 6) -> tuple[str, int]:
    seq = seq.upper()
    for period in range(1, min(max_period, len(seq)) + 1):
        motif = seq[:period]
        if (motif * math.ceil(len(seq) / period))[: len(seq)] == seq:
            return motif, period
    return seq, len(seq)


def max_homopolymer(seq: str) -> int:
    best = run = 0
    previous = ""
    for base in seq.upper():
        if base == previous:
            run += 1
        else:
            previous, run = base, 1
        best = max(best, run)
    return best


def max_tandem_copies(seq: str, max_period: int = 6) -> tuple[int, str]:
    seq = seq.upper()
    best_copies, best_motif = 1, ""
    for period in range(1, max_period + 1):
        for offset in range(period):
            i = offset
            while i + period <= len(seq):
                motif = seq[i : i + period]
                copies = 1
                j = i + period
                while j + period <= len(seq) and seq[j : j + period] == motif:
                    copies += 1
                    j += period
                if copies > best_copies:
                    best_copies, best_motif = copies, motif
                i = max(i + period, j)
    return best_copies, best_motif


def repeat_context(record_seq: str, pos: int, allele: str) -> dict:
    # Input coordinates are 1-based; indel representation is anchored at pos.
    left = max(0, pos - 21)
    right = min(len(record_seq), pos + 20)
    context = record_seq[left:right].upper()
    inserted_or_deleted = allele[1:].upper()
    motif, motif_period = primitive_motif(inserted_or_deleted)
    hp = max_homopolymer(context)
    tandem_copies, tandem_motif = max_tandem_copies(context)

    # Test whether the indel unit is already repeated immediately around its anchor.
    # VarScan-style indels are anchored on the reference base at `pos`;
    # inserted/deleted sequence begins immediately after that base.
    anchor0 = pos
    flank_left = record_seq[max(0, anchor0 - 30) : anchor0].upper()
    flank_right = record_seq[anchor0 : min(len(record_seq), anchor0 + 30)].upper()
    adjacent_left = flank_left.endswith(motif)
    adjacent_right = flank_right.startswith(motif)
    left_copies = 0
    cursor = len(flank_left)
    while cursor >= len(motif) and flank_left[cursor - len(motif) : cursor] == motif:
        left_copies += 1
        cursor -= len(motif)
    right_copies = 0
    cursor = 0
    while cursor + len(motif) <= len(flank_right) and flank_right[cursor : cursor + len(motif)] == motif:
        right_copies += 1
        cursor += len(motif)
    adjacent_repeat_copies = left_copies + right_copies

    homopolymer_artifact = motif_period == 1 and adjacent_repeat_copies >= 4
    short_tandem_artifact = (
        motif_period <= 6
        and adjacent_repeat_copies >= 3
    )
    reasons = []
    if homopolymer_artifact:
        reasons.append("homopolymer_run_ge5")
    if short_tandem_artifact:
        reasons.append("short_tandem_repeat_ge3")
    return {
        "sequence_context_41bp": context,
        "indel_unit": inserted_or_deleted,
        "primitive_motif": motif,
        "motif_period": motif_period,
        "max_homopolymer_run": hp,
        "max_tandem_copies": tandem_copies,
        "local_tandem_motif": tandem_motif,
        "motif_adjacent_left": adjacent_left,
        "motif_adjacent_right": adjacent_right,
        "adjacent_repeat_copies": adjacent_repeat_copies,
        "repeat_artifact": bool(reasons),
        "removal_reason": "|".join(reasons),
    }


def read_targets(path: Path) -> dict[str, str]:
    targets = {}
    for line in path.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 2 and fields[0].startswith("DCC"):
            targets[fields[0]] = fields[1]
    return {dcc: targets[dcc] for dcc in DCCS}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--indel-dir", type=Path, required=True)
    parser.add_argument("--tree", type=Path, required=True)
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--genbank", type=Path, required=True)
    parser.add_argument("--original-counts", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--audit-threshold", type=int, default=30)
    args = parser.parse_args()

    sys.path.insert(0, str(args.indel_dir))
    from analyze_DCC_indels import annotate_event, load_features, read_events
    from count_tree_indel_events_by_gene import sankoff_gain_count

    args.outdir.mkdir(parents=True, exist_ok=True)
    original = pd.read_csv(args.original_counts)
    high_genes = set(original.loc[original.total_event_count > args.audit_threshold, "locus_tag"])

    record = SeqIO.read(args.genbank, "genbank")
    refseq = str(record.seq).upper()
    _, _, features, starts = load_features(args.genbank)
    event_samples, _, assignments = read_events(args.indel_dir, {})
    available = set(assignments["sample"])

    targets = read_targets(args.targets)
    tree = Phylo.read(args.tree, "newick")
    named = {node.name: node for node in tree.find_clades() if node.name}
    roots = {dcc: named[node] for dcc, node in targets.items()}
    dcc_tips = {
        dcc: {tip.name for tip in root.get_terminals() if tip.name in available}
        for dcc, root in roots.items()
    }

    event_hits = {}
    gene_meta = {}
    strain_gene_events = defaultdict(set)
    for event, samples in event_samples.items():
        start, end, typ, length, hits = annotate_event(event, features, starts)
        tags = []
        for hit in hits:
            tag = hit["locus_tag"] or f"intergenic:{start}-{end}"
            tags.append(tag)
            gene_meta[tag] = hit
            for sample in samples:
                strain_gene_events[(sample, tag)].add(event)
        event_hits[event] = (tags, start, end, typ, length)

    # Preserve the prior QC rule: remove every sample-gene pair with >=2 distinct indels.
    excluded_samples_by_gene = defaultdict(set)
    for (sample, tag), events in strain_gene_events.items():
        if len(events) >= 2:
            excluded_samples_by_gene[tag].add(sample)

    audit_rows = []
    clean_counts = defaultdict(lambda: defaultdict(int))
    clean_meta = {}
    for event, samples in event_samples.items():
        pos, ref, allele = event
        tags, start, end, typ, length = event_hits[event]
        context = repeat_context(refseq, pos, allele)
        for tag in tags:
            clean_samples = samples - excluded_samples_by_gene[tag]
            event_counts = {}
            carrier_counts = {}
            for dcc in DCCS:
                carriers = clean_samples & dcc_tips[dcc]
                carrier_counts[dcc] = len(carriers)
                event_counts[dcc] = sankoff_gain_count(roots[dcc], carriers) if carriers else 0

            is_audited = tag in high_genes
            is_hotspot = is_audited and sum(event_counts.values()) >= 2
            remove = is_hotspot and context["repeat_artifact"]
            if not remove:
                clean_meta[tag] = gene_meta[tag]
                for dcc, count in event_counts.items():
                    clean_counts[tag][dcc] += count

            if is_audited:
                present_dccs = [dcc for dcc in DCCS if carrier_counts[dcc] > 0]
                row = {
                    "locus_tag": tag,
                    "gene": gene_meta[tag].get("gene", ""),
                    "position": pos,
                    "ref": ref,
                    "allele": allele,
                    "type": typ,
                    "event_length": length,
                    "affected_start": start,
                    "affected_end": end,
                    "n_clean_carriers": len(clean_samples),
                    "n_DCCs_with_carriers": len(present_dccs),
                    "DCCs_with_carriers": "|".join(present_dccs),
                    "cross_DCC_exact_overlap": len(present_dccs) >= 2,
                    "max_within_DCC_carriers": max(carrier_counts.values()),
                    "total_tree_event_count": sum(event_counts.values()),
                    "is_indel_hotspot_ge2_tree_events": is_hotspot,
                    "removed_from_clean_result": remove,
                    **context,
                }
                for dcc in DCCS:
                    row[f"{dcc}_carriers"] = carrier_counts[dcc]
                    row[f"{dcc}_tree_events"] = event_counts[dcc]
                audit_rows.append(row)

    rows = []
    for tag, counts in clean_counts.items():
        meta = clean_meta[tag]
        row = {"locus_tag": tag, "gene": meta.get("gene", "")}
        for dcc in DCCS:
            row[f"{dcc}_event_count"] = counts[dcc]
        row["total_event_count"] = sum(counts.values())
        if row["total_event_count"]:
            rows.append(row)
    clean = pd.DataFrame(rows).sort_values("total_event_count", ascending=False)
    clean = clean.merge(
        original[["locus_tag", "ref1_product", "ref2_product", "Type"]],
        on="locus_tag", how="left",
    )

    audit = pd.DataFrame(audit_rows).sort_values(
        ["removed_from_clean_result", "total_tree_event_count", "locus_tag", "position"],
        ascending=[False, False, True, True],
    )
    removed = audit[audit.removed_from_clean_result].copy()
    hotspots = audit[audit.is_indel_hotspot_ge2_tree_events].copy()
    clean.to_csv(args.outdir / "gene_indel_tree_event_counts_repeat_cleaned.csv", index=False)
    audit.to_csv(args.outdir / "high_event_genes_indel_repeat_audit.csv", index=False)
    hotspots.to_csv(args.outdir / "indel_hotspots_high_event_genes.csv", index=False)
    removed.to_csv(args.outdir / "removed_repeat_associated_indels.csv", index=False)

    old = original.set_index("locus_tag")["total_event_count"]
    new = clean.set_index("locus_tag")["total_event_count"]
    comparison = pd.DataFrame({"original_event_count": old, "clean_event_count": new}).fillna(0)
    comparison["events_removed"] = comparison.original_event_count - comparison.clean_event_count
    comparison = comparison.reset_index().sort_values("events_removed", ascending=False)
    comparison.to_csv(args.outdir / "gene_event_counts_before_after_repeat_filter.csv", index=False)

    summary = pd.DataFrame([
        {"metric": "genes_audited_original_count_gt_threshold", "value": len(high_genes)},
        {"metric": "audited_exact_indel_gene_records", "value": len(audit)},
        {"metric": "hotspot_records_ge2_tree_events", "value": len(hotspots)},
        {"metric": "repeat_associated_records_removed", "value": len(removed)},
        {"metric": "unique_repeat_associated_indels_removed", "value": removed[["position", "ref", "allele"]].drop_duplicates().shape[0]},
        {"metric": "original_total_tree_events", "value": int(original.total_event_count.sum())},
        {"metric": "clean_total_tree_events", "value": int(clean.total_event_count.sum())},
    ])
    summary.to_csv(args.outdir / "repeat_filter_summary.csv", index=False)
    print(summary.to_string(index=False))
    print("\nLargest reductions:")
    print(comparison.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
