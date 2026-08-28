#!/usr/bin/env python3

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from statistics import median

from Bio import Phylo, SeqIO


DEFAULT_NODES = {
    "DCC1": "Node_641",
    "DCC2": "Node_147",
    "DCC3": "Node_845",
    "DCC4": "Node_271",
    "DCC5": "Node_362",
    "DCC6": "Node_1049",
    "DCC7": "Node_931",
}


def merge_intervals(intervals, max_gap=0):
    if not intervals:
        return []
    ordered = sorted(intervals)
    merged = [list(ordered[0])]
    for start, end in ordered[1:]:
        if start <= merged[-1][1] + max_gap + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [tuple(x) for x in merged]


def load_pass_intervals(path):
    by_chrom = defaultdict(list)
    raw_rows = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["flank_status"] != "PASS":
                continue
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])
            row["relative_depth"] = float(row["relative_depth"])
            by_chrom[row["chrom"]].append((row["start"], row["end"]))
            raw_rows.append(row)
    return {chrom: merge_intervals(v) for chrom, v in by_chrom.items()}, raw_rows


def consensus_segments(sample_intervals, threshold_count, max_gap=0):
    changes = defaultdict(int)
    for intervals in sample_intervals:
        for start, end in intervals:
            changes[start] += 1
            changes[end + 1] -= 1
    if not changes:
        return []

    segments = []
    depth = 0
    prev = None
    for pos in sorted(changes):
        if prev is not None and prev < pos and depth >= threshold_count:
            segments.append((prev, pos - 1))
        depth += changes[pos]
        prev = pos
    return merge_intervals(segments, max_gap=max_gap)


def feature_name(feature):
    q = feature.qualifiers
    locus = (q.get("locus_tag") or [""])[0]
    gene = (q.get("gene") or [""])[0]
    product = (q.get("product") or [""])[0]
    label = gene or locus or feature.type
    return locus, gene, product, label


def load_annotations(genbank):
    annotations = defaultdict(list)
    for record in SeqIO.parse(genbank, "genbank"):
        for feature in record.features:
            # GenBank commonly contains both a gene feature and its matching CDS.
            # Keeping both would duplicate every protein-coding locus in the output.
            if feature.type not in {"CDS", "rRNA", "tRNA", "ncRNA", "tmRNA"}:
                continue
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            locus, gene, product, label = feature_name(feature)
            annotations[record.id].append({
                "start": start,
                "end": end,
                "strand": "+" if feature.location.strand == 1 else "-" if feature.location.strand == -1 else "",
                "type": feature.type,
                "locus_tag": locus,
                "gene": gene,
                "product": product,
                "label": label,
            })
        annotations[record.id].sort(key=lambda x: (x["start"], x["end"]))
    return annotations


def annotate_region(chrom, start, end, annotations):
    features = annotations.get(chrom, [])
    overlaps = [f for f in features if f["start"] <= end and f["end"] >= start]
    if overlaps:
        full = [f for f in overlaps if start <= f["start"] and end >= f["end"]]
        partial = [f for f in overlaps if f not in full]
        region_type = "gene_overlap"
        if full and partial:
            region_type = "multi_gene_full_and_partial"
        elif len(full) > 1:
            region_type = "multi_gene_full"
        elif full:
            region_type = "full_gene"
        elif len(partial) > 1:
            region_type = "multi_gene_partial"
        else:
            region_type = "partial_gene"
        return {
            "annotation_type": region_type,
            "affected_locus_tags": "|".join(f["locus_tag"] for f in overlaps if f["locus_tag"]),
            "affected_genes": "|".join(f["gene"] for f in overlaps if f["gene"]),
            "affected_products": "|".join(f["product"] for f in overlaps if f["product"]),
            "fully_deleted_locus_tags": "|".join(f["locus_tag"] for f in full if f["locus_tag"]),
            "partially_affected_locus_tags": "|".join(f["locus_tag"] for f in partial if f["locus_tag"]),
            "left_flanking_gene": "",
            "right_flanking_gene": "",
        }

    left = max((f for f in features if f["end"] < start), key=lambda f: f["end"], default=None)
    right = min((f for f in features if f["start"] > end), key=lambda f: f["start"], default=None)
    return {
        "annotation_type": "intergenic",
        "affected_locus_tags": "",
        "affected_genes": "",
        "affected_products": "",
        "fully_deleted_locus_tags": "",
        "partially_affected_locus_tags": "",
        "left_flanking_gene": left["gene"] or left["locus_tag"] if left else "",
        "right_flanking_gene": right["gene"] or right["locus_tag"] if right else "",
    }


def matching_intervals(intervals, core_start, core_end):
    return [(s, e) for s, e in intervals if s <= core_end and e >= core_start]


def main():
    parser = argparse.ArgumentParser(description="Summarize DCC deletions carried by at least a specified fraction of tree-defined members.")
    parser.add_argument("--result-dir", required=True, type=Path)
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--genbank", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--carriage", type=float, default=0.95)
    parser.add_argument("--min-length", type=int, default=20)
    parser.add_argument("--consensus-gap", type=int, default=0)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    tree = Phylo.read(args.tree, "newick")
    named = {c.name: c for c in tree.find_clades() if c.name}
    result_files = {p.stem: p for p in args.result_dir.glob("*.tsv")}
    annotations = load_annotations(args.genbank)

    all_rows = []
    qc_rows = []
    for dcc, node in DEFAULT_NODES.items():
        tips = [x.name for x in named[node].get_terminals()]
        available = [x for x in tips if x in result_files]
        missing = sorted(set(tips) - set(available))
        threshold_count = math.ceil(args.carriage * len(available))
        sample_data = {}
        raw_data = {}
        for sample in available:
            sample_data[sample], raw_data[sample] = load_pass_intervals(result_files[sample])

        chroms = sorted({chrom for data in sample_data.values() for chrom in data})
        dcc_rows = []
        event_number = 0
        for chrom in chroms:
            interval_sets = [sample_data[s].get(chrom, []) for s in available]
            cores = consensus_segments(interval_sets, threshold_count, args.consensus_gap)
            for core_start, core_end in cores:
                if core_end - core_start + 1 < args.min_length:
                    continue
                matched = {}
                for sample in available:
                    hits = matching_intervals(sample_data[sample].get(chrom, []), core_start, core_end)
                    if hits:
                        matched[sample] = max(hits, key=lambda x: min(x[1], core_end) - max(x[0], core_start) + 1)
                if len(matched) < threshold_count:
                    continue

                event_number += 1
                starts = [x[0] for x in matched.values()]
                ends = [x[1] for x in matched.values()]
                representative_start = int(median(starts))
                representative_end = int(median(ends))
                rel_depths = []
                for sample, (s, e) in matched.items():
                    for row in raw_data[sample]:
                        if row["chrom"] == chrom and row["start"] <= core_end and row["end"] >= core_start:
                            rel_depths.append(row["relative_depth"])

                row = {
                    "DCC": dcc,
                    "tree_node": node,
                    "event_id": f"{dcc}_DEL_{event_number:03d}",
                    "chrom": chrom,
                    "consensus_start": core_start,
                    "consensus_end": core_end,
                    "consensus_length": core_end - core_start + 1,
                    "representative_start": representative_start,
                    "representative_end": representative_end,
                    "representative_length": representative_end - representative_start + 1,
                    "observed_start_min": min(starts),
                    "observed_start_max": max(starts),
                    "observed_end_min": min(ends),
                    "observed_end_max": max(ends),
                    "start_breakpoint_span": max(starts) - min(starts),
                    "end_breakpoint_span": max(ends) - min(ends),
                    "n_tree_tips": len(tips),
                    "n_result_files": len(available),
                    "n_carriers": len(matched),
                    "carriage_pct": round(100 * len(matched) / len(available), 4),
                    "median_relative_depth": round(median(rel_depths), 6) if rel_depths else "",
                    "carrier_samples": "|".join(sorted(matched)),
                }
                tolerance = max(50, round(0.10 * max(1, representative_end - representative_start + 1)))
                row["breakpoint_consistency"] = (
                    "stable" if max(starts) - min(starts) <= tolerance and max(ends) - min(ends) <= tolerance
                    else "variable"
                )
                row.update(annotate_region(chrom, representative_start, representative_end, annotations))
                dcc_rows.append(row)
                all_rows.append(row)

        qc_rows.append({
            "DCC": dcc,
            "tree_node": node,
            "n_tree_tips": len(tips),
            "n_result_files": len(available),
            "n_missing_result_files": len(missing),
            "missing_samples": "|".join(missing),
            "minimum_carriers_for_95pct": threshold_count,
            "n_consensus_deletions": len(dcc_rows),
            "n_stable_breakpoint_deletions": sum(r["breakpoint_consistency"] == "stable" for r in dcc_rows),
        })
        write_csv(args.output_dir / f"{dcc}_deletions_carried_by_95pct.csv", dcc_rows)

    write_csv(args.output_dir / "all_DCC_deletions_carried_by_95pct.csv", all_rows)
    write_csv(args.output_dir / "DCC_deletion_analysis_qc_summary.csv", qc_rows)
    if all_rows:
        event_fields = list(all_rows[0])
        for dcc in DEFAULT_NODES:
            path = args.output_dir / f"{dcc}_deletions_carried_by_95pct.csv"
            if path.stat().st_size == 0:
                write_csv(path, [], fieldnames=event_fields)
    with (args.output_dir / "README.txt").open("w") as out:
        out.write(
            "Only deletion candidates with flank_status=PASS were used.\n"
            "DCC membership was defined by descendants of the target nodes in the labelled tree.\n"
            "A consensus segment is retained when every base is covered by a PASS deletion in at least 95% of DCC result files.\n"
            "Coordinates are 1-based inclusive. Representative boundaries are medians of the matching strain-level intervals.\n"
            "The consensus interval is the conservative shared deletion core; representative boundaries describe the typical full event.\n"
        )


def write_csv(path, rows, fieldnames=None):
    if not rows:
        if fieldnames:
            with path.open("w", newline="") as handle:
                csv.DictWriter(handle, fieldnames=fieldnames).writeheader()
        else:
            path.write_text("")
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
