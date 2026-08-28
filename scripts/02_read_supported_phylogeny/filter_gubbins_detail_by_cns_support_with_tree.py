#!/usr/bin/env python3
"""
Filter Gubbins gene_recombination_detail.txt by sample-level CNS support.

Inputs:
  1. core_problematic_sites.site_calls.csv
     Expected columns, either with header:
       sample, position, ..., assembly_base, cns_base, ..., call
     or headerless format similar to:
       sample,position,locus_tag,gene,product,ref,assembly_base,cns_base,coverage,call

  2. gene_recombination_detail.txt
     Required columns:
       node, genome_pos, alt

  3. Gubbins node-labelled tree
     Newick tree with internal names such as Node_641 and terminal isolate names.

Filtering logic:
  - Terminal branch rows:
      keep if (sample, genome_pos) has call == supported.

  - Internal node rows:
      find descendant terminal isolates of that internal node.
      keep if enough descendants have a supported call at genome_pos.
      Default: position-supported descendants / all descendants >= 0.80
      You can change this with --internal-min-fraction.

Outputs:
  - filtered gene_recombination_detail.txt
  - summary CSV
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from io import StringIO
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from Bio import Phylo


def clean_sample(name: str) -> str:
    name = str(name).strip()
    for suffix in (".contigs", ".fasta", ".fa", ".fna", ".aln", ".cns"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
    return name


def clean_base(value: str) -> str:
    value = str(value or "").strip().upper()
    if "/" in value:
        value = value.split("/")[0]
    return value


def detect_site_call_columns(header: List[str]) -> Optional[Dict[str, int]]:
    lower = [x.strip().lower() for x in header]
    aliases = {
        "sample": ["sample", "isolate", "strain"],
        "position": ["position", "genome_pos", "pos"],
        "assembly_base": ["assembly_base", "assembly", "aln_base"],
        "cns_base": ["cns_base", "cns"],
        "call": ["call", "status"],
    }
    out = {}
    for key, names in aliases.items():
        for name in names:
            if name in lower:
                out[key] = lower.index(name)
                break
    if set(out) == set(aliases):
        return out
    return None


def read_site_calls(path: Path, supported_label: str) -> Tuple[Dict[Tuple[str, str], List[dict]], Counter]:
    """
    Return (sample, position) -> list of supported call records.
    """
    supported_label = supported_label.lower()
    supported: Dict[Tuple[str, str], List[dict]] = defaultdict(list)
    counts = Counter()

    with path.open(newline="", errors="ignore") as handle:
        reader = csv.reader(handle)
        first = next(reader)
        cols = detect_site_call_columns(first)

        if cols is None:
            # Headerless format from screenshot:
            # sample, position, locus_tag, gene, product, ref, assembly_base, cns_base, depth, call
            cols = {
                "sample": 0,
                "position": 1,
                "assembly_base": 6,
                "cns_base": 7,
                "call": 9,
            }
            rows: Iterable[List[str]] = [first]
        else:
            rows = []

        for row in list(rows) + list(reader):
            if len(row) <= max(cols.values()):
                continue
            sample = clean_sample(row[cols["sample"]])
            pos = str(row[cols["position"]]).strip()
            assembly_base = clean_base(row[cols["assembly_base"]])
            cns_base = clean_base(row[cols["cns_base"]])
            call = str(row[cols["call"]]).strip().lower()
            counts[call] += 1

            if call != supported_label:
                continue

            supported[(sample, pos)].append(
                {
                    "sample": sample,
                    "position": pos,
                    "assembly_base": assembly_base,
                    "cns_base": cns_base,
                    "call": call,
                }
            )

    return supported, counts


def read_tree(tree_path: Path):
    text = tree_path.read_text().strip()
    return Phylo.read(StringIO(text), "newick")


def build_tree_maps(tree) -> Tuple[dict, dict]:
    clades_by_name = {}
    descendant_samples = {}
    for clade in tree.find_clades(order="preorder"):
        if clade.name:
            name = str(clade.name).strip()
            clades_by_name[name] = clade
            descendant_samples[name] = [clean_sample(t.name) for t in clade.get_terminals() if t.name]
    return clades_by_name, descendant_samples


def supported_position(records: List[dict]) -> bool:
    return bool(records)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--site-calls", required=True, type=Path)
    parser.add_argument("--detail", required=True, type=Path)
    parser.add_argument("--tree", required=True, type=Path)
    parser.add_argument("--output-detail", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--supported-label", default="supported")
    parser.add_argument("--internal-min-fraction", default=0.80, type=float,
                        help="Minimum fraction of descendant isolates with a supported call at the position.")
    parser.add_argument("--internal-min-supported", default=1, type=int,
                        help="Minimum number of descendant isolates with a supported call at the position.")
    parser.add_argument("--keep-unknown-node", action="store_true",
                        help="Keep rows whose node cannot be found in the tree instead of removing them.")
    args = parser.parse_args()

    args.output_detail.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)

    supported, site_status_counts = read_site_calls(args.site_calls, args.supported_label)
    tree = read_tree(args.tree)
    clades_by_name, descendant_samples = build_tree_maps(tree)

    counters = Counter()
    removal_reasons = Counter()
    by_region = Counter()
    by_node_type = Counter()

    with args.detail.open(errors="ignore") as fin, args.output_detail.open("w", newline="") as fout:
        reader = csv.DictReader(fin, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit("detail file has no header")

        required = {"node", "genome_pos"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise SystemExit("detail file missing required columns: " + ", ".join(sorted(missing)))

        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames, delimiter="\t")
        writer.writeheader()

        for row in reader:
            counters["detail_total_rows"] += 1
            node = str(row["node"]).strip()
            pos = str(row["genome_pos"]).strip()
            region = str(row.get("region", "NA")).strip()
            by_region[f"input_{region}"] += 1

            keep = False
            reason = ""

            if node.startswith("Node_"):
                by_node_type["input_internal"] += 1
                desc = descendant_samples.get(node)
                if desc is None:
                    if args.keep_unknown_node:
                        keep = True
                        reason = "kept_unknown_internal_node"
                    else:
                        reason = "removed_unknown_internal_node"
                else:
                    supported_n = 0
                    for sample in desc:
                        records = supported.get((sample, pos), [])
                        if supported_position(records):
                            supported_n += 1
                    frac = supported_n / len(desc) if desc else 0.0
                    if supported_n >= args.internal_min_supported and frac >= args.internal_min_fraction:
                        keep = True
                        reason = "kept_internal_supported"
                    else:
                        reason = "removed_internal_insufficient_supported_descendants"
            else:
                by_node_type["input_terminal"] += 1
                sample = clean_sample(node)
                records = supported.get((sample, pos), [])
                if supported_position(records):
                    keep = True
                    reason = "kept_terminal_supported"
                else:
                    reason = "removed_terminal_not_supported"

            if keep:
                writer.writerow(row)
                counters["kept_rows"] += 1
                by_region[f"kept_{region}"] += 1
            else:
                counters["removed_rows"] += 1
                by_region[f"removed_{region}"] += 1
                removal_reasons[reason] += 1

            counters[reason] += 1

    with args.summary.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        writer.writerow(["site_calls_supported_pairs", len(supported)])
        for key, value in sorted(site_status_counts.items()):
            writer.writerow([f"site_calls_status_{key}", value])
        writer.writerow(["internal_min_fraction", args.internal_min_fraction])
        writer.writerow(["internal_min_supported", args.internal_min_supported])
        for key, value in sorted(counters.items()):
            writer.writerow([key, value])
        for key, value in sorted(removal_reasons.items()):
            writer.writerow([f"removal_reason_{key}", value])
        for key, value in sorted(by_node_type.items()):
            writer.writerow([key, value])
        for key, value in sorted(by_region.items()):
            writer.writerow([key, value])

    print("Finished")
    print(args.output_detail)
    print(args.summary)


if __name__ == "__main__":
    main()
