#!/usr/bin/env python3
"""
Extract assembly-vs-CNS site-level calls inside selected core genes.

This script uses the same inputs as core_gene_assembly_vs_mapping_simple_summary.py:
  --reference-gb
  --alignment-fasta
  --cns-dir
  --gene-list

For each sample in all.aln, it finds A/C/G/T positions where the assembly
alignment differs from ATCC19977 inside selected CDS features. It then checks
the sample CNS file and classifies each site as:
  supported     CNS consensus equals the all.aln assembly base
  inconsistent CNS exists with enough coverage, but consensus differs
  missing       no CNS record or coverage < --min-cov

Default output keeps only inconsistent and missing sites for enrichment/QC.
Use --include-supported to also keep supported sites.

Outputs:
  <output-prefix>.site_calls.csv
  <output-prefix>.gene_summary.csv

In gene_summary.csv, the count columns are pooled counts across isolates, while
the rate columns are calculated per isolate first and then averaged across
isolates with at least one assembly-derived SNP in that gene. Pooled rates are
also reported separately for comparison.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Set, Tuple


VALID_BASES = set("ACGT")


def norm_base(value: str) -> str:
    value = (value or "").strip().upper()
    if not value or value == ".":
        return "."
    if "/" in value:
        value = value.split("/")[0]
    return value


def read_gene_list(path: Path) -> Set[str]:
    genes: Set[str] = set()
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            genes.add(line.split()[0])
    return genes


def read_genbank_origin(path: Path) -> str:
    seq: List[str] = []
    in_origin = False
    with path.open(errors="ignore") as handle:
        for line in handle:
            if line.startswith("ORIGIN"):
                in_origin = True
                continue
            if in_origin:
                if line.startswith("//"):
                    break
                parts = line.strip().split()
                if parts:
                    seq.extend(parts[1:])
    if not seq:
        raise ValueError(f"No ORIGIN sequence found in {path}")
    return "".join(seq).upper()


def iter_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    name: Optional[str] = None
    chunks: List[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def sample_stem(name: str) -> str:
    for suffix in (".contigs", ".fasta", ".fa", ".fna", ".aln"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name


def find_cns(sample_name: str, cns_dir: Path) -> Optional[Path]:
    candidates = [
        cns_dir / f"{sample_name}.cns",
        cns_dir / f"{sample_name}.contigs.cns",
        cns_dir / f"{sample_stem(sample_name)}.cns",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    matches = sorted(cns_dir.glob(f"{sample_name}*.cns"))
    if not matches:
        matches = sorted(cns_dir.glob(f"{sample_stem(sample_name)}*.cns"))
    return matches[0] if matches else None


def parse_location_ranges(location: str) -> List[Tuple[int, int]]:
    """Parse simple GenBank CDS locations into inclusive 1-based ranges."""
    nums = re.findall(r"(\d+)\.\.(\d+)", location)
    ranges: List[Tuple[int, int]] = []
    for start, end in nums:
        s, e = int(start), int(end)
        if s > e:
            s, e = e, s
        ranges.append((s, e))
    if not ranges:
        single = re.findall(r"(?<!\.)\b(\d+)\b(?!\.)", location)
        for value in single:
            pos = int(value)
            ranges.append((pos, pos))
    return ranges


def extract_qualifier(block: List[str], key: str) -> Optional[str]:
    pattern = f"/{key}="
    for line in block:
        text = line.strip()
        if text.startswith(pattern):
            value = text[len(pattern):].strip()
            return value.strip('"')
    return None


def selected_cds_position_map(
    gb_path: Path,
    selected_genes: Set[str],
) -> Tuple[Dict[int, List[dict]], Dict[str, dict]]:
    """
    Return:
      position_map: genomic position -> list of selected CDS annotations
      gene_info: locus_tag -> gene-level annotation and length
    """
    position_map: Dict[int, List[dict]] = defaultdict(list)
    gene_info: Dict[str, dict] = {}
    current: List[str] = []

    def flush_feature(feature_lines: List[str]) -> None:
        if not feature_lines:
            return
        first = feature_lines[0]
        if not first.startswith("     CDS"):
            return
        location = first[21:].strip()
        locus_tag = extract_qualifier(feature_lines, "locus_tag") or ""
        gene = extract_qualifier(feature_lines, "gene") or ""
        old_locus_tag = extract_qualifier(feature_lines, "old_locus_tag") or ""
        product = extract_qualifier(feature_lines, "product") or ""
        names = {x for x in (locus_tag, gene, old_locus_tag) if x}
        if not names.intersection(selected_genes):
            return

        ranges = parse_location_ranges(location)
        cds_len = sum(end - start + 1 for start, end in ranges)
        key = locus_tag or gene or old_locus_tag
        annot = {
            "locus_tag": locus_tag,
            "gene": gene,
            "old_locus_tag": old_locus_tag,
            "product": product,
            "location": location,
            "cds_length": cds_len,
        }
        gene_info[key] = annot
        for start, end in ranges:
            for pos in range(start, end + 1):
                position_map[pos].append(annot)

    with gb_path.open(errors="ignore") as handle:
        in_features = False
        for line in handle:
            if line.startswith("FEATURES"):
                in_features = True
                continue
            if in_features and line.startswith("ORIGIN"):
                flush_feature(current)
                break
            if not in_features:
                continue
            if line.startswith("     ") and not line.startswith("                     "):
                flush_feature(current)
                current = [line.rstrip("\n")]
            elif current:
                current.append(line.rstrip("\n"))

    return dict(position_map), gene_info


def read_cns_for_positions(
    cns_path: Path,
    positions: Set[int],
    min_cov: int,
) -> Dict[int, Tuple[Optional[str], int]]:
    """Return position -> (consensus base, coverage); consensus None means low coverage."""
    out: Dict[int, Tuple[Optional[str], int]] = {}
    with cns_path.open() as handle:
        next(handle, None)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            try:
                pos = int(parts[1])
            except ValueError:
                continue
            if pos not in positions:
                continue
            fields = parts[4].split(":")
            cons = norm_base(fields[0])
            cov = 0
            if len(fields) >= 2:
                try:
                    cov = int(float(fields[1]))
                except ValueError:
                    cov = 0
            out[pos] = (cons if cov >= min_cov else None, cov)
    return out


def classify_site(
    pos: int,
    aln_base: str,
    cns: Dict[int, Tuple[Optional[str], int]],
) -> Tuple[str, str, int]:
    if pos not in cns:
        return "missing", "", 0
    cns_base, cov = cns[pos]
    if cns_base is None:
        return "missing", "", cov
    if cns_base == aln_base:
        return "supported", cns_base, cov
    return "inconsistent", cns_base, cov


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-gb", required=True, type=Path)
    parser.add_argument("--alignment-fasta", required=True, type=Path)
    parser.add_argument("--cns-dir", required=True, type=Path)
    parser.add_argument("--gene-list", required=True, type=Path,
                        help="One selected core gene/locus tag per line.")
    parser.add_argument("--output-prefix", required=True, type=Path)
    parser.add_argument("--min-cov", default=10, type=int)
    parser.add_argument("--include-supported", action="store_true",
                        help="Also write supported sites to site_calls.csv.")
    args = parser.parse_args()

    selected_genes = read_gene_list(args.gene_list)
    ref = read_genbank_origin(args.reference_gb)
    position_map, gene_info = selected_cds_position_map(args.reference_gb, selected_genes)
    selected_positions = set(position_map)
    print(f"[INFO] selected genes in list: {len(selected_genes)}", flush=True)
    print(f"[INFO] selected CDS features in GenBank: {len(gene_info)}", flush=True)
    print(f"[INFO] selected CDS positions in GenBank: {len(selected_positions)}", flush=True)

    site_fields = [
        "sample",
        "position",
        "locus_tag",
        "gene",
        "old_locus_tag",
        "product",
        "ref_base",
        "assembly_base",
        "cns_base",
        "cns_coverage",
        "call",
    ]
    site_path = args.output_prefix.with_suffix(".site_calls.csv")
    gene_summary_path = args.output_prefix.with_suffix(".gene_summary.csv")

    gene_counts: Dict[str, Counter] = defaultdict(Counter)
    gene_sample_counts: Dict[str, Dict[str, Counter]] = defaultdict(lambda: defaultdict(Counter))
    gene_meta: Dict[str, dict] = {}

    with site_path.open("w", newline="") as site_handle:
        writer = csv.DictWriter(site_handle, fieldnames=site_fields)
        writer.writeheader()

        for sample, seq in iter_fasta(args.alignment_fasta):
            print(f"[INFO] {sample}", flush=True)
            n = min(len(seq), len(ref))
            mutations: Dict[int, str] = {}
            for pos in selected_positions:
                if pos < 1 or pos > n:
                    continue
                ref_base = ref[pos - 1]
                aln_base = seq[pos - 1]
                if ref_base in VALID_BASES and aln_base in VALID_BASES and ref_base != aln_base:
                    mutations[pos] = aln_base

            cns_path = find_cns(sample, args.cns_dir)
            if cns_path is None:
                cns = {}
            else:
                cns = read_cns_for_positions(cns_path, set(mutations), args.min_cov)

            for pos, aln_base in sorted(mutations.items()):
                ref_base = ref[pos - 1]
                call, cns_base, cov = classify_site(pos, aln_base, cns)

                for annot in position_map[pos]:
                    locus_tag = annot["locus_tag"] or annot["gene"] or annot["old_locus_tag"]
                    gene_meta[locus_tag] = annot
                    gene_counts[locus_tag]["assembly_mutations_vs_reference"] += 1
                    gene_counts[locus_tag][call] += 1
                    gene_sample_counts[locus_tag][sample]["assembly_mutations_vs_reference"] += 1
                    gene_sample_counts[locus_tag][sample][call] += 1
                    if call == "supported" and not args.include_supported:
                        continue
                    writer.writerow({
                        "sample": sample,
                        "position": pos,
                        "locus_tag": annot["locus_tag"],
                        "gene": annot["gene"],
                        "old_locus_tag": annot["old_locus_tag"],
                        "product": annot["product"],
                        "ref_base": ref_base,
                        "assembly_base": aln_base,
                        "cns_base": cns_base,
                        "cns_coverage": cov,
                        "call": call,
                    })

    summary_fields = [
        "locus_tag",
        "gene",
        "old_locus_tag",
        "product",
        "location",
        "cds_length",
        "assembly_mutations_vs_reference",
        "supported",
        "inconsistent",
        "missing",
        "problematic",
        "samples_with_assembly_mutations",
        "problematic_rate",
        "inconsistent_rate",
        "missing_rate",
        "pooled_problematic_rate",
        "pooled_inconsistent_rate",
        "pooled_missing_rate",
    ]
    with gene_summary_path.open("w", newline="") as summary_handle:
        writer = csv.DictWriter(summary_handle, fieldnames=summary_fields)
        writer.writeheader()
        for locus_tag, counts in sorted(
            gene_counts.items(),
            key=lambda item: (-(item[1]["inconsistent"] + item[1]["missing"]), item[0]),
        ):
            meta = gene_meta[locus_tag]
            total = counts["assembly_mutations_vs_reference"]
            problematic = counts["inconsistent"] + counts["missing"]
            sample_rates = []
            for sample_counts in gene_sample_counts[locus_tag].values():
                sample_total = sample_counts["assembly_mutations_vs_reference"]
                if sample_total == 0:
                    continue
                sample_problematic = sample_counts["inconsistent"] + sample_counts["missing"]
                sample_rates.append((
                    sample_problematic / sample_total,
                    sample_counts["inconsistent"] / sample_total,
                    sample_counts["missing"] / sample_total,
                ))
            if sample_rates:
                mean_problematic_rate = sum(x[0] for x in sample_rates) / len(sample_rates)
                mean_inconsistent_rate = sum(x[1] for x in sample_rates) / len(sample_rates)
                mean_missing_rate = sum(x[2] for x in sample_rates) / len(sample_rates)
            else:
                mean_problematic_rate = ""
                mean_inconsistent_rate = ""
                mean_missing_rate = ""
            writer.writerow({
                "locus_tag": meta["locus_tag"],
                "gene": meta["gene"],
                "old_locus_tag": meta["old_locus_tag"],
                "product": meta["product"],
                "location": meta["location"],
                "cds_length": meta["cds_length"],
                "assembly_mutations_vs_reference": total,
                "supported": counts["supported"],
                "inconsistent": counts["inconsistent"],
                "missing": counts["missing"],
                "problematic": problematic,
                "samples_with_assembly_mutations": len(sample_rates),
                "problematic_rate": f"{mean_problematic_rate:.6f}" if sample_rates else "",
                "inconsistent_rate": f"{mean_inconsistent_rate:.6f}" if sample_rates else "",
                "missing_rate": f"{mean_missing_rate:.6f}" if sample_rates else "",
                "pooled_problematic_rate": f"{problematic / total:.6f}" if total else "",
                "pooled_inconsistent_rate": f"{counts['inconsistent'] / total:.6f}" if total else "",
                "pooled_missing_rate": f"{counts['missing'] / total:.6f}" if total else "",
            })

    print(f"[INFO] wrote {site_path}", flush=True)
    print(f"[INFO] wrote {gene_summary_path}", flush=True)


if __name__ == "__main__":
    main()
