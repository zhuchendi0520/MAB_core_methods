#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
from pathlib import Path

from scipy.stats import fisher_exact


COG_NAMES = {
    "J": "Translation, ribosomal structure and biogenesis",
    "A": "RNA processing and modification",
    "K": "Transcription",
    "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics",
    "D": "Cell cycle control, cell division, chromosome partitioning",
    "Y": "Nuclear structure",
    "V": "Defense mechanisms",
    "T": "Signal transduction mechanisms",
    "M": "Cell wall/membrane/envelope biogenesis",
    "N": "Cell motility",
    "Z": "Cytoskeleton",
    "W": "Extracellular structures",
    "U": "Intracellular trafficking, secretion, and vesicular transport",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "X": "Mobilome: prophages, transposons",
    "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism",
    "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism",
    "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism",
    "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only",
    "S": "Function unknown",
}


TERM_COLUMNS = [
    "COG_category",
    "GOs",
    "EC",
    "KEGG_ko",
    "KEGG_Pathway",
    "KEGG_Module",
    "KEGG_Reaction",
    "BRITE",
    "KEGG_TC",
    "CAZy",
    "PFAMs",
]


def parse_args():
    p = argparse.ArgumentParser(description="eggNOG term enrichment for increased/decreased r/m genes.")
    p.add_argument("--sig", required=True, help="Significant r/m stats CSV.")
    p.add_argument("--background", required=True, help="All-gene background CSV.")
    p.add_argument("--eggnog", required=True, help="eggNOG emapper annotations file.")
    p.add_argument("--outdir", required=True, help="Output directory.")
    p.add_argument("--min-bg", type=int, default=3, help="Minimum background genes for a term.")
    return p.parse_args()


def read_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def read_sig(path):
    rows = []
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        header = next(reader)
        for values in reader:
            row = {h: values[i] if i < len(values) else "" for i, h in enumerate(header)}
            if len(values) >= 17:
                row["rm_delta"] = values[-2]
                row["direction"] = values[-1].strip().lower()
            rows.append(row)
    return rows


def gene_from_row(row):
    return (row.get("ATCC19977") or row.get("gene") or row.get("#query") or row.get("query") or "").strip()


def hmm_call_from_row(row):
    """Return the HMM call without depending on parenthetical citation text."""
    for key, value in row.items():
        if key.strip().lower().startswith("hmm call"):
            return value
    return ""


def read_eggnog(path):
    header = None
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for line in f:
            if line.startswith("##"):
                continue
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#query"):
                header = line.lstrip("#").split("\t")
                continue
            if header is None:
                continue
            values = line.split("\t")
            row = {h: values[i] if i < len(values) else "" for i, h in enumerate(header)}
            rows.append(row)
    return rows


def split_terms(value, col):
    value = (value or "").strip()
    if not value or value == "-":
        return []
    if col == "COG_category":
        return [c for c in value if c and c != "-"]
    return [x.strip() for x in value.split(",") if x.strip() and x.strip() != "-"]


def term_label(term_type, term):
    if term_type == "COG_category":
        return f"{term}: {COG_NAMES.get(term, 'Unknown COG category')}"
    return term


def bh_adjust(pvals):
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    adj = [1.0] * n
    prev = 1.0
    for rank_from_end, idx in enumerate(reversed(order), start=1):
        rank = n - rank_from_end + 1
        val = pvals[idx] * n / rank
        prev = min(prev, val)
        adj[idx] = min(prev, 1.0)
    return adj


def enrich_for_type(direction, target_genes, universe_genes, gene_terms, term_type, min_bg):
    term_bg = defaultdict(set)
    term_target = defaultdict(set)
    universe = set(universe_genes)
    target = set(target_genes) & universe
    for gene in universe:
        for term in gene_terms.get(gene, []):
            term_bg[term].add(gene)
            if gene in target:
                term_target[term].add(gene)

    rows = []
    for term, bg_genes in term_bg.items():
        if len(bg_genes) < min_bg:
            continue
        a = len(term_target.get(term, set()))
        if a == 0:
            continue
        b = len(target) - a
        c = len(bg_genes) - a
        d = len(universe) - len(target) - c
        odds, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({
            "direction": direction,
            "term_type": term_type,
            "term": term,
            "term_label": term_label(term_type, term),
            "target_genes_in_term": a,
            "target_genes_total": len(target),
            "background_genes_in_term": len(bg_genes),
            "background_genes_total": len(universe),
            "odds_ratio": odds,
            "p_value": p,
            "genes": ";".join(sorted(term_target.get(term, set()))),
        })
    adj = bh_adjust([r["p_value"] for r in rows]) if rows else []
    for row, fdr in zip(rows, adj):
        row["fdr_BH_within_term_type"] = fdr
    rows.sort(key=lambda r: (r["fdr_BH_within_term_type"], r["p_value"], -r["target_genes_in_term"]))
    return rows


def write_csv(path, rows, fields):
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sig_rows = read_sig(args.sig)
    bg_rows = read_csv(args.background)
    eggnog_rows = read_eggnog(args.eggnog)

    background_genes_raw = {gene_from_row(r) for r in bg_rows if gene_from_row(r)}
    eggnog_by_gene = {gene_from_row(r): r for r in eggnog_rows if gene_from_row(r)}
    universe_genes = sorted(background_genes_raw & set(eggnog_by_gene))

    direction_genes = defaultdict(set)
    for row in sig_rows:
        gene = gene_from_row(row)
        direction = row.get("direction", "").strip().lower()
        if gene and direction in {"increased", "decreased"}:
            direction_genes[direction].add(gene)

    all_rows = []
    for term_type in TERM_COLUMNS:
        gene_terms = {}
        for gene in universe_genes:
            gene_terms[gene] = split_terms(eggnog_by_gene[gene].get(term_type, ""), term_type)
        for direction in ["increased", "decreased"]:
            all_rows.extend(enrich_for_type(
                direction,
                direction_genes[direction],
                universe_genes,
                gene_terms,
                term_type,
                args.min_bg,
            ))

    fields = [
        "direction", "term_type", "term", "term_label",
        "target_genes_in_term", "target_genes_total",
        "background_genes_in_term", "background_genes_total",
        "odds_ratio", "p_value", "fdr_BH_within_term_type", "genes",
    ]
    write_csv(outdir / "rm_direction_eggnog_enrichment_all.csv", all_rows, fields)
    write_csv(
        outdir / "rm_direction_eggnog_enrichment_FDR0.05.csv",
        [r for r in all_rows if r["fdr_BH_within_term_type"] <= 0.05],
        fields,
    )

    sig_anno_fields = [
        "direction", "ATCC19977", "Preferred_name_input", "HMM_call",
        "Preferred_name_eggnog", "COG_category", "Description", "GOs", "EC",
        "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "BRITE", "PFAMs",
    ]
    sig_anno = []
    for row in sig_rows:
        gene = gene_from_row(row)
        if not gene:
            continue
        e = eggnog_by_gene.get(gene, {})
        sig_anno.append({
            "direction": row.get("direction", ""),
            "ATCC19977": gene,
            "Preferred_name_input": row.get("Preferred_name", ""),
            "HMM_call": hmm_call_from_row(row),
            "Preferred_name_eggnog": e.get("Preferred_name", ""),
            "COG_category": e.get("COG_category", ""),
            "Description": e.get("Description", ""),
            "GOs": e.get("GOs", ""),
            "EC": e.get("EC", ""),
            "KEGG_ko": e.get("KEGG_ko", ""),
            "KEGG_Pathway": e.get("KEGG_Pathway", ""),
            "KEGG_Module": e.get("KEGG_Module", ""),
            "BRITE": e.get("BRITE", ""),
            "PFAMs": e.get("PFAMs", ""),
        })
    write_csv(outdir / "rm_direction_significant_genes_with_eggnog.csv", sig_anno, sig_anno_fields)

    summary = []
    for direction in ["increased", "decreased"]:
        target = direction_genes[direction] & set(universe_genes)
        summary.append({
            "direction": direction,
            "target_genes_in_annotated_universe": len(target),
            "target_genes_raw": len(direction_genes[direction]),
            "annotated_background_genes": len(universe_genes),
            "raw_background_genes": len(background_genes_raw),
            "FDR0.05_terms": sum(1 for r in all_rows if r["direction"] == direction and r["fdr_BH_within_term_type"] <= 0.05),
        })
    write_csv(outdir / "rm_direction_eggnog_enrichment_summary.csv", summary, [
        "direction", "target_genes_in_annotated_universe", "target_genes_raw",
        "annotated_background_genes", "raw_background_genes", "FDR0.05_terms",
    ])

    print(f"Raw background genes: {len(background_genes_raw)}")
    print(f"Annotated universe genes: {len(universe_genes)}")
    for direction in ["increased", "decreased"]:
        print(f"{direction}: {len(direction_genes[direction] & set(universe_genes))}/{len(direction_genes[direction])} target genes in annotated universe")
    print(f"Outputs written to: {outdir}")


if __name__ == "__main__":
    main()
