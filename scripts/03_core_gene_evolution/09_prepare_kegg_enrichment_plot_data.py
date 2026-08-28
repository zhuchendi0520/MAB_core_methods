#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path


KEGG_NAMES = {
    "ko01100": "Metabolic pathways",
    "ko00550": "Peptidoglycan biosynthesis",
    "ko04112": "Cell cycle - Caulobacter",
    "ko00190": "Oxidative phosphorylation",
    "ko03050": "Proteasome",
    "ko01053": "Biosynthesis of siderophore group nonribosomal peptides",
    "ko02020": "Two-component system",
    "ko00230": "Purine metabolism",
    "ko00240": "Pyrimidine metabolism",
    "ko00984": "Steroid degradation",
    "ko00860": "Porphyrin metabolism",
    "ko04922": "Glucagon signaling pathway",
    "ko00500": "Starch and sucrose metabolism",
    "ko02010": "ABC transporters",
    "ko03420": "Nucleotide excision repair",
    "M00144": "NADH:quinone oxidoreductase, prokaryotes",
    "M00342": "Cytochrome bc1 complex",
    "M00019": "Valine/isoleucine biosynthesis",
    "M00046": "Pyrimidine degradation",
    "M00082": "Fatty acid biosynthesis",
    "M00570": "Isoleucine biosynthesis",
    "M00089": "Triacylglycerol biosynthesis",
    "M00121": "Heme biosynthesis",
    "M00307": "Pyruvate oxidation",
    "M00302": "Dicarboxylate-hydroxybutyrate cycle",
    "M00026": "Histidine biosynthesis",
    "M00210": "ABC-2 type transport system",
    "M00669": "Multidrug resistance efflux pump",
}


def parse_args():
    p = argparse.ArgumentParser(description="Prepare KEGG enrichment dot plot data.")
    p.add_argument("--enrichment", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--top-n", type=int, default=10)
    p.add_argument("--min-count", type=int, default=3)
    return p.parse_args()


def as_float(x, default=math.nan):
    try:
        return float(x)
    except Exception:
        return default


def clean_term(row):
    term = row["term"]
    term_type = row["term_type"]
    if term_type == "KEGG_Pathway":
        if not term.startswith("ko"):
            return None
        return term
    if term_type == "KEGG_Module":
        return term
    return None


def main():
    args = parse_args()
    rows = []
    with open(args.enrichment, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            if row["term_type"] not in {"KEGG_Pathway", "KEGG_Module"}:
                continue
            term = clean_term(row)
            if term is None:
                continue
            count = int(row["target_genes_in_term"])
            if count < args.min_count:
                continue
            fdr = as_float(row["fdr_BH_within_term_type"], 1.0)
            pvalue = as_float(row["p_value"], 1.0)
            bg_count = int(row["background_genes_in_term"])
            target_total = int(row["target_genes_total"])
            bg_total = int(row["background_genes_total"])
            rows.append({
                "direction": row["direction"],
                "term_type": "Pathway" if row["term_type"] == "KEGG_Pathway" else "Module",
                "term": term,
                "term_name": KEGG_NAMES.get(term, term),
                "label": f"{KEGG_NAMES.get(term, term)} ({term})",
                "target_genes_in_term": count,
                "target_genes_total": target_total,
                "background_genes_in_term": bg_count,
                "background_genes_total": bg_total,
                "gene_ratio": count / target_total,
                "background_ratio": bg_count / bg_total,
                "odds_ratio": as_float(row["odds_ratio"]),
                "p_value": pvalue,
                "fdr_BH": fdr,
                "neg_log10_fdr": -math.log10(max(fdr, 1e-300)),
                "genes": row["genes"],
            })

    selected = []
    for direction in ["increased", "decreased"]:
        sub = [r for r in rows if r["direction"] == direction]
        sig = sorted([r for r in sub if r["fdr_BH"] <= 0.05], key=lambda r: (r["fdr_BH"], r["p_value"]))
        nonsig = sorted([r for r in sub if r["fdr_BH"] > 0.05], key=lambda r: (r["p_value"], r["fdr_BH"]))
        pick = (sig + nonsig)[:args.top_n]
        selected.extend(pick)

    fields = [
        "direction", "term_type", "term", "term_name", "label",
        "target_genes_in_term", "target_genes_total",
        "background_genes_in_term", "background_genes_total",
        "gene_ratio", "background_ratio", "odds_ratio",
        "p_value", "fdr_BH", "neg_log10_fdr", "genes",
    ]
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in selected:
            w.writerow(row)

    print(f"Wrote {len(selected)} rows to {out}")


if __name__ == "__main__":
    main()
