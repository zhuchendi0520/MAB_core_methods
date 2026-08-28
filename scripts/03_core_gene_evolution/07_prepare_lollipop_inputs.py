#!/usr/bin/env python3
from pathlib import Path
import argparse
import glob
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Prepare pre-expansion mutation lollipop inputs")
parser.add_argument("input_dir", type=Path)
parser.add_argument("genbank", type=Path)
parser.add_argument("uniprot_features", type=Path)
parser.add_argument("output_dir", type=Path)
args = parser.parse_args()
BASE = args.input_dir
GB = args.genbank
UNIPROT = args.uniprot_features
OUT = args.output_dir

TARGETS = [
    "MAB_3029", "MAB_0674", "MAB_1638", "MAB_0173", "MAB_4596c",
    "MAB_3404c", "MAB_1678c", "MAB_1915", "MAB_1499", "MAB_0401",
    "MAB_4595c", "MAB_1060", "MAB_0673", "MAB_2928", "MAB_4470c",
    "MAB_0101", "MAB_4294", "MAB_1076", "MAB_2819", "MAB_4075",
    "MAB_4559", "MAB_0114", "MAB_1218", "MAB_0019", "MAB_3045c",
    "MAB_1635", "MAB_1140", "MAB_3675", "MAB_2689",
]

OUT.mkdir(parents=True, exist_ok=True)
frames = [pd.read_csv(p) for p in glob.glob(str(BASE / "new_outlier_cleaned/DCC*.raw_stem_post_events.csv"))]
raw = pd.concat(frames, ignore_index=True)
raw = raw[(raw["stage"] == "pre") & (raw["region"] == "mutation") & raw["gene"].isin(TARGETS)].copy()

raw["ATCC19977_GENE"] = raw["gene"]
raw["AA_OR_CDS_POSITION"] = raw["codon"]
raw["ATCC19977_LOC"] = raw["genome_pos"]
raw["EVENT_REF_ATCC_ORIENTATION"] = raw["ref"]
raw["EVENT_ALT_ATCC_ORIENTATION"] = raw["alt"]
raw["SNP_TYPE"] = raw["mutation_type"].map({"nonsynonymous": "NSY", "synonymous": "SYN"}).fillna("SYN")
raw["CONSEQUENCE"] = raw.apply(
    lambda x: ("Nonsynonymous" if x["SNP_TYPE"] == "NSY" else "Synonymous")
    + f"-{x['aa_ref']}/{x['aa_ref']}-{x['aa_alt']}/{x['aa_alt']}", axis=1
)
raw["PRODUCT"] = ""
raw["DCC"] = raw["DCC"].astype(str)

mutation_cols = [
    "DCC", "ATCC19977_GENE", "AA_OR_CDS_POSITION", "ATCC19977_LOC",
    "EVENT_REF_ATCC_ORIENTATION", "EVENT_ALT_ATCC_ORIENTATION", "SNP_TYPE",
    "CONSEQUENCE", "PRODUCT", "target_node", "node",
]
raw[mutation_cols].to_csv(OUT / "pre_29_genes_mutations.csv", index=False)

metadata = {}
for record in SeqIO.parse(str(GB), "genbank"):
    for feature in record.features:
        if feature.type != "CDS":
            continue
        tags = feature.qualifiers.get("locus_tag", []) + feature.qualifiers.get("gene", [])
        gene = next((x for x in tags if x in TARGETS), None)
        if gene:
            aa = feature.qualifiers.get("translation", [""])[0]
            metadata[gene] = {
                "gene": gene,
                "protein_length_aa": len(aa) if aa else int(len(feature.location) / 3),
                "product": feature.qualifiers.get("product", [""])[0],
            }

pd.DataFrame([metadata[g] for g in TARGETS]).to_csv(OUT / "pre_29_genes_metadata.csv", index=False)

if UNIPROT.exists():
    pd.read_csv(UNIPROT, sep="\t").to_csv(
        OUT / "pre_29_genes_uniprot_features.tsv", sep="\t", index=False
    )
else:
    pd.DataFrame(columns=[
        "Entry", "Entry Name", "Gene Names", "Protein names", "Length", "Domain [FT]",
        "Region", "Repeat", "Coiled coil", "Transmembrane", "Topological domain",
        "Motif", "Binding site", "Active site",
    ]).to_csv(OUT / "pre_29_genes_uniprot_features.tsv", sep="\t", index=False)

print(f"pre mutation rows: {len(raw)}")
print(f"genes: {raw['ATCC19977_GENE'].nunique()}")
print(OUT)
