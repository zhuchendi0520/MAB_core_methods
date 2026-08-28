#!/usr/bin/env python3

import sys
import os
import pandas as pd
from Bio import AlignIO

if len(sys.argv) not in (2, 3):
    print("Usage: python 02_evaluate_core_gene_alignments.py ALIGNMENT_DIR [OUTPUT_CSV]")
    sys.exit(1)

dir_path = sys.argv[1]
output_csv = sys.argv[2] if len(sys.argv) == 3 else "aligned_gene_quality.csv"

results = []

# Process Panaroo per-gene alignments.
for f in os.listdir(dir_path):
    if not f.endswith(".aln.fas"):
        continue

    path = os.path.join(dir_path, f)

    # Preserve the Panaroo cluster identifier.
    gene = f.replace(".aln.fas", "")

    try:
        aln = AlignIO.read(path, "fasta")
    except Exception as e:
        print(f"Skip {f}, parse error: {e}")
        continue

    aln_len = aln.get_alignment_length()
    n_seq = len(aln)

    gap_cols = 0

    # Count columns dominated by gaps.
    for i in range(aln_len):
        col = aln[:, i]
        gap_ratio = col.count("-") / n_seq

        # A column is gap dominated when more than half of sequences contain gaps.
        if gap_ratio > 0.5:
            gap_cols += 1

    overall_gap_ratio = gap_cols / aln_len

    # Classify genes by the fraction of gap-dominated columns.
    if overall_gap_ratio <= 0.1:
        status = "high_quality_core"
    elif overall_gap_ratio <= 0.3:
        status = "medium_core"
    else:
        status = "bad_core"

    results.append([
        gene,
        overall_gap_ratio,
        status,
        n_seq,
        aln_len
    ])

df = pd.DataFrame(
    results,
    columns=[
        "gene",
        "gap_ratio",
        "status",
        "num_sequences",
        "alignment_length"
    ]
)

df.to_csv(output_csv, index=False)

print("Done")
print(f"Total genes processed: {len(df)}")
print(f"Output: {output_csv}")
