#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np
from Bio import Phylo
from collections import defaultdict
import os
import io
import re

# =========================
# INPUT
# =========================
if len(sys.argv) != 4:
    print("Usage: python script.py <tree.nwk> <rtab> <prefix>")
    sys.exit(1)

tree_file = sys.argv[1]
rtab_file = sys.argv[2]
prefix = sys.argv[3]

print("Loading tree...")
tree = Phylo.read(tree_file, "newick")

# =========================
# Assign names to unnamed internal nodes.
# =========================
def assign_names_to_internal_nodes(tree):
    internal_node_count = 1
    for clade in tree.find_clades():
        if not clade.is_terminal() and clade.name is None:
            clade.name = f"Node_{internal_node_count}"
            internal_node_count += 1

assign_names_to_internal_nodes(tree)

print("Loading Rtab...")
df = pd.read_csv(rtab_file, sep="\t", index_col=0).fillna(0).astype(int)
genes = df.index

# =========================
# STORAGE
# =========================
edge_gain = defaultdict(int)
edge_loss = defaultdict(int)

node_gene_gain = defaultdict(set)
node_gene_loss = defaultdict(set)

gene_event = []

# =========================
# TREE traversal
# =========================
def preorder(node):
    yield node
    for c in node.clades:
        yield from preorder(c)

# =========================
# FITCH
# =========================
def fitch_bottom_up(node, tip_states):
    if node.is_terminal():
        node.fitch_set = {tip_states[node.name]}
        return node.fitch_set

    child_sets = [fitch_bottom_up(c, tip_states) for c in node.clades]
    inter = set.intersection(*child_sets)
    node.fitch_set = inter if inter else set.union(*child_sets)
    return node.fitch_set


def fitch_top_down(node, parent_state=None):
    if parent_state is None:
        node.state = 1 if 1 in node.fitch_set else 0
    else:
        node.state = parent_state if parent_state in node.fitch_set else next(iter(node.fitch_set))

    for c in node.clades:
        fitch_top_down(c, node.state)

# =========================
# MAIN LOOP
# =========================
print("Processing genes...")

for gene in genes:

    tip_states = df.loc[gene].to_dict()

    for node in tree.find_clades():
        node.fitch_set = None
        node.state = None

    fitch_bottom_up(tree.root, tip_states)
    fitch_top_down(tree.root)

    for node in preorder(tree.root):
        for child in node.clades:

            p = node.name
            c = child.name

            if node.state == 0 and child.state == 1:
                edge_gain[(p, c)] += 1
                node_gene_gain[c].add(gene)
                gene_event.append([p, c, gene, "gain"])

            if node.state == 1 and child.state == 0:
                edge_loss[(p, c)] += 1
                node_gene_loss[c].add(gene)
                gene_event.append([p, c, gene, "loss"])

# =========================
# OUTPUT 1: EDGE TABLE
# =========================
print("Writing edge table...")

edge_rows = []
edges = set(edge_gain.keys()) | set(edge_loss.keys())

for (p, c) in edges:

    g = edge_gain[(p, c)]
    l = edge_loss[(p, c)]

    ratio = g / (l + 1)
    logb = np.log((g + 1) / (l + 1))

    edge_rows.append({
        "parent": p,
        "child": c,
        "gain": g,
        "loss": l,
        "gain_loss_ratio": ratio,
        "log_bias": logb
    })

edge_df = pd.DataFrame(edge_rows)
edge_df.to_csv(f"{prefix}_edge_gain_loss.tsv", sep="\t", index=False)

# =========================
# OUTPUT 2: NEWICK TREE
# =========================
print("Writing annotated newick tree...")

for node in preorder(tree.root):
    for child in node.clades:

        p = node.name
        c = child.name

        g = edge_gain[(p, c)]
        l = edge_loss[(p, c)]

        ratio = g / (l + 1)
        logb = np.log((g + 1) / (l + 1))

        child.comment = f"&gain={g},loss={l},ratio={ratio:.4f},logb={logb:.4f}"

out_tree = f"{prefix}_annotated.tree"
Phylo.write(tree, out_tree, "newick")

# =========================
# OUTPUT 3: NODE SUMMARY
# =========================
print("Writing node summary...")

node_rows = []
for n in set(node_gene_gain.keys()) | set(node_gene_loss.keys()):
    node_rows.append({
        "node": n,
        "gained_genes": ";".join(sorted(node_gene_gain[n])),
        "lost_genes": ";".join(sorted(node_gene_loss[n]))
    })

pd.DataFrame(node_rows).to_csv(
    f"{prefix}_node_gene_summary.tsv", sep="\t", index=False
)

# =========================
# OUTPUT 4: GENE EVENTS
# =========================
print("Writing gene events...")

pd.DataFrame(
    gene_event,
    columns=["parent", "child", "gene", "event"]
).to_csv(
    f"{prefix}_gene_event.tsv", sep="\t", index=False
)

# =========================
# OUTPUT 5: NEXUS tree.
# =========================
print("Writing NEX tree (SAFE VERSION)...")

handle = io.StringIO()
Phylo.write(tree, handle, "newick")
newick_str = handle.getvalue().strip()

def annotate_node(newick, node, annotation):
    pattern = rf'(?<![A-Za-z0-9_])({re.escape(node)})(?=[:\),])'
    return re.sub(pattern, rf'\1{annotation}', newick)

for _, row in edge_df.iterrows():

    node = row["child"]
    g = row["gain"]
    l = row["loss"]
    ratio = row["gain_loss_ratio"]
    logb = row["log_bias"]

    annotation = f'[&label="{node}",gain={g},loss={l},ratio={ratio:.4f},logb={logb:.4f}]'

    newick_str = annotate_node(newick_str, node, annotation)

with open(f"{prefix}_annotated.nex", "w") as f:
    f.write("#NEXUS\n")
    f.write("Begin trees;\n")
    f.write(f"Tree tree1 = {newick_str}\n")
    f.write("End;\n")

# =========================
# DONE
# =========================
print("DONE")
print("1 edge:", f"{prefix}_edge_gain_loss.tsv")
print("2 tree:", out_tree)
print("3 node:", f"{prefix}_node_gene_summary.tsv")
print("4 gene:", f"{prefix}_gene_event.tsv")
print("5 nex :", f"{prefix}_annotated.nex")
