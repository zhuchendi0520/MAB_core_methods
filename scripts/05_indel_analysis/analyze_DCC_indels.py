#!/usr/bin/env python3
"""DCC fixed and phylogeny-aware convergent indel analysis."""

from __future__ import annotations

import argparse
import bisect
import csv
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo, SeqIO


DCCS = [f"DCC{i}" for i in range(1, 8)]
INFERRED_GROUPS = {
    "SRR15713784": "DCC1",
    "SRR32604902": "DCC1",
    "SRR32604920": "DCC4",
    "SRR18969534": "DCC6",
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--indel-dir", required=True, type=Path)
    p.add_argument("--groups", required=True, type=Path)
    p.add_argument("--tree", required=True, type=Path)
    p.add_argument("--genbank", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--fixed-threshold", type=float, default=0.95)
    p.add_argument("--other-max", type=float, default=0.05)
    p.add_argument("--top-genes", type=int, default=40)
    return p.parse_args()


def event_interval(pos, allele):
    seq = allele[1:]
    if allele.startswith("-"):
        return pos + 1, pos + len(seq), "deletion", len(seq)
    return pos, pos, "insertion", len(seq)


def load_features(path):
    record = SeqIO.read(path, "genbank")
    rows = []
    for f in record.features:
        if f.type not in {"CDS", "rRNA", "tRNA", "ncRNA"}:
            continue
        q = f.qualifiers
        rows.append({
            "start": int(f.location.start) + 1,
            "end": int(f.location.end),
            "strand": f.location.strand,
            "feature_type": f.type,
            "locus_tag": q.get("locus_tag", [""])[0],
            "gene": q.get("gene", [""])[0],
            "product": q.get("product", [""])[0],
        })
    rows.sort(key=lambda x: x["start"])
    return record.id, len(record.seq), rows, [x["start"] for x in rows]


def annotate_interval(start, end, features, starts):
    i = max(0, bisect.bisect_right(starts, end) - 1)
    hits = []
    while i >= 0 and features[i]["end"] >= start:
        if features[i]["start"] <= end:
            hits.append(features[i])
        i -= 1
    if hits:
        hits.sort(key=lambda x: x["start"])
        return hits
    j = bisect.bisect_left(starts, start)
    nearest = []
    if j > 0:
        nearest.append(features[j - 1])
    if j < len(features):
        nearest.append(features[j])
    if not nearest:
        return []
    d = lambda x: min(abs(start - x["end"]), abs(x["start"] - end))
    return [min(nearest, key=d) | {"feature_type": "intergenic_nearest"}]


def read_groups(path):
    df = pd.read_csv(path)
    return dict(zip(df.iloc[:, 0].astype(str), df.iloc[:, 1].astype(str)))


def read_events(indel_dir, groups):
    event_samples = defaultdict(set)
    sample_events = {}
    assignments = []
    for path in sorted(indel_dir.glob("*.indel_filter")):
        sample = path.name.removesuffix(".indel_filter")
        group = groups.get(sample, INFERRED_GROUPS.get(sample, "Unassigned"))
        events = set()
        with path.open() as fh:
            for line in fh:
                parts = line.split()
                if len(parts) < 3 or not parts[0].isdigit():
                    continue
                pos, ref, allele = int(parts[0]), parts[1].upper(), parts[2].upper()
                if not allele.startswith(("+", "-")):
                    continue
                event = (pos, ref, allele)
                events.add(event)
                event_samples[event].add(sample)
        sample_events[sample] = events
        assignments.append({"sample": sample, "group": group,
                            "group_source": "table" if sample in groups else "tree_inferred",
                            "n_indels": len(events)})
    return event_samples, sample_events, pd.DataFrame(assignments)


def annotate_event(event, features, starts):
    pos, ref, allele = event
    s, e, typ, length = event_interval(pos, allele)
    hits = annotate_interval(s, e, features, starts)
    if not hits:
        hits = [{"locus_tag": "", "gene": "", "product": "", "feature_type": "intergenic",
                 "start": np.nan, "end": np.nan, "strand": np.nan}]
    return s, e, typ, length, hits


def event_prevalence(event_samples, assignments, features, starts, threshold, other_max):
    sample_group = dict(zip(assignments["sample"], assignments["group"]))
    denominators = assignments.query("group in @DCCS").groupby("group")["sample"].nunique().to_dict()
    rows = []
    for event, samples in event_samples.items():
        c = Counter(sample_group.get(s) for s in samples)
        s, e, typ, length, hits = annotate_event(event, features, starts)
        fixed = [g for g in DCCS if denominators.get(g, 0) and c[g] / denominators[g] >= threshold]
        if not fixed:
            continue
        base = {"position": event[0], "ref": event[1], "allele": event[2], "type": typ,
                "event_length": length, "affected_start": s, "affected_end": e,
                "fixed_DCCs": "|".join(fixed), "n_fixed_DCCs": len(fixed)}
        for g in DCCS:
            base[f"{g}_n"] = c[g]
            base[f"{g}_pct"] = c[g] / denominators.get(g, np.nan)
        base["classification"] = "shared_fixed" if len(fixed) > 1 else "fixed_unique"
        if len(fixed) == 1:
            target = fixed[0]
            base["strict_specific"] = all(base[f"{g}_pct"] < other_max for g in DCCS if g != target)
        else:
            base["strict_specific"] = False
        for hit in hits:
            rows.append(base | hit)
    return pd.DataFrame(rows), denominators


def gene_prevalence(sample_events, assignments, features, starts, threshold):
    sample_group = dict(zip(assignments["sample"], assignments["group"]))
    den = assignments.query("group in @DCCS").groupby("group")["sample"].nunique().to_dict()
    gene_samples = defaultdict(set)
    gene_meta = {}
    for sample, events in sample_events.items():
        if sample_group.get(sample) not in DCCS:
            continue
        seen = set()
        for event in events:
            s, e, _, _, hits = annotate_event(event, features, starts)
            for h in hits:
                tag = h["locus_tag"] or f"intergenic:{s}-{e}"
                seen.add(tag)
                gene_meta[tag] = h
        for tag in seen:
            gene_samples[tag].add(sample)
    rows = []
    for tag, samples in gene_samples.items():
        c = Counter(sample_group[s] for s in samples)
        row = {"locus_tag": tag, "gene": gene_meta[tag].get("gene", ""),
               "product": gene_meta[tag].get("product", ""),
               "feature_type": gene_meta[tag].get("feature_type", "")}
        fixed = []
        for g in DCCS:
            row[f"{g}_n"] = c[g]
            row[f"{g}_pct"] = c[g] / den.get(g, np.nan)
            if row[f"{g}_pct"] >= threshold:
                fixed.append(g)
        row["fixed_DCCs"] = "|".join(fixed)
        row["n_fixed_DCCs"] = len(fixed)
        rows.append(row)
    return pd.DataFrame(rows)


def prepare_tree(tree_path, valid_samples):
    tree = Phylo.read(tree_path, "newick")
    parent = {}
    children = {}
    tips = {}
    for clade in tree.find_clades(order="preorder"):
        children[clade] = list(clade.clades)
        for child in clade.clades:
            parent[child] = clade
        if clade.is_terminal() and clade.name in valid_samples:
            tips[clade.name] = clade
    return tree, parent, children, tips


def maximal_carrier_clades(carriers, parent, children, tips):
    active = {tips[s] for s in carriers if s in tips}
    queue = list({parent[x] for x in active if x in parent})
    queued = set(queue)
    while queue:
        node = queue.pop()
        queued.discard(node)
        kids = children[node]
        if kids and all(k in active for k in kids):
            active.difference_update(kids)
            active.add(node)
            if node in parent and parent[node] not in queued:
                queue.append(parent[node]); queued.add(parent[node])
    return active


def descendant_groups(node, children, tip_group, cache):
    if node in cache:
        return cache[node]
    if not children[node]:
        cache[node] = {tip_group.get(node.name, "Unassigned")}
    else:
        cache[node] = set().union(*(descendant_groups(k, children, tip_group, cache) for k in children[node]))
    cache[node].discard("Unassigned")
    return cache[node]


def convergence(event_samples, assignments, tree_path, features, starts):
    tip_group = dict(zip(assignments["sample"], assignments["group"]))
    tree, parent, children, tips = prepare_tree(tree_path, set(tip_group))
    cache = {}
    descendant_groups(tree.root, children, tip_group, cache)
    event_rows, gene_acc = [], defaultdict(lambda: {"events": set(), "samples": set(), "gains": Counter()})
    for event, samples in event_samples.items():
        dcc_samples = {s for s in samples if tip_group.get(s) in DCCS and s in tips}
        if not dcc_samples:
            continue
        clusters = maximal_carrier_clades(dcc_samples, parent, children, tips)
        gains = Counter()
        for node in clusters:
            gs = cache.get(node) or descendant_groups(node, children, tip_group, cache)
            if len(gs) == 1:
                g = next(iter(gs))
                if g in DCCS:
                    gains[g] += 1
        s, e, typ, length, hits = annotate_event(event, features, starts)
        base = {"position": event[0], "ref": event[1], "allele": event[2], "type": typ,
                "event_length": length, "affected_start": s, "affected_end": e,
                "total_min_independent_gains": sum(gains.values()), "n_carriers": len(dcc_samples)}
        for g in DCCS:
            base[f"{g}_gains"] = gains[g]
        for hit in hits:
            event_rows.append(base | hit)
            tag = hit["locus_tag"] or f"intergenic:{s}-{e}"
            acc = gene_acc[tag]
            acc["meta"] = hit
            acc["events"].add(event)
            acc["samples"].update(dcc_samples)
            acc["gains"].update(gains)
    gene_rows = []
    for tag, acc in gene_acc.items():
        row = {"locus_tag": tag, "gene": acc["meta"].get("gene", ""),
               "product": acc["meta"].get("product", ""),
               "feature_type": acc["meta"].get("feature_type", ""),
               "distinct_indel_events": len(acc["events"]), "carrier_strains": len(acc["samples"]),
               "total_min_independent_gains": sum(acc["gains"].values())}
        for g in DCCS:
            row[f"{g}_gains"] = acc["gains"][g]
        row["DCCs_with_parallel_gains"] = sum(acc["gains"][g] >= 2 for g in DCCS)
        gene_rows.append(row)
    return pd.DataFrame(event_rows), pd.DataFrame(gene_rows)


def plot_fixed_heatmap(df, out):
    if df.empty:
        return
    x = df.drop_duplicates(["position", "ref", "allele"]).copy()
    x = x.sort_values(["n_fixed_DCCs", "type", "position"], ascending=[False, True, True]).head(120)
    mat = np.array([[r[f"{g}_pct"] for g in DCCS] for _, r in x.iterrows()])
    fig, ax = plt.subplots(figsize=(6.2, max(3.2, len(x) * .055)))
    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
    ax.set_xticks(range(7), DCCS)
    ax.set_yticks([])
    ax.set_xlabel("DCC")
    ax.set_ylabel("Fixed indel events")
    cb = fig.colorbar(im, ax=ax, pad=.02); cb.set_label("Carriage proportion")
    fig.tight_layout(); fig.savefig(out, dpi=400); plt.close(fig)


def plot_parallel(df, out, top_n):
    if df.empty:
        return
    x = df.sort_values(["total_min_independent_gains", "DCCs_with_parallel_gains"], ascending=False).head(top_n)
    labels = [r.gene if isinstance(r.gene, str) and r.gene else r.locus_tag for _, r in x.iterrows()]
    mat = x[[f"{g}_gains" for g in DCCS]].to_numpy()
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(8.5, max(4.5, len(x) * .22)),
                                  gridspec_kw={"width_ratios": [7, 1.25]}, sharey=True)
    im = ax.imshow(mat, aspect="auto", cmap="Blues", vmin=0)
    ax.set_xticks(range(7), DCCS, rotation=45, ha="right")
    ax.set_yticks(range(len(labels)), labels, fontsize=8)
    ax.set_xlabel("Minimum independent gains")
    cb = fig.colorbar(im, ax=ax, pad=.02); cb.set_label("Gain count")
    ax2.barh(range(len(x)), x["total_min_independent_gains"], color="#E59672")
    ax2.set_xlabel("Total")
    ax2.tick_params(axis="y", left=False, labelleft=False)
    fig.tight_layout(); fig.savefig(out, dpi=400); plt.close(fig)


def main():
    a = parse_args(); a.outdir.mkdir(parents=True, exist_ok=True)
    ref_id, ref_len, features, starts = load_features(a.genbank)
    groups = read_groups(a.groups)
    event_samples, sample_events, assignments = read_events(a.indel_dir, groups)
    assignments.to_csv(a.outdir / "sample_group_assignments.csv", index=False)

    fixed, den = event_prevalence(event_samples, assignments, features, starts,
                                  a.fixed_threshold, a.other_max)
    fixed.to_csv(a.outdir / "DCC_fixed_indels_95pct.csv", index=False)
    fixed.query("classification == 'shared_fixed'").to_csv(a.outdir / "DCC_shared_fixed_indels.csv", index=False)
    fixed.query("classification == 'fixed_unique'").to_csv(a.outdir / "DCC_fixed_unique_indels.csv", index=False)
    fixed.query("strict_specific == True").to_csv(a.outdir / "DCC_strict_specific_indels.csv", index=False)

    genes = gene_prevalence(sample_events, assignments, features, starts, a.fixed_threshold)
    genes.to_csv(a.outdir / "DCC_gene_level_indel_prevalence.csv", index=False)
    genes.query("n_fixed_DCCs >= 2").to_csv(a.outdir / "DCC_shared_fixed_indel_genes.csv", index=False)

    parallel_events, parallel_genes = convergence(event_samples, assignments, a.tree, features, starts)
    parallel_events.sort_values("total_min_independent_gains", ascending=False).to_csv(
        a.outdir / "DCC_parallel_exact_indel_events.csv", index=False)
    parallel_genes.sort_values("total_min_independent_gains", ascending=False).to_csv(
        a.outdir / "DCC_parallel_indel_gene_summary.csv", index=False)

    summary = [
        ["reference", ref_id], ["reference_length", ref_len], ["samples", len(assignments)],
        ["unique_exact_indels", len(event_samples)], ["fixed_indel_rows", len(fixed)],
        ["shared_fixed_exact_indels", fixed.query("classification == 'shared_fixed'")[["position", "ref", "allele"]].drop_duplicates().shape[0]],
        ["strict_DCC_specific_exact_indels", fixed.query("strict_specific == True")[["position", "ref", "allele"]].drop_duplicates().shape[0]],
        ["shared_fixed_indel_genes", int((genes.n_fixed_DCCs >= 2).sum())],
        ["parallel_exact_events_ge2_gains", int((parallel_events.total_min_independent_gains >= 2).sum())],
        ["parallel_genes_ge2_gains", int((parallel_genes.total_min_independent_gains >= 2).sum())],
    ] + [[f"{g}_samples", den.get(g, 0)] for g in DCCS]
    pd.DataFrame(summary, columns=["metric", "value"]).to_csv(a.outdir / "analysis_summary.csv", index=False)
    plot_fixed_heatmap(fixed, a.outdir / "DCC_fixed_indel_carriage_heatmap.png")
    plot_parallel(parallel_genes, a.outdir / "DCC_parallel_indel_gene_heatmap.png", a.top_genes)


if __name__ == "__main__":
    main()
