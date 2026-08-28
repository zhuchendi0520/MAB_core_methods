#!/usr/bin/env python3
"""Count independent indel gains inside target-node-defined DCC subtrees."""

from collections import defaultdict
import argparse
from pathlib import Path

import pandas as pd
from Bio import Phylo

from analyze_DCC_indels import annotate_event, load_features, read_events

def read_targets(path):
    targets = {}
    with path.open() as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 2 and p[0].startswith('DCC'):
                targets[p[0]] = p[1]
    return {dcc: targets[dcc] for dcc in sorted(targets, key=lambda x: int(x[3:]))}


def sankoff_gain_count(root, carriers):
    """Minimum 0->1 gains, with the parent of the target node fixed to state 0."""
    cost = {}
    for node in root.find_clades(order='postorder'):
        if node.is_terminal():
            observed = 1 if node.name in carriers else 0
            cost[node] = (0 if observed == 0 else float('inf'),
                          0 if observed == 1 else float('inf'))
        else:
            vals = []
            for state in (0, 1):
                total = 0
                for child in node.clades:
                    total += min(cost[child][0] + (state != 0),
                                 cost[child][1] + (state != 1))
                vals.append(total)
            cost[node] = tuple(vals)

    root_state = min((0, 1), key=lambda s: (cost[root][s] + (s != 0), s != 0))
    gains = int(root_state == 1)
    stack = [(root, root_state)]
    while stack:
        node, state = stack.pop()
        for child in node.clades:
            child_state = min(
                (0, 1),
                key=lambda s: (cost[child][s] + (state != s), state != s, s),
            )
            gains += int(state == 0 and child_state == 1)
            stack.append((child, child_state))
    return gains


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('indel_dir', type=Path)
    parser.add_argument('tree', type=Path)
    parser.add_argument('target_nodes', type=Path)
    parser.add_argument('genbank', type=Path)
    parser.add_argument('output_dir', type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    targets = read_targets(args.target_nodes)
    tree = Phylo.read(args.tree, 'newick')
    name_to_clade = {c.name: c for c in tree.find_clades() if c.name}

    # DCC labels are intentionally omitted: states are assigned solely by target-node descendants.
    event_samples, _, assignments = read_events(args.indel_dir, {})
    available = set(assignments['sample'])
    _, _, features, starts = load_features(args.genbank)

    roots, dcc_tips = {}, {}
    for dcc, node_name in targets.items():
        if node_name not in name_to_clade:
            raise ValueError(f'{node_name} for {dcc} not found in tree')
        root = name_to_clade[node_name]
        roots[dcc] = root
        dcc_tips[dcc] = {x.name for x in root.get_terminals() if x.name in available}

    gene_counts = defaultdict(lambda: defaultdict(int))
    gene_meta = {}
    for event, samples in event_samples.items():
        s, e, _, _, hits = annotate_event(event, features, starts)
        event_counts = {}
        for dcc, root in roots.items():
            carriers = samples & dcc_tips[dcc]
            event_counts[dcc] = sankoff_gain_count(root, carriers) if carriers else 0
        if not any(event_counts.values()):
            continue
        for hit in hits:
            tag = hit['locus_tag'] or f'intergenic:{s}-{e}'
            gene_meta[tag] = hit
            for dcc, count in event_counts.items():
                gene_counts[tag][dcc] += count

    dccs = list(targets)
    rows = []
    for tag, counts in gene_counts.items():
        row = {'locus_tag': tag, 'gene': gene_meta[tag].get('gene', '')}
        for dcc in dccs:
            row[f'{dcc}_event_count'] = counts[dcc]
        row['total_event_count'] = sum(counts.values())
        rows.append(row)
    result = pd.DataFrame(rows).sort_values('total_event_count', ascending=False)
    result.to_csv(args.output_dir/'gene_indel_event_counts_by_target_node_tree.csv', index=False)

    qc = pd.DataFrame([
        {'DCC': dcc, 'target_node': targets[dcc],
         'tree_tips': len(roots[dcc].get_terminals()),
         'tips_with_indel_file': len(dcc_tips[dcc])}
        for dcc in dccs
    ])
    qc.to_csv(args.output_dir/'target_node_indel_sample_counts.csv', index=False)
    print(qc.to_string(index=False))
    print('\nTop genes:')
    print(result.head(30).to_string(index=False))


if __name__ == '__main__':
    main()
