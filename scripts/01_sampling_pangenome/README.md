# DCC phylogeny, TreeMmer sampling, and pangenome reconstruction

## Workflow order

1. Reconstruct an unsampled phylogeny for every DCC from the complete isolate set.
2. Apply TreeMmer independently to each DCC tree using the prespecified retained-tree-length target.
3. Sample Non-DCC isolates by subspecies and country or region.
4. Annotate the 1,130 selected assemblies with Prokka and run Panaroo in strict mode.
5. Define the initial core at >=99% prevalence in the representative set.
6. Evaluate every initial core-gene alignment and remove poorly aligned loci.
7. Validate retained genes in the independent, nonsampled isolate set and retain genes present in >99% of genomes.
8. Use the resulting 3,001 genes as a conservative lower-bound stable core-genome set.

## Numbered scripts

- `01_plot_treemmer_sampling_schematic.R`: figure-only schematic of diversity-preserving sampling; it does not perform TreeMmer pruning.
- `02_evaluate_core_gene_alignments.py`: evaluates Panaroo per-gene alignments and classifies alignment quality.
- `03_summarize_DCC_gene_categories.py`: calculates core, soft-core, shell, and cloud gene counts within DCCs.
- `04_plot_pangenome_filtering_workflow.R`: plots the core-genome filtering workflow and retained gene counts.
- `05_plot_DCC_pangenome_composition_polar.R`: plots DCC-level pangenome composition.
- `06_plot_DCC_pangenome_composition_absolute.py`: plots absolute pangenome category counts.
- `07_summarize_DCC_pangenome_categories.py`: summarizes DCC-level pangenome categories from Panaroo outputs.

## Core filtering rule implemented by script 02

For each alignment column, a gap fraction >0.5 defines a gap-dominated column. A gene is classified as high quality when gap-dominated columns account for <=0.1 of its alignment length. This is more precise than describing the rule simply as "gene gap proportion <10%".

## Inputs

- Per-DCC unsampled phylogenetic trees.
- TreeMmer-selected isolate lists.
- Prokka GFF files for 1,130 representative genomes.
- Panaroo `gene_presence_absence.Rtab` and per-gene alignments.
- Presence/absence results for the independent validation set.

## Principal outputs

- 1,130 representative genomes.
- 34,571 Panaroo gene clusters.
- 3,436 initial core genes.
- 3,195 alignment-QC-passed genes.
- 3,001 independently validated stable core genes.

## Missing reproducibility components

The original TreeMmer execution command, isolate-selection table, Prokka batch command, Panaroo command, and independent 10,184-genome validation script were not found in the local script inventory. Their result files are preserved, but these commands/scripts should be recovered from the server before public GitHub release.
