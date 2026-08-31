# Evolutionary genomics of dominant circulating clones in *Mycobacterium abscessus*

This repository contains the custom analysis and visualization scripts used to
investigate the genomic evolution of dominant circulating clones (DCCs) in the
*Mycobacterium abscessus* complex (MAB).

The workflow was designed to address three related questions:

1. What constitutes a conservative and stable MAB core genome?
2. How do mutation, recombination, and selection change before and after DCC
   emergence?
3. Does accessory-gene turnover provide convergent signals associated with DCC
   formation?

The repository focuses on the custom downstream analyses developed for this
study. Large sequencing datasets and intermediate files are not distributed
here. Three phylogenetic trees used by the analyses are included under `data/`.

## Analysis overview

### 1. Representative sampling and pangenome analysis

Isolates were sampled within DCCs using diversity-preserving phylogenetic
subsampling. Non-DCC isolates were stratified by subspecies and geographic
origin. Panaroo-derived gene clusters were subsequently evaluated to define a
conservative stable core-genome set and to summarize core, soft-core, shell,
and cloud components.

Scripts: [`scripts/01_sampling_pangenome`](scripts/01_sampling_pangenome)

### 2. Read-supported phylogenetic reconstruction

Assembly-derived variants were compared with read-mapping consensus calls and
classified as supported, inconsistent, or missing. These analyses quantify the
reliability of assembly-derived SNPs, identify problematic core-genome sites,
and support construction of the reference-coordinate alignment used for
Gubbins and DCC phylogenetic analyses.

Scripts: [`scripts/02_read_supported_phylogeny`](scripts/02_read_supported_phylogeny)

### 3. Core-gene evolutionary dynamics

Mutation-derived and recombination-derived SNPs were extracted from
node-labelled Gubbins results. Gene-level recombination-to-mutation ratios
(`r/m`) and mutation-only nonsynonymous-to-synonymous ratios (`pN/pS`) were
compared between DCC stem branches and post-emergence descendant branches.
The module also includes posterior selection-switch inference, convergence
summaries across DCCs, functional enrichment, branch-quality control, and
protein-position lollipop plots.

Scripts: [`scripts/03_core_gene_evolution`](scripts/03_core_gene_evolution)

### 4. Accessory-genome evolution

Accessory-gene presence/absence patterns were used to reconstruct gene gains
and losses on DCC stem branches, identify DCC-associated losses, and summarize
the predicted taxonomic origins of accessory genes.

Scripts: [`scripts/04_accessory_genome`](scripts/04_accessory_genome)

## Repository structure

```text
MAB_core_methods/
|-- data/
|   |-- read_supported_phylogeny/
|   `-- sampling_pangenome/
|-- scripts/
|   |-- 01_sampling_pangenome/
|   |-- 02_read_supported_phylogeny/
|   |-- 03_core_gene_evolution/
|   `-- 04_accessory_genome/
|-- .gitattributes
`-- README.md
```

Module-specific README files describe the analysis order, input expectations,
and principal outputs in greater detail.

## Included phylogenetic data

The repository includes only lightweight tree files required by the archived
analyses:

- `data/read_supported_phylogeny/MAB.node_labelled.final_tree.tre`
- `data/sampling_pangenome/core.tree`
- `data/sampling_pangenome/accessory_matrix.fasta.tree`

Raw reads are available from the public accessions reported in the associated
study. Large alignments, BAM files, Panaroo outputs, Gubbins intermediate files,
and derived spreadsheets are excluded because of their size.

## Software requirements

The custom scripts require Python 3 and R. The main Python dependencies are:

- `numpy`
- `pandas`
- `scipy`
- `biopython`
- `ete3`
- `openpyxl`
- `matplotlib`
- `seaborn`

The main R dependencies are:

- `tidyverse`
- `ape`
- `ggplot2`
- `patchwork`
- `scales`
- `openxlsx`
- `trackViewer`
- `GenomicRanges`

Upstream analyses additionally used external tools including SPAdes, Prokka,
Panaroo, read-mapping and variant-calling software, Gubbins, IQ-TREE, TreeMmer,
eggNOG-mapper, and standard phylogenetic utilities. Exact versions and command
parameters should be taken from the Materials and Methods of the associated
manuscript.

## Usage

All archived scripts are written in English and avoid hard-coded local paths.
Most scripts accept input and output paths as command-line arguments. Display a
script's usage information by running it without arguments, for example:

```bash
python3 scripts/03_core_gene_evolution/dcc_stem_post_rm_by_gene.py
```

or:

```bash
Rscript scripts/03_core_gene_evolution/13_plot_pns_switch_scatter.R
```

Because the analyses operate on different upstream result formats, this is a
collection of modular workflows rather than a single end-to-end pipeline.
Consult the README in each module before running its scripts.

## Key analytical conventions

- DCC stem events represent changes assigned to the branch leading to the DCC
  ancestral node; post-emergence events include descendant internal and
  terminal branches.
- Recombination-derived SNPs contribute to `r/m` but are excluded from the
  mutation-only `pN/pS` analysis.
- Repeated changes on separate branches are retained as independent
  evolutionary events.
- Selection switches are evaluated using posterior probabilities derived from
  observed synonymous and nonsynonymous events and adjusted mutational
  opportunities.
- The stable core-genome definition is deliberately conservative and should be
  interpreted as a lower-bound core set rather than an exhaustive list of all
  biologically conserved genes.

## Reproducibility scope

This repository preserves the custom scripts used for the main downstream
analyses and figures. It does not currently include the original batch commands
or environment files for every upstream software stage. In particular, the
original TreeMmer execution command, independent stable-core validation script,
and consensus-to-Gubbins alignment construction script should be recovered from
the analysis environment for a fully automated raw-read-to-figure workflow.

## Citation

If you use these scripts, please cite the associated manuscript. Full citation
information will be added after publication.

## Contact

Questions and reproducibility issues can be submitted through the GitHub issue
tracker.
