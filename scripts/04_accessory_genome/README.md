# Accessory-genome analyses

1. `01_reconstruct_gene_gain_loss_on_tree.py`
   - Reconstructs binary gene presence/absence states on a rooted phylogeny and reports gains and losses on specified DCC stem branches.

2. `find_DCC_root_gene_losses.py`
   - Identifies strict and relaxed DCC-root gene-loss candidates from the gene presence/absence matrix and phylogeny.

3. `summarize_plot_DCC_accessory_seed_taxonomy.R`
   - Summarizes eggNOG seed assignments by broad genus and by mycobacterial species and draws DCC-level heatmaps.

Core- and accessory-genome tree inference and functional annotation are external
Panaroo, phylogenetic, and eggNOG-mapper stages. Their resulting tree files are
stored under `data/sampling_pangenome`.
