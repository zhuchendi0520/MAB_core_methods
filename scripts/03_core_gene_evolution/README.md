# Core-gene evolutionary analyses

This module contains the downstream analyses of mutation-derived and recombination-derived SNPs, DCC stem-versus-post-expansion r/m estimates, mutation-only pN/pS estimates, selection-switch classification, DCC-level convergence summaries, and gene-level lollipop plots.

## Numbered pN/pS workflow

1. `01_calculate_pns_from_raw_stem_post_events.py`
   - Input: `DCC*.raw_stem_post_events.csv` generated from the node-labelled Gubbins results.
   - Retains only `region == mutation`; SNPs assigned to recombinant regions are excluded.
   - Counts synonymous and nonsynonymous events independently on the DCC stem branch (`pre`) and all descendant branches (`post`).
   - Corrects mutational opportunities using the reference codon and the empirical nucleotide substitution matrix embedded in the script.
   - Output: per-DCC, merged, and DCC1-7 pooled gene-level pN/pS tables.

2. `02_add_pns_fdr_to_wide_xlsx.py`
   - Uses a two-sided exact binomial test for deviation of the observed nonsynonymous count from its codon- and mutation-spectrum-adjusted expectation.
   - Applies the Benjamini-Hochberg correction independently to each analysis block.
   - Output: the input workbook with `pNS_p_value` and `pNS_FDR` columns.

3. `03_calculate_pns_switch_posterior.py`
   - Models the nonsynonymous-event fraction in each phase with a Beta posterior using a uniform Beta(1,1) prior.
   - Converts posterior draws to codon- and mutation-spectrum-adjusted omega values.
   - Reports posterior medians, 95% credible intervals, directional probabilities, and the probability that omega crosses 1 between the pre- and post-expansion phases.
   - The principal high-confidence switch threshold used in Figure 3 is posterior switch probability >0.90.

4. `04_calculate_and_plot_DCC_TBL_TRL.R`
   - Reads the node-labelled Gubbins final tree and identifies all descendant tips of each predefined DCC ancestral node.
   - TBL is the length of the single terminal edge leading to each isolate.
   - TRL is the sum of all edge lengths from the DCC ancestral node to each descendant tip.
   - The final analysis uses `Node_640` for DCC1 and `Node_844` for DCC3; the remaining DCC roots are read from `target_nodes.txt`.
   - Outputs isolate-level measurements, DCC-level summary statistics, and TBL/TRL boxplots.

5. `05_calculate_dcc_pns_switch_matrix_fullpost.py`
   - Maps full-cohort post-expansion pN/pS results to ATCC 19977 genes and classifies DCC-specific selection switches.

6. `06_rm_direction_eggnog_enrichment.py`
   - Joins r/m direction classes to eggNOG annotations and performs pathway/module enrichment analyses.

7. `07_prepare_lollipop_inputs.py` and `08_plot_trackviewer_lollipops.R`
   - Prepare mutation-event tables and draw protein-position lollipop plots for selected genes.

8. `09_prepare_kegg_enrichment_plot_data.py` and `10_plot_kegg_enrichment_dotplot.R`
   - Prepare and plot KEGG enrichment results for genes with increased or decreased r/m.

9. `11_plot_global_vs_dcc_internal_rm_by_gene_class.R`
   - Compares global and DCC-internal r/m values within the four HMM gene classes.

10. `12_plot_rm_pre_post_by_gene_class.R`
   - Draws paired pre/post-expansion r/m distributions stratified by gene class and direction.

11. `13_plot_pns_switch_scatter.R`
   - Plots pre- versus post-expansion pN/pS and highlights high-confidence directional switches.

12. `14_clean_branch_outliers_and_recalculate.py`
   - Removes predefined anomalous post-expansion branches and rebuilds r/m and pN/pS summaries.

## Additional analyses

- `dcc_stem_post_rm_by_gene.py`: extracts stem and descendant branch events and summarizes gene-level r/m.
- `calculate_all_core_dcc_pns_direction_consistency.py`: quantifies direction consistency across DCCs.
- `plot_dcc_pns_switch_convergence_heatmap_fullpost_horizontal.R`: plots the DCC convergence heatmap.
- `plot_coupled_rm_pns_alluvial.py` and `plot_coupled_rm_pns_circular.py`: summarize coupled r/m and pN/pS states.
- `pns_long_to_pre_post_wide.py`: reshapes long-format pN/pS output.

## Important conventions

- Recombination-derived SNPs are analyzed through r/m and are not included in pN/pS.
- Repeated occurrences on different branches are retained as independent evolutionary events.
- When no synonymous event is observed, the calculation script uses one synonymous event only as a denominator stabilizer and reports this in `syn_set_to_1`; raw counts are retained separately.
- Genes without usable mutation events remain missing rather than being interpreted as evidence of neutrality or purifying selection.
