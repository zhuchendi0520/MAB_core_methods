# Read-supported phylogenetic reconstruction

This module compares assembly-derived variants with read-mapping consensus calls, identifies supported, inconsistent, and missing sites, and reconstructs the final reference-coordinate alignment used for Gubbins analysis. Sites passing the read-depth and allele-frequency thresholds are retained; unresolved positions are represented as `N` in the consensus alignment. The resulting alignment is used for recombination inference and DCC phylogenetic reconstruction.

The exact mapper, variant caller, reference genome, Gubbins version, and tree-building parameters must be recorded when this Methods subsection is finalized.

