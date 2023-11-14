
This directory stores scripts and data for processing the input and output
of [cellsnp-lite](https://github.com/single-cell-genetics/cellsnp-lite).

The "csp_utils.R" is an R script. It contains several functions:
- `update_cellsnp_matrices()`: update three sparse matrices `AD`, `DP`, and 
  `OTH` given a VCF file storing a list of SNPs that passed filtering.

