=======
History
=======

Release v0.0.5 (24/09/2018)
===========================
* pileup a list of positions with `pysam-fetch`, which may returns more
  reads than `pysam-pileup`. This feature requires further check
* change vcf file header to be more compatible with bcftools
* support turning cell-barcode off to return a sample level only

Release v0.0.4 (25/08/2018)
===========================
* pileup the whole genome for 10x single-cell RNA-seq data
* Note, post-filetering is needed as the current filtering doesn't 
  consider the heterozygous genotype for all donors.

