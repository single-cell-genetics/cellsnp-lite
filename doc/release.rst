=======
History
=======

Release v0.0.6 (29/09/2018)
===========================
* fix the bug in pileup a list of positions with ``pysam-fetch``: 
  input wrong REF and ALT bases.
* support pileup a list of positions for multiple bulk samples
* check liftOver works fine: the last part of the SNPs have matched
  REF in fasta file.
* polish the printout log: label the three modes: 
  
  * Mode 1: Pileup a list of positions for single cells (most common)
  * Mode 2: Pileup whole genome for single cells
  * Mode 3: Pileup a list of positions for (multiple) bulk sample(s)

Release v0.0.5 (24/09/2018)
===========================
* pileup a list of positions with ``pysam-fetch``, which may returns more
  reads than ``pysam-pileup``. This feature requires further check
* change vcf file header to be more compatible with bcftools
* support turning cell-barcode off to return a sample level only

Release v0.0.4 (25/08/2018)
===========================
* pileup the whole genome for 10x single-cell RNA-seq data
* Note, post-filetering is needed as the current filtering doesn't 
  consider the heterozygous genotype for all donors.

