=======
History
=======

Release v0.0.8 (21/12/2018)
===========================
* Support output file in the same path of command line
* Support input file in format of .cram
* Update readme file, especially for compiled common variants from 1000 genome 
  project (https://sourceforge.net/projects/cellsnp/files/SNPlist/)

Release v0.0.7 (04/10/2018)
===========================
* change the header of the VCF file to be more suitable for bcftools
* realise the issue of heavy memory consuming, which even kills the 
  jobs in cluster. The menory taken increase linearly to the number 
  of processors used. When using 20 CUPs, >20G memory is recomended 
  for >5K cells. Solution for higher memory efficiency will be 
  proposed in future.

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

