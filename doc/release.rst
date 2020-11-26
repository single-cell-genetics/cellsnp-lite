=======
History
=======

Release v1.1.0 (26/11/2020)
===========================
* split cellsnp_utils.h & general_utils.h & cellsnp.c into separate modules
  and scripts (small .h and .c files)
* bgzip temporary files

Release v1.0.1 (17/11/2020)
===========================
* fix a bug of m2 & k2 in csp_infer_allele() which could lead to error AD
  and error MAF calculation.
* fix a bug of not allocating space for sample_id
* add badges

Release v1.0.0 (15/10/2020)
===========================
* add mode 2
* replace --maxFLAG with --inclFLAG and --exclFLAG
* always filter unmapped reads

Release v0.3.1 (22/07/2020)
===========================
* turn off the PCR duplicate filtering by default (--maxFLAG), as it is not 
  well flagged in CellRanger, hence may result in loss of a substantial 
  fraction of SNPs.

Release v0.3.0 (05/06/2020)
===========================
* fix a bug of using read.qqual and read.query_alignment_sequence in pileup_bases() and 
  fetch_bases(), which could cause error when CIGAR string includes the 'I' (Insertion) op
* fix a bug of repetitive UMIs existing in different cells when grouping and counting UMIs
* the pos in the output vcf of mode 2 is switched from 0-based to 1-based

Release v0.1.8 (06/02/2020)
===========================
* Mode 2 now supports pile up chromosomes for a single bulk sample
* Mode 3 now supports multiple bam files in a list file

Release v0.1.7 (04/10/2019)
===========================
* fix a bug when chromosome is not in the bam file
* support barcodes.tsv.gz
* liftOver supports bgzip compress
* add vcf format to cellSNP.base.vcf.gz

Release v0.1.6 (14/07/2019)
===========================
* support saving to sparse matrices:
  Please use ``-O`` for out directory instead of ``-o`` for VCF output only. 
  Also, you can use ``sparseVCF.py`` to convert existing VCF.gz into sparse 
  matrices
* turn off breaking from warnings
* change P_error with BQ ranges [0.25, 45]
* h5py is not a required dependent package anymore

Release v0.1.5 (02/07/2019)
===========================
* fix a bug in qual_vector when the Quality (Phred) scores is 0, i.e., ASCII 
  Code "!", and it will give a P_error as 1, hence fail with log transformation.
  Now, I limited the P_error to 0.9999.

Release v0.1.4 (24/06/2019)
===========================
* use bgzip by default if bgzip is executable, otherwise use gzip
* change GL to PL: Phred-scaled genotype likelihoods and move to after OTH tag
* filter_reads is not in use anymore and the filtering is combined into 
  fetch_bases or pileup_bases
* slightly optimise the memory by not keeping all reads but only positions

Release v0.1.3 (12/06/2019)
===========================
* Fix a minor bug for when loading unzipped vcf file.

Release v0.1.2 (10/06/2019)
===========================
* turn off the defaul HDF5 file output, but keep it optional.

Release v0.1.1 (09/06/2019)
===========================
* support output in hdf5 format for sparse matrix. To convert .vcf.gz to hdf5 
  file, you can use this script: 
  https://github.com/huangyh09/cellSNP/blob/master/test/VCF_convert.py

Release v0.1.0 (21/05/2019)
===========================
* support the estimate the genotype and genotype likelihood for each cell.
  The GT is for 0/0, 1/0, 1/1, while the genotype likelihood is for 0/0, 1/0,
  1/1, and 0/0+1/0, 1/1+1/0.
  The genotype estimate is based on the this paper (table 1; same as supp table
  S3 in Demuxlet paper): https://doi.org/10.1016/j.ajhg.2012.09.004
* cell tag changed from CR to CB and the lane info is kept
* pileup whole genome uses the same reads filtering as pile up positions
* add test files (note, the bam file is 19G)
* require pysam>=0.15.2 to get the qqual for each base call in the reads


Release v0.0.8 (21/12/2018)
===========================
* update the default setting that UMItag is not in use in bulk RNA-seq, as UMI 
  is cell specific in pseudo-bulk RNA-seq, hence better turn it UMI off by
  default 
* support output file in the same path of command line
* support cram input file, besides bam/sam 
* update readme file, especially for processed common variants from 1000 genome 
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

