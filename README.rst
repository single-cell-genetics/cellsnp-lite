=======
cellSNP
=======

cellSNP aims to pileup the expressed alleles in single-cell or bulk RNA-seq 
data, which can be directly used for donor deconvolution in multiplexed single-
cell RNA-seq data, particularly with vireo_ in cardelino_ R package, which 
assigns cells to donors and detects doublets, even without genotyping reference.

cellSNP heavily depends on pysam_, a Python interface for samtools and bcftools. 
This program should give very similar results as samtools/bcftools mpileup. 
Also, there are two major differences comparing to bcftools mpileup:

1. cellSNP can pileup either the whole genome or a list of positions, with 
   directly splitting into a list of cell barcodes, e.g., for 10x genome. With 
   bcftools, you may need to manipulate the RG tag in the bam file if you want 
   to divide reads into cell barcode groups.
2. cellSNP uses simple filtering for outputting SNPs, i.e., total UMIs or counts
   and minor alleles fractions. The idea here is to keep most information of 
   SNPs and the downstream statistical model can take the full use of it.


Installation
------------

cellSNP is available through `pypi`_. To install, type the following command 
line, and add ``-U`` for upgrading:

.. code-block:: bash

  pip install cellSNP

Alternatively, you can download or clone this repository and type 
``python setup.py install`` to install. In either case, add ``--user`` if you 
don't have the permission as a root or for your Python environment.

**Note**, cellSNP (>=0.1.0 requires pysam>=0.15.2), so make sure you are using 
the right version of `pysam`. Try `pip uninstall pysam` and then reinstall 
`pip install -U pysam`


Quick usage
-----------

Once installed, check all arguments by type ``cellSNP -h`` (see a snapshot_)
There are three modes of cellSNP:

* **Mode 1: pileup a list of SNPs for single cells in a big BAM/SAM file**

Require: a single BAM/SAM file, e.g., from cellranger, a VCF file for 
a list of common SNPs. This mode is recommended comparing to mode 2, if a 
list of common SNP is known, e.g., human (see Candidate SNPs below)

.. code-block:: bash

  cellSNP -s $BAM -b $BARCODE -o $OUT_FILE -R $REGION_VCF -p 20
  
Recommend filtering SNPs with <20UMIs or <10% minor alleles for downstream 
donor deconvolution, by adding ``--minMAF 0.1 --minCOUNT 20``


* **Mode 2: pileup the whole genome for single cells in a big BAM/SAM file**

.. code-block:: bash

  cellSNP -s $BAM -b $BARCODE -o $OUT_FILE -p 22
  
Recommend filtering SNPs with <100UMIs or <10% minor alleles for saving space
and speed up inference when pileup whole genome: ``--minMAF 0.1 --minCOUNT 100``

Note, this mode may output false positive SNPs, for example somatic variants or 
falses caussed by RNA editing. These false SNPs are probably not consistent in 
all cells within one individual, hence confounding the demultiplexing. 
Nevertheless, for species, e.g., zebrafish, without a good list of common SNPs, 
this strategy is still worth a good try, and it does not take much more time 
than mode 1.

* **Mode 3: pileup a list of SNPs for one or multiple bulk BAM/SAM files**

Require: one or multiple BAM/SAM files, their according sample ids, and a VCF 
file for a list of common SNPs.

.. code-block:: bash

  cellSNP -s $BAM1,$BAM2,$BAM3 -I sample_id1,sample_id2,sample_id3 -o $OUT_FILE -R $REGION_VCF -p 20
  
Set filtering thresholds according to the downstream analysis.


List of candidate SNPs
----------------------

A quality list of candidate SNPs (ususally common SNPs) are important for mode 1
and mode 3. If a list of genotyped SNPs is available, it can be used to pile up.
Alternatively, for human, common SNPs in population that have been idenetified 
from consortiums can also be very good candidates, e.g., gnomAD_ and 
1000_Genome_Project_. For the latter, we have compiled a list of 37 million 
common variants with this bash script_ and stored in this folder_.

In case you want to lift over SNP positions in vcf file from one genome build 
to another, see our `LiftOver_vcf`_ wrap function.


Release Notes
-------------

All releases are included in pypi_. Notes for each release are recorded in
`release.rst`_.

.. _vireo: https://rawgit.com/PMBio/cardelino/master/inst/doc/vignette-donorid.html
.. _cardelino: https://github.com/PMBio/cardelino
.. _snapshot: https://github.com/huangyh09/cellSNP/blob/master/doc/manual.rst
.. _pysam: https://github.com/pysam-developers/pysam
.. _pypi: https://pypi.org/project/cellSNP/
.. _gnomAD: http://gnomad.broadinstitute.org
.. _1000_Genome_Project: http://www.internationalgenome.org
.. _script: https://github.com/huangyh09/cellSNP/blob/master/SNPlist_1Kgenome.sh
.. _folder: https://sourceforge.net/projects/cellsnp/files/SNPlist/
.. _LiftOver_vcf: https://github.com/huangyh09/cellSNP/tree/master/liftOver
.. _release.rst: https://github.com/huangyh09/cellSNP/blob/master/doc/release.rst
