=======
cellSNP
=======

cellSNP aims to pileup the expressed alleles in single-cell or bulk RNA-seq 
data, which can be directly used for donor deconvolution in multiplexed single-
cell RNA-seq data, particularly with cardelino_, an R package, which assigns 
cells to donors and detects doublets, even without genotyping the given donors.

cellSNP heavily depends on pysam_, a Python interface for samtools and bcftools. 
This program should give very similar results as samtools/bcftools mpileup, if 
it isn't the same. Also, there are two major differences comparing to bcftools 
mpileup:

1. cellSNP can pileup either the whole genome or a list of positions, with 
   directly splitting into a list of cell barcodes, e.g., for 10x genome. With 
   bcftools, you may need to manipulate the RG tag in the bam file first.
2. cellSNP uses simple filtering for outputting SNPs, i.e., total UMIs or counts
   and minor alleles fractions. The idea here is to keep most information of 
   SNPs and the downstream statistical model can handle adaptively.


Installation
------------

cellSNP is available through `pypi`_. To install, type the following command 
line, and add ``-U`` for upgrading:

.. code-block:: bash

  pip install cellSNP

Alternatively, you can download or clone this repository and type 
``python setup.py install`` to install. In either case, add ``--user`` if you 
don't have the permission as a root or for your Python environment.


Quick usage
-----------

Once installed, check all arguments by type ``cellSNP -h``. There are three 
modes of cellSNP:

**Mode 1**: pileup a list of common SNPs for single cells in a big BAM/SAM 
file. Require: a single BAM/SAM file, e.g., from cellranger, a VCF file for 
a list of common SNPs. This mode is recommended comparing to mode 2, if a 
list of common SNP is known, e.g., human.

.. code-block:: bash

  cellSNP -s $BAM -b $BARCODE -o $OUT_FILE -R $REGION_VCF -p 20
  
Recommend filtering SNPs with <20UMIs or <10% minor alleles for downstream 
donor deconvolution, by adding ``--minMAF 0.1 --minCOUNT 20``


**Mode 2**: pileup the whole genome for single cells in a big BAM/SAM file. 
This mode may give uninformative SNPs, but can be useful when the data set 
is highly sparse.

.. code-block:: bash

  cellSNP -s $BAM -b $BARCODE -o $OUT_FILE -p 22
  
Recommend filtering SNPs with <100UMIs or <10% minor alleles for saving space
and speed up inference when pileup whole genome: ``--minMAF 0.1 --minCOUNT 100``


**Mode 3**: pileup a list of common SNPs for one or multiple bulk BAM/SAM files.
Require: one or multiple BAM/SAM files, their according sample ids, and a VCF 
file for a list of common SNPs.

.. code-block:: bash

  cellSNP -s $BAM1,$BAM2,$BAM3 -I sample_id1,sample_id2,sample_id3 -o $OUT_FILE -R $REGION_VCF -p 20
  
Set filtering thresholds according to the downstream analysis.



.. note::

   - For lift over SNP positions in vcf file from one genome build to another, 
     see our `LiftOver_vcf`_ wrap function.
   - For release notes, see `release.rst`_.

.. _cardelino: https://github.com/PMBio/cardelino
.. _pysam: https://github.com/pysam-developers/pysam
.. _pypi: https://pypi.org/project/cellSNP/
.. _LiftOver_vcf: https://github.com/huangyh09/cellSNP/tree/master/liftOver
.. _release.rst: https://github.com/huangyh09/cellSNP/blob/master/doc/release.rst

