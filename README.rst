============
cellsnp-lite
============

cellsnp-lite aims to pileup the expressed alleles in single-cell or bulk RNA-seq 
data, which can be directly used for donor deconvolution in multiplexed 
single-cell RNA-seq data, particularly with vireo_, which assigns cells to 
donors and detects doublets, even without genotyping reference.

cellsnp-lite heavily depends on htslib_. 
This program should give very similar results as samtools/bcftools mpileup. 
Also, there are two major differences comparing to bcftools mpileup:

1. cellsnp-lite can now pileup a list of positions, with 
   directly splitting into a list of cell barcodes, e.g., for 10x genome. With 
   bcftools, you may need to manipulate the RG tag in the bam file if you want 
   to divide reads into cell barcode groups.
2. cellsnp-lite uses simple filtering for outputting SNPs, i.e., total UMIs or counts
   and minor alleles fractions. The idea here is to keep most information of 
   SNPs and the downstream statistical model can take the full use of it.

cellsnp-lite is the C version of cellSNP_, which is implemented in Python. Compared to cellSNP, cellsnp-lite is basically more efficient with 
higher speed and less memory usage. But it can only pileup SNPs with a list of positions and cannot pileup 
whole chromosome(s), which is mode 2 of cellSNP. If you would like to pileup whole chromosome(s) for a single 
BAM/SAM file, please refer to cellSNP.

News
----
We have turn off the PCR duplicate filtering by default (--maxFLAG), as it is not well flagged in CellRanger, hence may result in loss of a substantial fraction of SNPs. Please use v0.3.1 (or above) or setting --maxFLAG to large number. Credits to issue13_ of cellSNP.

All release notes can be found in `doc/release.rst`_.

For computational efficiency, we initialised comments on this: `doc/speed.rst`_

.. _issue13: https://github.com/single-cell-genetics/cellSNP/issues/13
.. _doc/release.rst: https://github.com/single-cell-genetics/cellsnp-lite/blob/master/doc/release.rst
.. _doc/speed.rst: https://github.com/single-cell-genetics/cellsnp-lite/blob/master/doc/speed.rst

Installation
------------

cellsnp-lite is implemented in C. You can install it via conda_ or from this github repo.

* **Method 1: Install via conda (latest stable version)**

Step 1: add config

.. code-block:: bash

  conda config --add channels bioconda
  conda config --add channels conda-forge
  
Step 2: install  

to your current environment:

.. code-block:: bash

  conda install cellsnp-lite
  
or to a new environment:

.. code-block:: bash

  conda create -n CSP cellsnp-lite     # you can replace 'CSP' with another env name.

.. _conda: https://docs.conda.io/en/latest/

* **Method 2: Install from this Github Repo (latest stable/dev version)**

cellsnp-lite depends on `zlib`_ and `htslib`_. The two libs should have been installed in the system before
installing cellsnp-lite. Then to install cellsnp-lite,  

.. code-block:: bash

  git clone https://github.com/single-cell-genetics/cellsnp-lite.git;
  cd cellsnp-lite; 
  make;
  sudo make install;
  
By default, this will build against an HTSlib source tree in ../htslib. You can alter this to a source tree elsewhere or to a 
previously-installed HTSlib by running ``make htslib_dir=<path_to_htslib_dir>``.  

Besides, if you met the error ``error while loading shared libraries: libhts.so.3`` when running cellsnp-lite, you could fix this 
by setting environment variable ``LD_LIBRARY_PATH`` to proper value,

.. code-block:: bash

  echo 'export LD_LIBRARY_PATH=<abspath_to_htslib_dir>:$LD_LIBRARY_PATH' >> ~/.bashrc;
  source ~/.bashrc;
  
Quick usage
-----------

Once installed, check all arguments by type ``cellsnp-lite -h`` (see a snapshot_)
There are three modes of cellsnp-lite:

* **Mode 1: pileup a list of SNPs for a single BAM/SAM file**

Use both `-R` and `-b`. 

Require: a single BAM/SAM file, e.g., from cellranger, a list of cell barcodes,
a VCF file for common SNPs. This mode is recommended comparing to mode 2, if a 
list of common SNP is known, e.g., human (see Candidate SNPs below)

.. code-block:: bash

  cellsnp-lite -s $BAM -b $BARCODE -O $OUT_DIR -R $REGION_VCF -p 20 --minMAF 0.1 --minCOUNT 20 --genotype --gzip
  
As shown in the above command line, we recommend filtering SNPs with <20UMIs  
or <10% minor alleles for downstream donor deconvolution, by adding 
``--minMAF 0.1 --minCOUNT 20``

Besides, special care needs to be taken when filtering PCR duplicates for scRNA-seq data by 
setting maxFLAG to a small value, for the upstream pipeline may mark each extra read sharing 
the same CB/UMI pair as PCR duplicate, which will result in most variant data being lost. 
Due to the reason above, cellsnp-lite by default uses a large maxFLAG value to include PCR 
duplicates for scRNA-seq data when UMItag is turned on.

* **Mode 2: pileup whole chromosome(s) for a single BAM/SAM file**

This mode requires inputting a single bam file with either cell barcoded 
(add `-b`) or a bulk sample. It is not available now and will be supported in future. 
If you would like to use this mode, please refer to cellSNP_.

* **Mode 3: pileup a list of SNPs for one or multiple BAM/SAM files**

Use `-R` but not `-b`.

Require: one or multiple BAM/SAM files (bulk or smart-seq), their according 
sample ids (optional), and a VCF file for a list of common SNPs. BAM/SAM files 
can be input in comma separated way (`-s`) or in a list file (`-S`). 

.. code-block:: bash

  cellsnp-lite -s $BAM1,$BAM2,$BAM3 -I sample_id1,sample_id2,sample_id3 -O $OUT_DIR -R $REGION_VCF -p 20 --UMItag None --genotype --gzip

  cellsnp-lite -S $BAM_list_file -I sample_list_file -O $OUT_DIR -R $REGION_VCF -p 20 --UMItag None --genotype --gzip

Set filtering thresholds according to the downstream analysis. Please add 
``--UMItag None`` if your bam file does not have UMIs, e.g., smart-seq and bulk 
RNA-seq.


List of candidate SNPs
----------------------

A quality list of candidate SNPs (ususally common SNPs) are important for mode 1
and mode 3. If a list of genotyped SNPs is available, it can be used to pile up.
Alternatively, for human, common SNPs in population that have been idenetified 
from consortiums can also be very good candidates, e.g., gnomAD_ and 
1000_Genome_Project_. For the latter, we have compiled a list of 7.4 million 
common variants (AF>5%) with this bash script_ and stored in this folder_.

In case you want to lift over SNP positions in vcf file from one genome build 
to another, see our `LiftOver_vcf`_ wrap function.

.. _script: https://github.com/single-cell-genetics/cellsnp-lite/blob/master/SNPlist_1Kgenome.sh
.. _folder: https://sourceforge.net/projects/cellsnp/files/SNPlist/

FAQ and releases
----------------
For troubleshooting, please have a look of `FAQ.rst`_, and we welcome reporting 
any issue_.

.. _cellSNP: https://github.com/single-cell-genetics/cellSNP
.. _vireo: https://github.com/huangyh09/vireo
.. _zlib: http://zlib.net/
.. _htslib: https://github.com/samtools/htslib
.. _snapshot: https://github.com/single-cell-genetics/cellsnp-lite/blob/master/doc/manual.rst
.. _gnomAD: http://gnomad.broadinstitute.org
.. _1000_Genome_Project: http://www.internationalgenome.org
.. _LiftOver_vcf: https://github.com/single-cell-genetics/cellsnp-lite/blob/master/liftOver/liftOver_vcf.py
.. _FAQ.rst: https://github.com/single-cell-genetics/cellsnp-lite/blob/master/doc/FAQ.rst
.. _issue: https://github.com/single-cell-genetics/cellsnp-lite/issues
