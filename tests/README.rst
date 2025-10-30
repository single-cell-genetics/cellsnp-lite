===============
Testing cellSNP
===============


10x Genomics data
=================

Run testing data
----------------
* Testing bash script: `test_10x.sh`_
* Run script for cellSNP in mode 1 & 2 with or without downloading data (125Mb).
  Note, for mode 2, it uses 22 CPUs by default, i.e., one CPU per chromosome in 
  parallel. If your machine doesn't have that many CPUs, make a smaller one, 
  e.g., `-p 4`.

  .. code-block:: bash

     bash test_10x.sh TRUE # if you need to download
     bash test_10x.sh # if you already downloaded
     
     
Generating test files
---------------------

* Script for generating testing files (for developer only): `data_maker_10x.sh`_
* Testing data generated and stored: http://ufpr.dl.sourceforge.net/project/cellsnp/testdata/cellSNP_testdata_10x.zip

Alternatives
------------

Here are example of alternative methods, while a systematic comparison to 
cellSNP is still to come.

* Mode 1 alternative with `VarTrix`_ to pileup each cell:

  .. code-block:: bash
     
     vatrix --umi --mapq 30 -b $BAM -c $BARCODE --scoring-method coverage --threads 1 \
         --ref-matrix ref.mtx --out-matrix alt.mtx -v $REGION --fasta $FASTA --no-duplicates
         
         
* Mode 2 alternative with `freebayes`_ (+ `VarTrix`_) to identify heterozygous SNPs:

  .. code-block:: bash
  
     freebayes -C 0 -F 0 --fasta-reference $faFile $BAM > freebayes.vcf
     # vcffilter -f "QUAL > 20" freebayes.vcf | bgzip -c > freebayes.sorted.vcf.gz

.. _test_10x.sh: https://github.com/single-cell-genetics/cellSNP/blob/master/test/test_10x.sh
.. _data_maker_10x.sh: https://github.com/single-cell-genetics/cellSNP/blob/master/test/data_maker_10x.sh
.. _VarTrix: https://github.com/10XGenomics/vartrix
.. _freebayes: https://github.com/ekg/freebayes


Smart-seq data
==============

* Testing data is coming soon..

