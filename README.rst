=======
cellSNP
=======

Analysis of expressed alleles in single cells with pysam


Quick Start
-----------

**Installation**: 

- download or clone this repository
- type ``python setup.py install`` or  ``python setup.py develop`` if 
  upgrades frequently, e.g., for active development period
- add ``--user`` if you don't have root permission

**Command line**

- ``cellSNP -s $SAM -b $BARCODE -o $OUT_FILE -p 22``
- Type command line ``cellSNP -h``
- Recommand filtering SNPs with <100UMIs or <10% minor alleles for saving space
  and speed up inference when pileup whole genome: ``--minMAF 0.1 --minCOUNT 100``

**Notes**

- Fetch a list of SNPs is still under test
- `LiftOver_vcf`_ for SNP positions in vcf file between genome builds
- For lease notes, see `release.rst`_

.. _LiftOver_vcf: https://github.com/huangyh09/cellSNP/tree/master/liftOver
.. _release.rst: https://github.com/huangyh09/cellSNP/blob/master/doc/release.rst

