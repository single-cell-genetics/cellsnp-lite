List of candidate SNPs
======================

A quality list of candidate SNPs (ususally common SNPs) are important for mode 1. If a list of genotyped SNPs is available, it can be used to pile up.
Alternatively, for human, common SNPs in population that have been idenetified
from consortiums can also be very good candidates, e.g., gnomAD_ and
1000_Genome_Project_. For the latter, we have compiled a list of 7.4 million
common variants (AF>5%) with this bash script_ and stored in this folder_.

In case you want to lift over SNP positions in vcf file from one genome build
to another, see our `LiftOver_vcf`_ wrap function.

.. _gnomAD: http://gnomad.broadinstitute.org
.. _1000_Genome_Project: http://www.internationalgenome.org
.. _script: https://github.com/single-cell-genetics/cellsnp-lite/blob/master/scripts/SNPlist_1Kgenome.sh
.. _folder: https://sourceforge.net/projects/cellsnp/files/SNPlist/
.. _LiftOver_vcf: https://github.com/single-cell-genetics/cellsnp-lite/blob/master/scripts/liftOver/liftOver_vcf.py

