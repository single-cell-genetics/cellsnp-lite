..
   Data
   ====


Data
====

.. _data List of common SNPs:

List of candidate SNPs
----------------------
A quality list of candidate SNPs (ususally common SNPs) are important for 
mode 1. 
If a list of genotyped SNPs is available, it can be used to pile up.
Alternatively, for human, common SNPs in population that have been idenetified
from consortiums can also be very good candidates, e.g., gnomAD_ and
`1000 Genome Project`_. 
For the latter, we have compiled a list of 7.4 million common variants (AF>5%)
with script SNPlist_1Kgenome.sh_ and stored in `folder SNPlist`_.

In case you want to lift over SNP positions in vcf file from one genome build
to another, see our `LiftOver_vcf`_ wrap function.


.. _1000 Genome Project: http://www.internationalgenome.org
.. _folder SNPlist: https://sourceforge.net/projects/cellsnp/files/SNPlist/
.. _gnomAD: http://gnomad.broadinstitute.org
.. _LiftOver_vcf: https://github.com/single-cell-genetics/cellsnp-lite/blob/master/scripts/liftOver/liftOver_vcf.py
.. _SNPlist_1Kgenome.sh: https://github.com/single-cell-genetics/cellsnp-lite/blob/master/scripts/SNPlist_1Kgenome.sh

