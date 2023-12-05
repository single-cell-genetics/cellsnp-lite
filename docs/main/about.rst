..
   About
   =====


Cellsnp-lite: Efficient Genotyping Bi-Allelic SNPs on Single Cells
------------------------------------------------------------------
.. include:: /main/breief_introduction.rst


Implementation
--------------

About Genotyping Model
~~~~~~~~~~~~~~~~~~~~~~
For genotyping in single cells, cellsnp-lite first needs to know the ``REF`` 
and ``ALT`` alleles. 
These two alleles can be either specified by users (``-R`` option in mode 1), 
or de novo inferred from data (in mode 2). 
After that, cellsnp-lite will perform genotyping, to select the genotype with 
the maximum likelihood in each single cell, with the error model as presented 
in ``Table 1`` in `Jun et al, 2012`_.

The two options ``--minMAF`` and ``--minCOUNT`` are used for filtering SNPs 
in a pseudo-bulk manner, not in each single cell. 
The corresponding "MAF and COUNT" values are calculated based on aggregated 
read/UMI counts of all cells.

Depth of ``REF`` and ``ALT`` alleles indeed can be used to infer genotype. 
However, the accuracy of inference could be low in this way, 
due to the high noise in the sequencing data (e.g., sequencing error). 
For instance, in your first example (``"1/0:1:5:1:31,29,147:4,0,1,1,0"``), 
if the one read supporting ``ALT`` allele is an artifact arising from 
sequencing error, then the truth genotype could be ``0/0``.

To account for sequencing error in genotyping, the error model presented in 
``Table 1`` of `Jun et al, 2012`_ uses a parameter ``e`` indicating occurrence
of "Base Calling Error Event". 
Likihood can be simply treated as possibility. 
For each SNP, the likelihoods of three genotypes (``"0/0"``, ``"1/0"``, 
``"1/1"``) are calculated by aggregating the information provided by all 
bases/alleles (from pileup all supporting reads/UMIs) and their corresponding 
sequencing qualities (reflecting probability of sequencing error), 
modified from ``Equation 1`` in `Jun et al, 2012`_. 
The final reported genotype is the one with maximum likelihood.

See also: issue #109.


Benchmark
---------
Cellsnp-lite shows substantial improvement in running speed and memory 
usage with highly concordant results compared to existing 
methods.
For details, please checkout our paper:
https://doi.org/10.1093/bioinformatics/btab358


Limitations
-----------


.. _Jun et al, 2012: https://doi.org/10.1016%2Fj.ajhg.2012.09.004

