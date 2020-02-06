========================
Computational efficiency
========================

In theory, the computational complexity (i.e., running time) of cellSNP is O(n) 
for number of variants and O(n*log(n)) for number of cells, which means it is 
more sensitive to the cell counts.

Roughly, for a common 10x sample with 15K cells, cellSNP genotypes ~7 million 
variants with 15 CPUs in around 20 hours. In case you have more cells or more 
variants to genotype, you could consider split the variants into multiple sets 
and run it on a cluster server.

For `human SNP list`_, we suggest using the version with AF5e2 (i.e., AF>5%, 7.4M 
SNPS), instead of AF5e4 (i.e., AF>0.05%, 36.6M SNPs).

.. _human SNP list: https://sourceforge.net/projects/cellsnp/files/SNPlist/

