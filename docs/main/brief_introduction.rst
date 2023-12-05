..
   Brief Introduction
   ==================


Cellsnp-lite is a C/C++ tool for efficient genotyping bi-allelic SNPs on
single cells.
You can use *cellsnp-lite* after read alignment to obtain the
*snp x cell* pileup UMI or read count matrices for each alleles of given or
detected SNPs.

The output from *cellsnp-lite* can be directly used for downstream analysis
such as:

#. Donor deconvolution in multiplexed single-cell RNA-seq data
   (e.g., with vireo_).
#. Allele-specific CNV analysis in single-cell or spatial transcriptomics data
   (e.g., with Numbat_ or XClone_).
#. Clonal substructure discovery using single cell mitochondrial variants
   (e.g., with MQuad_).

Cellsnp-lite has following features:

* **Wide applicability:** *cellsnp-lite* can take data from various omics as
  input, including RNA-seq, DNA-seq, ATAC-seq, either in bulk or single cells.
* **Simplified user interface** that supports parallel computing, cell barcode
  and UMI tags.
* **High efficiency** in terms of running speed and memory usage with highly
  concordant results compared to existing methods.

For details of the tool, please checkout our paper:

    Xianjie Huang, Yuanhua Huang, Cellsnp-lite: an efficient tool for
    genotyping single cells,
    Bioinformatics, Volume 37, Issue 23, December 2021, Pages 4569–4571,
    https://doi.org/10.1093/bioinformatics/btab358


.. _MQuad: https://github.com/single-cell-genetics/MQuad
.. _Numbat: https://github.com/kharchenkolab/numbat
.. _vireo: https://github.com/huangyh09/vireo
.. _XClone: https://github.com/single-cell-genetics/XClone

