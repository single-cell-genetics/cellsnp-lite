============
Cellsnp-lite
============

|conda| |platforms| |license|

.. |conda| image:: https://anaconda.org/bioconda/cellsnp-lite/badges/version.svg
    :target: https://bioconda.github.io/recipes/cellsnp-lite/README.html
.. |platforms| image:: https://anaconda.org/bioconda/cellsnp-lite/badges/platforms.svg
   :target: https://bioconda.github.io/recipes/cellsnp-lite/README.html
.. |license| image:: https://anaconda.org/bioconda/cellsnp-lite/badges/license.svg
   :target: https://bioconda.github.io/recipes/cellsnp-lite/README.html



Cellsnp-lite: Efficient Genotyping Bi-Allelic SNPs on Single Cells
------------------------------------------------------------------
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
* **High efficiency** in terms of running speed and memory usage, with highly
  concordant results compared to existing methods.

For details of the tool, please checkout our paper:

    Xianjie Huang, Yuanhua Huang, Cellsnp-lite: an efficient tool for 
    genotyping single cells, 
    Bioinformatics, Volume 37, Issue 23, December 2021, Pages 4569–4571, 
    https://doi.org/10.1093/bioinformatics/btab358


Installation
------------
Cellsnp-lite depends on several external libraries such as htslib_. 
We highly recommend installing cellsnp-lite via conda_ to avoid potential 
issues regarding dependency.

.. code-block:: bash

  conda install -c bioconda cellsnp-lite
  

Alternatively, you may also compile from source code. For details, please 
check `install from this github repo`_ in the user guide.


Manual
------
The full manual is available in the user guide at 
https://cellsnp-lite.readthedocs.io


FAQ and feedback
----------------
For troubleshooting, please have a look of `FAQ.rst`_, and we welcome reporting 
any issue_ for bugs, questions and new feature requests.


Acknowledgement
---------------
Cellsnp-lite heavily depends on htslib_ for accessing high-throughput 
sequencing data. 
In addition, it uses the ``kvec.h`` file (from klib_) for dynamic array
usage and the ``thpool.{h,c}`` files (from C-Thread-Pool_) for
thread pool management.


.. _C-Thread-Pool: https://github.com/Pithikos/C-Thread-Pool
.. _conda: https://docs.conda.io/en/latest/
.. _FAQ.rst: https://github.com/single-cell-genetics/cellsnp-lite/blob/master/docs/main/FAQ.rst
.. _htslib: https://github.com/samtools/htslib
.. _install from this github repo: https://cellsnp-lite.readthedocs.io/en/latest/install.html#install-from-this-github-repo-latest-stable-dev-version
.. _issue: https://github.com/single-cell-genetics/cellsnp-lite/issues
.. _klib: https://github.com/attractivechaos/klib
.. _MQuad: https://github.com/single-cell-genetics/MQuad
.. _Numbat: https://github.com/kharchenkolab/numbat
.. _vireo: https://github.com/huangyh09/vireo
.. _XClone: https://github.com/single-cell-genetics/XClone

