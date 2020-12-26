============
cellsnp-lite
============

|conda| |platforms| |license|

.. |conda| image:: https://anaconda.org/bioconda/cellsnp-lite/badges/version.svg
    :target: https://bioconda.github.io/recipes/cellsnp-lite/README.html
.. |platforms| image:: https://anaconda.org/bioconda/cellsnp-lite/badges/platforms.svg
   :target: https://bioconda.github.io/recipes/cellsnp-lite/README.html
.. |license| image:: https://anaconda.org/bioconda/cellsnp-lite/badges/license.svg
   :target: https://bioconda.github.io/recipes/cellsnp-lite/README.html

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

cellsnp-lite is the C version of cellSNP_, which is implemented in Python. Compared to 
cellSNP, cellsnp-lite is basically more efficient with higher speed and less memory usage. 

News
----

All release notes can be found in `doc/release.rst`_.

For computational efficiency, we initialised comments on this: `doc/speed.rst`_

A pre-compiled candidate SNP list for human is at `Candidate_SNPs`_.

.. _doc/release.rst: https://github.com/single-cell-genetics/cellsnp-lite/blob/master/doc/release.rst
.. _doc/speed.rst: https://github.com/single-cell-genetics/cellsnp-lite/blob/master/doc/speed.rst
.. _Candidate_SNPs: https://cellsnp-lite.readthedocs.io/en/latest/snp_list.html

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

cellsnp-lite depends on `zlib`_ and `htslib`_. The two libs should have been installed in 
the system before installing cellsnp-lite. Then to install cellsnp-lite,  

.. code-block:: bash

  git clone https://github.com/single-cell-genetics/cellsnp-lite.git;
  cd cellsnp-lite; 
  make;
  sudo make install;
  
By default, this will build against an HTSlib source tree in ../htslib. You can alter this 
to a source tree elsewhere or to a previously-installed HTSlib by running 
``make htslib_dir=<path_to_htslib_dir>``.  

Besides, if you met the error ``error while loading shared libraries: libhts.so.3`` when 
running cellsnp-lite, you could fix this by setting environment variable ``LD_LIBRARY_PATH`` 
to proper value,

.. code-block:: bash

  echo 'export LD_LIBRARY_PATH=<abspath_to_htslib_dir>:$LD_LIBRARY_PATH' >> ~/.bashrc;
  source ~/.bashrc;
  
Manual
------

The full manual is at `https://cellsnp-lite.readthedocs.io`_.

Also, type ``cellsnp-lite -h`` for all arguments with the version you are using.

.. _`https://cellsnp-lite.readthedocs.io`: https://cellsnp-lite.readthedocs.io

FAQ and releases
----------------
For troubleshooting, please have a look of `FAQ.rst`_, and we welcome reporting 
any issue_.

.. _cellSNP: https://github.com/single-cell-genetics/cellSNP
.. _vireo: https://github.com/huangyh09/vireo
.. _htslib: https://github.com/samtools/htslib
.. _FAQ.rst: https://github.com/single-cell-genetics/cellsnp-lite/blob/master/doc/FAQ.rst
.. _issue: https://github.com/single-cell-genetics/cellsnp-lite/issues
