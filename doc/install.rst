Install
===================

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

.. _zlib: http://zlib.net/
.. _htslib: https://github.com/samtools/htslib

