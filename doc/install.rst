Install
===================

cellsnp-lite is implemented in C. You can install it via conda_ or from this github repo.

Install via conda (latest stable version)
-----------------------------------------

This is the recommended way to install cellsnp-lite.

Step 1: add config
^^^^^^^^^^^^^^^^^^

.. code-block:: bash

  conda config --add channels bioconda
  conda config --add channels conda-forge

Step 2: install
^^^^^^^^^^^^^^^

to your current environment:

.. code-block:: bash

  conda install cellsnp-lite

or to a new environment:

.. code-block:: bash

  conda create -n CSP cellsnp-lite     # you can replace 'CSP' with another env name.

.. _conda: https://docs.conda.io/en/latest/

Install from this Github Repo (latest stable/dev version)
---------------------------------------------------------

We recommend installing cellsnp-lite via conda, as described above. The method of compiling
from source code (ie., installing from this repo) is described below:

For better demonstration, the whole installation process is assumed to be performed under 
a root dir ``~/tools/`` (you could change this dir based on your demands)

Step 0: install dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

cellsnp-lite depends on `zlib`_ and `htslib`_. The two libs should have been installed into
the system before installing cellsnp-lite. 

You could often find that zlib has already been installed in your system. If not, it's easy to 
install it with package management tools (eg., yum for CentOS). 

Following the `htslib_instruction`_, you could easily install htslib by

.. code-block:: bash

  cd ~/tools     # The root dir
  git clone https://github.com/samtools/htslib.git
  cd htslib
  autoheader     # If using configure, generate the header template...
  autoconf       # ...and configure script (or use autoreconf to do both)
  ./configure    # Optional but recommended, for choosing extra functionality
  make
  
  sudo make install   # Optional

If htslib is successfully installed, there should be some executable files (e.g., bgzip) 
and library files (eg., libhts.so) in htslib dir.

Step 1: compiling cellsnp-lite
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once zlib and htslib have been installed, it's ready to compile cellsnp-lite,

.. code-block:: bash

  cd ~/tools           # The root dir
  git clone https://github.com/single-cell-genetics/cellsnp-lite.git;
  cd cellsnp-lite;
  make;
  
  sudo make install;   # Optional

By default, this will build against an HTSlib source tree in ``../htslib``. You can alter this 
to a source tree elsewhere or to a previously-installed HTSlib by running 
``make htslib_dir=<path_to_htslib_dir>``.

Step 2: one more step config
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After compling cellsnp-lite, if you met the error ``error while loading shared libraries: libhts.so.3`` when running cellsnp-lite, you could fix this by setting environment variable ``LD_LIBRARY_PATH``
to proper value,

.. code-block:: bash

  # in this example, abspath_to_htslib_dir is ~/tools/htslib
  echo 'export LD_LIBRARY_PATH=<abspath_to_htslib_dir>:$LD_LIBRARY_PATH' >> ~/.bashrc;
  source ~/.bashrc;

.. _zlib: http://zlib.net/
.. _htslib: https://github.com/samtools/htslib
.. _htslib_instruction: https://github.com/samtools/htslib#building-htslib

