Install
===================
Cellsnp-lite is implemented in C. You can install it via conda_ or from this github repo.

Install via conda (latest stable version)
-----------------------------------------
This is the recommended way to install cellsnp-lite. Lacking the potential issues of 
dependency, it's simple and fast if conda is available on the machine.

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

Step 0: install dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Cellsnp-lite mainly depends on `zlib`_ and `HTSlib`_ (v1.10.2+). Note that HTSlib 
has some extra dependencies: ``liblzma``, ``libbz2``, ``libcurl``, and 
``libcrypto``. The whole list of dependencies of building cellsnp-lite is:

- gcc (have tested v4.8.5 on CentOS 7)
- autoreconf
- zlib
- HTSlib >= 1.10.2
- liblzma
- libbz2
- libcurl
- libcrypto

All dependencies should have been installed into the system before installing cellsnp-lite.

If you already have a pre-installed HTSlib then go to next step. Otherwise, a common way 
to install HTSlib is building from its Github repo following this `HTSlib_instruction`_.
When HTSlib is successfully installed, there should be some executable files (e.g., bgzip) 
and library files (eg., libhts.a and libhts.so) in htslib dir.

Step 1: compiling cellsnp-lite
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once all dependencies have been installed, it's ready to compile cellsnp-lite,

.. code-block:: bash

  git clone https://github.com/single-cell-genetics/cellsnp-lite.git
  cd cellsnp-lite
  autoreconf -iv
  ./configure     # Needed for choosing optional functionality
  make
  
  make install    # Optional, may need sudo privilege

By default, this will either build against a pre-installed HTSlib in ``../htslib`` 
or build with a HTSlib in a system library path (eg., ``/usr/lib``). You can alter 
this to a pre-installed HTSlib elsewhere by configuring with ``--with-htslib=DIR``.
The ``DIR`` should be either a dir containing a source tree with libhts.a/libhts.so
being in ``DIR``, or a dir containing ``include`` and ``lib`` subdirs with 
libhts.a/libhts.so being in ``DIR/lib``. Note that ``DIR`` must be absolute path and
please use ``/home/<user>`` instead of ``~`` if needed.

Possible issue
**************
Compilation in Step 1 prefers ``libhts.a`` than ``libhts.so`` for linking HTSlib. In rare
cases that the libhts.a does not exist and libhts.so has to be used for linking, the 
issue ``error while loading shared libraries: libhts.so.3`` could happen when running 
cellsnp-lite, if HTSlib is not in system library path (eg., /usr/lib).

The issue, if happed, could be fixed by adding abspath to the dir containing libhts.so
to the environment variable ``LD_LIBRARY_PATH``,

.. code-block:: bash

  echo 'export LD_LIBRARY_PATH=<abspath_to_htslib_dir>:$LD_LIBRARY_PATH' >> ~/.bashrc
  source ~/.bashrc

.. _zlib: http://zlib.net/
.. _HTSlib: https://github.com/samtools/htslib
.. _HTSlib_instruction: https://github.com/samtools/htslib#building-htslib

