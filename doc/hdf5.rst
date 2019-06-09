Notes for using HDF5 in Python
==============================

VCF is the standard file format for genetic data, but single cells data is 
usually highly sparse, and loading sparse VCF is not supported by common 
packages, e.g. vcfR, hence saving genetic data into HDF5 format can provide a 
way to load sparse matrix directly.

String
------
An issue for many years: https://github.com/h5py/h5py/issues/289
See more here: http://docs.h5py.org/en/latest/strings.html

Object size
-----------
https://docs.python.org/3/library/sys.html#sys.getsizeof
https://pypi.org/project/objsize/


Some links
----------
https://www.pythonforthelab.com/blog/how-to-use-hdf5-files-in-python/

