"""
cellSNP - Analysis of expressed alleles in single cells
See: https://github.com/huangyh09/cellSNP
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

from setuptools.extension import Extension
from Cython.Build import cythonize
from sys import version_info as sys_version

here = path.abspath(path.dirname(__file__))

# Set __version__ for the project.
exec(open("./cellSNP/version.py").read())

# Get the long description from the relevant file
with open(path.join(here, 'README.rst'), encoding='utf-8') as f:
	long_description = f.read()
	
reqs = ['numpy>=1.9.0', 'pysam>=0.15.2', 'cython>=0.29.16'] #, 'h5py'

# refer to the setup.py in pysam.
extra_compile_args = [
    "-Wno-unused",
    "-Wno-sign-compare",
]

# List cython extensions in order.
ext_modules = [
	dict(name = "cellSNP.utils.base_utils",
		language = "c",
		sources = [path.join('cellSNP', 'utils', 'base_utils.pyx')],
		libraries = []),
	dict(name = "cellSNP.utils.pileup_utils",
		language = "c",
		sources = [path.join('cellSNP', 'utils', 'pileup_utils.pyx')],
		libraries = []),
	dict(name = "cellSNP.utils.pileup_regions",
		language = "c",
		sources = [path.join('cellSNP', 'utils', 'pileup_regions.pyx')],
		libraries = []),
	dict(name = "cellSNP.utils.vcf_utils",
		language = "c",
		sources = [path.join('cellSNP', 'utils', 'vcf_utils.pyx')],
		libraries = [])
]

ext_module_common_options = dict(
	extra_compile_args = extra_compile_args,
	include_dirs = [path.abspath(x) for x in ["cellSNP", "."]],
	library_dirs = []
)

for module in ext_modules:
	module.update(**ext_module_common_options)

setup(
	name='cellSNP',

	# Versions should comply with PEP440.  For a discussion on single-sourcing
	# the version across setup.py and the project code, see
	# https://packaging.python.org/en/latest/single_source_version.html
	version=__version__,

	description='cellSNP - Analysis of expressed alleles in single cells',
	long_description=long_description,

	# The project's main homepage.
	url='https://github.com/huangyh09/cellSNP',

	# Author details
	author='Yuanhua Huang',
	author_email='yuanhua@ebi.ac.uk',

	# Choose your license
	license='Apache-2.0',

	# What does your project relate to?
	keywords=['allelic expression', 'single-cell RNA-seq'],

	# You can just specify the packages manually here if your project is
	# simple. Or you can use find_packages().
	packages=find_packages(),

	entry_points={
		  'console_scripts': [
			  'cellSNP = cellSNP.cellSNP:main'
			  ],
		  }, 

	# List run-time dependencies here.  These will be installed by pip when
	# your project is installed. For an analysis of "install_requires" vs pip's
	# requirements files see:
	# https://packaging.python.org/en/latest/requirements.html
	
	install_requires=reqs,

	py_modules = ['cellSNP'],

	# Cython extensions
	ext_modules = cythonize([Extension(**opts) for opts in ext_modules], language_level = sys_version.major),

	# Setting 'use_2to3 = True' to provide a facility to invoke 2to3 on the code as a part of the build process.
	use_2to3 = True

	# buid the distribution: python setup.py sdist
	# upload to pypi: twine upload dist/...

)