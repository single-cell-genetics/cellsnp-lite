Before releasing a new version (e.g., v100.0.1),
- update CSP_VERSION in src/config.h
- update doc/release.rst
- update doc/manual.rst if needed
- update doc/conf.py for readthedoc
- update AC_INIT in configure.ac

Files generated by autotools (`autoreconf -iv`) and related to bioconda-build
- autom4te.cache (this dir has been added into .gitignore and can be removed) 
- aclocal.m4
- conf.h.in  
- configure
- depcomp
- install-sh
- Makefile.in
- missing

