
## Donot use @PACKAGE_NAME@ here otherwise error
bin_PROGRAMS = cellsnp-lite

## Donot use a variable, like @some_dir@, to stand for `src` dir;
## the @some_dir@ would not be expanded in xxx_SOURCES
cellsnp_lite_SOURCES = src/cellsnp.c src/csp.c \
    src/csp_fetch.c src/csp_pileup.c src/jfile.c \
    src/jsam.c src/jstring.c src/mplp.c \
    src/snp.c src/thpool.c

## Differences among CPPFLAGS, AM_CPPFLAGS, and mumble_CPPFLAGS
## (rules are similar for CFLAGS, CPPFLAGS, and LDFLAGS etc.)
## Refer to 
##   https://www.gnu.org/software/automake/manual/html_node/Flag-Variables-Ordering.html
## Summary
##   Automake always uses two of these variables when compiling C sources files. 
##   When compiling an object file for the mumble target, the first variable 
##   will be mumble_CPPFLAGS if it is defined, or AM_CPPFLAGS otherwise. The 
##   second variable is always CPPFLAGS.
## NOTE
##   Better not try defining or modifying CPPFLAGS, CFLAGS, or LDFLAGS etc. in
##   Makefile* so that the end users are allowed to change these variables 
##   with configure options
cellsnp_lite_CPPFLAGS = -Isrc @HTSLIB_CPPFLAGS@

cellsnp_lite_CFLAGS = -g -Wall -O2 -Wno-unused-function

## The -pthread option is forced to be added to _LDFLAGS by autoreconf
cellsnp_lite_LDFLAGS = @HTSLIB_LDFLAGS@ -pthread

## xxx_LDADD is some additional libraries besides $LIBS which could collect
## the libraries added by AC_CHECK_LIB
cellsnp_lite_LDADD = @HTSLIB_LIB@ 

## Define which files should be removed by `make clean` target
CLEANFILES = *.o a.out src/*.o cellsnp-lite
