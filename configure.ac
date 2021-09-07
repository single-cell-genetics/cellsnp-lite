#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#

dnl Refer to 
dnl - github:htslib https://github.com/samtools/htslib/blob/550c6b1b98414ef60eeb665cbfda9f6350d8907c/configure.ac
dnl - github:samtools https://github.com/samtools/samtools/blob/develop/configure.ac

AC_PREREQ([2.63])
AC_INIT([cellsnp-lite], [1.2.1], [hxj5@hku.hk])

dnl Use subdir-objects option so that the source files in subdir `src`
dnl could be compiled successfully
AM_INIT_AUTOMAKE([subdir-objects foreign])

AC_CONFIG_SRCDIR([src/cellsnp.c])

dnl As cellsnp-lite has a header file named config.h, conf.h is used
dnl instead of config.h to store the outputed macros
AC_CONFIG_HEADERS([conf.h])

m4_include([ax_with_htslib.m4])

# Checks for programs.
AC_PROG_AWK
AC_PROG_CC

# Checks for libraries.
dnl Several dependencies of HTSlib & cellsnp-lite would be checked.
dnl The default action-if-found of AC_CHECK_LIB is to add -l<library> to $LIBS.

dnl Check libm
dnl libm is needed by both HTSlib and cellsnp-lite
libm_devel=ok
AC_CHECK_HEADER([math.h], [], [libm_devel=missing], [;])
AC_CHECK_LIB(m, log,  [], [libm_devel=missing])
if test $libm_devel != ok; then
  AC_MSG_ERROR([log() not found
Cellsnp-lite requires a working floating-point math library.

FAILED.  This error must be resolved in order to build cellsnp-lite successfully.])
fi

dnl Check zlib
dnl zlib is needed by both HTSlib and cellsnp-lite
zlib_devel=ok
AC_CHECK_HEADER([zlib.h], [], [zlib_devel=missing], [;])
AC_CHECK_LIB(z, inflate,  [], [zlib_devel=missing])
if test $zlib_devel != ok; then
  AC_MSG_ERROR([zlib development files not found

Cellsnp-lite uses compression routines from the zlib library <http://zlib.net>.
Building HTSlib requires zlib development files to be installed on the build
machine; you may need to ensure a package such as zlib1g-dev (on Debian or
Ubuntu Linux) or zlib-devel (on RPM-based Linux distributions or Cygwin)
is installed.

FAILED.  This error must be resolved in order to build cellsnp-lite successfully.])
fi

dnl Check libbz2
dnl libbz2 is needed by HTSlib when statically linking HTSlib
bz2_devel=ok
AC_CHECK_HEADER([bzlib.h], [], [bz2_devel=missing], [;])
AC_CHECK_LIB([bz2], [BZ2_bzBuffToBuffCompress], [], [bz2_devel=missing])
if test $bz2_devel != ok; then
  AC_MSG_ERROR([libbzip2 development files not found

The CRAM format may use bzip2 compression, which is implemented in HTSlib
by using compression routines from libbzip2 <http://www.bzip.org/>.
Linking HTSlib requires libbzip2 development files to be installed on the
build machine; you may need to ensure a package such as libbz2-dev (on Debian
or Ubuntu Linux) or bzip2-devel (on RPM-based Linux distributions or Cygwin)
is installed.

FAILED.  Please resolve this error to build cellsnp-lite.])
fi

dnl Check liblzma
dnl liblzma is needed by HTSlib when statically linking HTSlib
lzma_devel=ok
AC_CHECK_HEADERS([lzma.h], [], [lzma_devel=header-missing], [;])
AC_CHECK_LIB([lzma], [lzma_easy_buffer_encode], [], [lzma_devel=missing])
if test $lzma_devel != ok; then
  AC_MSG_ERROR([liblzma development files not found

The CRAM format may use LZMA2 compression, which is implemented in HTSlib
by using compression routines from liblzma <http://tukaani.org/xz/>.
Linking HTSlib requires liblzma development files to be installed on the
build machine; you may need to ensure a package such as liblzma-dev (on Debian
or Ubuntu Linux), xz-devel (on RPM-based Linux distributions or Cygwin), or
xz (via Homebrew on macOS) is installed; or build XZ Utils from source.

FAILED.  Please resolve this error to build cellsnp-lite.])
fi

dnl Check libcurl
dnl libcurl is needed by HTSlib when statically linking HTSlib
curl_devel=ok
AC_CHECK_LIB([curl], [curl_easy_pause], [],
  [curl_devel=missing
   AC_CHECK_LIB([curl], [curl_easy_init],
     [message="library is too old (7.18+ required)"],
     [message="library not found"])])
if test $curl_devel != ok; then
  AC_MSG_ERROR([libcurl $message

Support for HTTPS and other SSL-based URLs requires routines from the libcurl
library <http://curl.se/libcurl/>.  Linking HTSlib requires libcurl 
development files to be installed on the build machine; you may need to ensure 
a package such as libcurl4-{gnutls,nss,openssl}-dev (on Debian or Ubuntu Linux) 
or libcurl-devel (on RPM-based Linux distributions or Cygwin) is installed.

FAILED.  Please resolve this error to build cellsnp-lite.])
fi

dnl Check libcrypto
dnl libcrypto is needed by HTSlib when statically linking HTSlib
crypto_devel=ok
AC_CHECK_LIB([crypto], [HMAC], [], [crypto_devel=missing])
if test $crypto_devel != ok; then
  AC_MSG_ERROR([SSL development files not found

Support for AWS S3 URLs requires routines from an SSL library.  Building
HTSlib with libcurl enabled requires SSL development files to be installed
on the build machine; you may need to ensure a package such as libgnutls-dev,
libnss3-dev, or libssl-dev (on Debian or Ubuntu Linux, corresponding to the
libcurl4-*-dev package installed), or openssl-devel (on RPM-based Linux
distributions or Cygwin) is installed.

FAILED.  Please resolve this error to build cellsnp-lite.])
fi

dnl Check HTSlib
AX_WITH_HTSLIB
if test "$ax_cv_htslib" != yes; then
  AC_MSG_ERROR([HTSlib $ax_cv_htslib_msg

Cellsnp-lite uses HTSlib to parse bioinformatics file formats etc. Building it
requires a previously-installed HTSlib. You may need to configure 
--with-htslib=DIR to locate the appropriate HTSlib.

FAILED.  You must supply an HTSlib in order to build cellsnp-lite successfully.])
fi

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT8_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset mkdir pow strdup])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT

