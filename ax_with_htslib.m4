# modified from github:samtools
# Link: https://github.com/samtools/samtools/blob/develop/m4/ax_with_htslib.m4
#
# ===========================================================================
#      http://www.gnu.org/software/autoconf-archive/ax_with_htslib.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_WITH_HTSLIB
#
# DESCRIPTION
#
#   This macro checks whether HTSlib <http://www.htslib.org/> is installed
#   or nearby, and adds a --with-htslib=DIR option to the configure script
#   for specifying the location.  It locates either an installation prefix
#   (with 'include' and 'lib' subdirectories) or an HTSlib source tree, as
#   HTSlib is fast-moving and users may wish to use an in-development tree.
#
#   Different checks occur depending on the --with-htslib argument given:
#
#   With --with-htslib=DIR, checks whether DIR is a source tree or contains
#     a working installation.
#   By default, searches for a source tree (with a name matching htslib*)
#     within or alongside $srcdir.  Produces AC_MSG_ERROR if there are
#     several equally-likely candidates.  If there are none, checks for
#     a working default installation.
#   With --with-htslib=system, checks for a working default installation.
#   Perfer compling with libhts.a than libhts.so
#
#   The following output variables are set by this macro:
#
#     HTSDIR              Root directory of HTSlib
#     HTSLIB_CPPFLAGS     Preprocessor flags for compiling with HTSlib
#     HTSLIB_LDFLAGS      Linker flags for linking with HTSlib
#     HTSLIB_LIB          HTSLib library file
#
#   The following shell variables may be defined:
#
#     ax_cv_htslib        Set to "yes" if HTSlib was found
#     ax_cv_htslib_which  Set to "source", "system", "install", or "none"
#     ax_cv_htslib_type   Set to "static", "shared", or "none"
#     ax_cv_htslib_msg    Set to <message> for missing or old HTSlib
#
# LICENSE
#
#   Copyright (C) 2015,2017 Genome Research Ltd
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.  This file is offered as-is, without any
#   warranty.

#serial 1

AC_DEFUN([AX_WITH_HTSLIB],
[AC_ARG_WITH([htslib],
  [AS_HELP_STRING([--with-htslib=DIR],
    [use the HTSlib installation in DIR])
dnl Not indented, to avoid extra whitespace outwith AS_HELP_STRING()
AS_HELP_STRING([--with-htslib=system],
    [use only a system HTSlib installation])],
  [], [with_htslib=search])

ax_cv_htslib_type=none

dnl In this step, determine HTSDIR and ax_cv_htslib_which
dnl ax_cv_htslib_which has four possible values:
dnl - source: 
dnl   the installation is in a source tree dir (mostly installed 
dnl   from the github repo)
dnl - system: 
dnl   the installation is system-wide and library should be in the
dnl   system search path
dnl - install:
dnl   the installation is in a dir containing `include` and `lib`
dnl   subdirs
dnl - none:
dnl   when missing or old library
case $with_htslib in
yes|search)
  AC_MSG_CHECKING([location of HTSlib installation])
  case $srcdir in
    .) srcp= ;;
    *) srcp=$srcdir/ ;;
  esac
  found=
  for dir in ${srcp}htslib* -- ${srcp}../htslib -- ${srcp}../htslib*
  do
    if test "$dir" = "--"; then
      test -n "$found" && break
    elif test -f "$dir/libhts.a" || test -f "$dir/libhts.so"; then
      found="${found}1"
      HTSDIR=$dir
    fi
  done
  if test -z "$found"; then
    AC_MSG_RESULT([none found])
    ax_cv_htslib_which=system
    HTSDIR=
  elif test "$found" = 1; then
    AC_MSG_RESULT([$HTSDIR])
    ax_cv_htslib_which=source
  else
    AC_MSG_RESULT([several directories found])
    AC_MSG_ERROR([use --with-htslib=DIR to select which HTSlib to use])
  fi
  ;;
no) ax_cv_htslib_which=none ;;
system) 
  ax_cv_htslib_which=system 
  HTSDIR=
  ;;
*)
  HTSDIR=$with_htslib
  if test -f "$HTSDIR/libhts.a" || test -f "$HTSDIR/libhts.so"; then
    ax_cv_htslib_which=source
  elif test -f "$HTSDIR/lib/libhts.a" || test -f "$HTSDIR/lib/libhts.so"; then
    ax_cv_htslib_which=install
  else
    ax_cv_htslib_which=none
  fi
  ;;
esac

dnl In this step, determine ax_cv_htslib_type and HTSLIB_LIB, HTSLIB_CPPFLAGS, 
dnl HTSLIB_LDFLAGS. The dependencies of htslib.a would be added to $LIBS after
dnl calling AC_CHECK_LIB().
case $ax_cv_htslib_which in
source)
  HTSLIB_CPPFLAGS="-I$HTSDIR"
  HTSLIB_LDFLAGS="-L$HTSDIR"
  if test -f "$HTSDIR/libhts.a"; then
    ax_cv_htslib_type=static
    HTSLIB_LIB=$HTSDIR/libhts.a
  else
    ax_cv_htslib_type=shared
    HTSLIB_LIB="-lhts"
  fi
  ;;
system)
  HTSLIB_CPPFLAGS=
  HTSLIB_LDFLAGS=
  HTSLIB_LIB="-lhts"
  ax_cv_htslib_type=shared
  ;;
install)
  HTSLIB_CPPFLAGS="-I$HTSDIR/include"
  HTSLIB_LDFLAGS="-L$HTSDIR/lib"
  if test -f "$HTSDIR/lib/libhts.a"; then
    ax_cv_htslib_type=static
    HTSLIB_LIB=$HTSDIR/lib/libhts.a
  else
    ax_cv_htslib_type=shared
    HTSLIB_LIB="-lhts"
  fi
  ;;
none)
  ax_cv_htslib=no
  ;;
esac

dnl Check HTSlib and determine ax_cv_htslib, ax_cv_htslib_msg
dnl NOTE: 
dnl - AC_CHECK_LIB tends to? check the shared library (ie. the libhts.so). 
dnl   In this case, it's probably safe to use either libhts.a or libhts.so as 
dnl   the libhts.a and libhts.so are usually built simultaneously, though we
dnl   prefer using libhts.a
dnl - Cellsnp-lite requires HTSlib >= 1.10.2 as it at least uses the 
dnl   KS_INITIALIZE which was added to HTSlib since v1.10. We use function 
dnl   sam_hdr_init, which was added to HTSlib since v1.10, to test version.
if test "$ax_cv_htslib" != no; then
  ax_saved_CPPFLAGS=$CPPFLAGS
  ax_saved_LDFLAGS=$LDFLAGS
  CPPFLAGS="$CPPFLAGS $HTSLIB_CPPFLAGS"
  LDFLAGS="$LDFLAGS $HTSLIB_LDFLAGS"
  AC_CHECK_HEADER([htslib/sam.h],
    [AC_CHECK_LIB(hts, sam_hdr_init, [ax_cv_htslib=yes], 
      [ax_cv_htslib=no
       AC_CHECK_LIB(hts, hts_version,
         [ax_cv_htslib_msg="library is too old (1.10.2+ required)"], 
         [ax_cv_htslib_msg="library not found"])])],
    [ax_cv_htslib=no
     ax_cv_htslib_msg="library not found"], [;])
  CPPFLAGS=$ax_saved_CPPFLAGS
  LDFLAGS=$ax_saved_LDFLAGS
else
  ax_cv_htslib_msg="library not found"
fi

AC_SUBST([HTSDIR])
AC_SUBST([HTSLIB_CPPFLAGS])
AC_SUBST([HTSLIB_LDFLAGS])
AC_SUBST([HTSLIB_LIB])

])

