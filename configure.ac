###
################################################################
#
# CIO - Cartesian Input / Output library
#
# Copyright (c) RIKEN AICS, Japan. All right reserved. 2013
#
################################################################
###

#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.57])
AC_INIT([CIOlib],[1.0],[keno@riken.jp],[CIOlib])
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/cio_DFI.C])
AM_CONFIG_HEADER([config.h])

# revision No.
CIO_REVISION="20130606_1200"
AC_SUBST(CIO_REVISION)

#
# default prefix setting
#
AC_PREFIX_DEFAULT([/usr/local/CIOlib])
if test "$prefix" = "NONE" ; then
  CIO_INST_DIR=/usr/local/CIOlib
else
  CIO_INST_DIR="$prefix"
fi

#
# MPI environment
#
AC_SUBST(MPI_DIR)
AC_SUBST(MPI_CFLAGS)
AC_ARG_WITH(mpi, [AC_HELP_STRING([--with-mpi=dir],[Specify MPI install directory])], , with_mpi=none)
if test "$with_mpi" != "none" ; then
  MPI_DIR=$with_mpi;
  MPI_CFLAGS=-I$MPI_DIR/include
else
  echo "  Warning: MPI library not specified."
  echo "           if compilation error about MPI occurs, set --with-mpi option"
fi

#
# TextParser environment
#
AC_ARG_WITH(parser, [AC_HELP_STRING([--with-parser=dir],[Specify TextParser install directory])], , with_parser=none)
AC_SUBST(TP_DIR)
AC_SUBST(TP_CFLAGS)
if test "$with_parser" != "none" ; then
  TP_DIR=$with_parser;
  TP_CFLAGS=-I$TP_DIR/include
else
  echo "  Error: TextParser library not specified."
  echo "         set --with-parser option"
  exit 1
fi

#
# CIO special flags
#
CIO_CFLAGS="-I$CIO_INST_DIR/include"
AC_SUBST(CIO_CFLAGS)
CIO_LDFLAGS="-L$CIO_INST_DIR/lib"
AC_SUBST(CIO_LDFLAGS)
CIO_LIBS="-lciolib"
AC_SUBST(CIO_LIBS)
CIO_CXX="$CXX"
AC_SUBST(CIO_CXX)

#
# MAKE_SUB_DIRS
#
AC_SUBST(MAKE_SUB_DIRS)
MAKE_SUB_DIRS="src"

#
# Checks for programs.
#
AC_PROG_CXX
AC_PROG_RANLIB

#
# Checks for header files.
#
AC_CHECK_HEADERS([stdlib.h unistd.h])

#
# Checks for typedefs, structures, and compiler characteristics.
#
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

#
# Checks for library functions.
#
AC_CHECK_FUNCS([gethostname memset mkdir strcasecmp strerror])

AC_CONFIG_FILES([Makefile \
                 src/Makefile \
                 cio-config \
                ])

AC_OUTPUT

chmod +x ./cio-config

