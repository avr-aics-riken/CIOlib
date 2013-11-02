###
#########################################################################
#
# CIO - Cartesian Input / Output library
#
# Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN. 
# All rights reserved.
#
#########################################################################
###

#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.57])
AC_INIT([CIOlib],[1.4.2],[keno@riken.jp],[CIOlib])
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/cio_DFI.C])
AM_CONFIG_HEADER([config.h])

# revision No.
CIO_REVISION="20131102_1500"
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
# FORTRAN
#

FC=$F90
FCFLAGS=$F90FLAGS




#
# MPI environment
#
enable_mpi="none"


CIO_CXX="$CXX"
AC_SUBST(CIO_CXX)

case "$CIO_CXX" in
  mpicxx)   enable_mpi="wrappermpi" ;;
  mpicc)    enable_mpi="wrappermpi" ;;
  mpiCC)    enable_mpi="wrappermpi" ;;
  mpic++)   enable_mpi="wrappermpi" ;;
  mpiFCCpx) enable_mpi="wrappermpi" ;;
esac




AC_SUBST(MPI_DIR)
AC_SUBST(MPI_CFLAGS)
AC_SUBST(MPI_LDFLAGS)
AC_SUBST(MPI_LIBS)


# 
# MPI library
#
AC_ARG_WITH(MPI, [AC_HELP_STRING([--without-MPI],[do not use MPI library])], , with_MPI=yes)
if test "$with_MPI" = "no" ; then
  CXXFLAGS=-D_CIO_WITHOUT_MPI_
else

  # openmpi
  AC_ARG_WITH( ompi, [AC_HELP_STRING([--with-ompi=dir],[Specify OpenMPI install directory])], , with_ompi=none)
  if test "$with_ompi" != "none" ; then
    MPI_DIR=$with_ompi;
    MPI_CFLAGS=-I$MPI_DIR/include
    MPI_LDFLAGS=-L$MPI_DIR/lib
    MPI_LIBS=-lmpi
    enable_mpi="ompi"
  fi
fi


#
# Warning
#

if test "$with_MPI" = "yes" ; then
  if test "$enable_mpi" != "wrappermpi" ; then
    if test "$with_ompi" = "none" ; then
      echo "  Warning: OpenMPI library is not specified."
      echo "           If compilation error occurs about MPI, specify --with-ompi option."
    fi
  fi
fi




#
# TextParser environment
#
AC_ARG_WITH(parser, [AC_HELP_STRING([--with-parser=dir],[Specify TextParser install directory])], , with_parser=none)
AC_SUBST(TP_DIR)
AC_SUBST(TP_CFLAGS)
AC_SUBST(TP_LDFLAGS)

if test "$with_parser" != "none" ; then
  TP_DIR=$with_parser;
  TP_CFLAGS=`$TP_DIR/bin/tp-config --cflags`
  TP_LDFLAGS=`$TP_DIR/bin/tp-config --libs`
else
  echo "  Error: TextParser library is not specified."
  echo "         Specify --with-parser option."
  exit 1
fi




#
# CIO special flags
#
CIO_CFLAGS="-I$CIO_INST_DIR/include"
AC_SUBST(CIO_CFLAGS)

CIO_LDFLAGS="-L$CIO_INST_DIR/lib"
AC_SUBST(CIO_LDFLAGS)

CIO_LIBS="-lCIO"
AC_SUBST(CIO_LIBS)



#
# MAKE_SUB_DIRS
#
AC_SUBST(MAKE_SUB_DIRS)
MAKE_SUB_DIRS="src doc"

# K/FX
if test "$host" != "sparc64-unknown-linux-gnu" ; then
  MAKE_SUB_DIRS="$MAKE_SUB_DIRS"" tools/frm/src"
fi



#
# Checks for programs.
#
AC_PROG_CXX
AC_PROG_FC
AC_PROG_RANLIB


#
# Checks for header files.
#
AC_HEADER_STDC
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
                 doc/Makefile \
                 cio-config \
                 include/cio_Version.h \
                 ])

# K/FX
if test "$host" != "sparc64-unknown-linux-gnu" ; then
  AC_CONFIG_FILES([tools/frm/src/Makefile])
fi



AC_OUTPUT

chmod +x ./cio-config

