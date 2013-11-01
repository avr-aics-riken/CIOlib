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
AC_INIT([CIOlib],[1.4.1],[keno@riken.jp],[CIOlib])
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([src/cio_DFI.C])
AM_CONFIG_HEADER([config.h])

# revision No.
CIO_REVISION="20131101_1600"
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
# platform
#
TGT_DEVICE="$($SHELL ./chk-uname)"
AC_SUBST(TGT_DEVICE)


#
# FORTRAN f90 compiler check
#

if test "_$FC" = "_"; then
  if test "_$F90" != "_"; then
    FC=$F90
  fi
fi
if test "_$FCFLAGS" = "_"; then
  if test "_$F90FLAGS" != "_"; then
    FCFLAGS=$F90FLAGS
  fi
fi



#
# compiler
#
AC_ARG_WITH(comp, [AC_HELP_STRING([--with-comp=(INTEL|FJ)],[Specify Compiler type])], , with_comp=none)



#
# REAL type
#
AC_ARG_WITH(precision, [AC_HELP_STRING([--with-precision=(single|double)],[Specify REAL type [single]])], , with_precision=single)
AC_SUBST(REALOPT)
AC_SUBST(FREALOPT)

if test "$with_precision" = "double" ; then
  REALOPT=-D_REAL_IS_DOUBLE_
  if test "$with_comp" = "INTEL" ; then
    FREALOPT="-r8"
  elif test "$with_comp" = "FJ" ; then
    FREALOPT="-CcdRR8"
  fi
else
  REALOPT=
  FREALOPT=
fi



#
# MPI setting. wrappers like mpicxx, mpicc, mpiCC, mpic++, mpiFCCpx are taking 
# into consideration, if not add in case.
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



#
# MPI environment
#

AC_SUBST(MPI_DIR)
AC_SUBST(MPI_CFLAGS)
AC_SUBST(MPI_LDFLAGS)
AC_SUBST(MPI_LIBS)
AC_SUBST(MPI_FCLIBS)


# 
# MPI library
#
AC_ARG_WITH(MPI, [AC_HELP_STRING([--without-MPI],[do not use MPI library])], , with_MPI=yes)
if test "$with_MPI" = "no" ; then
  CXXFLAGS=-D_CIO_WITHOUT_MPI_
  enable_mpi="none"
fi


#
# if wrapper compiler is specified, skip this part
#

if test "$enable_mpi" = "none" ; then

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

if test "$enable_mpi" != "without" ; then
  if test "$enable_mpi" = "none" ; then
    echo "  Warning: MPI library is not specified."
    echo "           If compilation error occurs about MPI, specify --with-ompi option."
  fi
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

CIO_FC="$FC"
AC_SUBST(CIO_FC)

CIO_FC_LD="$FC"
AC_SUBST(CIO_FC_LD)



if test "$with_comp" = "FJ"; then
  CIO_FC_LD="$CXX"" --linkfortran"
elif test "$with_comp" = "INTEL" ; then
  CIO_LIBS="$CIO_LIBS"" -lstdc++"
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
# MAKE_SUB_DIRS
#
AC_SUBST(MAKE_SUB_DIRS)
MAKE_SUB_DIRS="src doc"

if test "$host" != "sparc*" ; then
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

if test "$host" != "sparc*" ; then
  AC_CONFIG_FILES([tools/frm/src/Makefile])
fi

AC_OUTPUT

chmod +x ./cio-config

