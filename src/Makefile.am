###
###################################################################
#
# CIOlib - Cartesian Input / Output library
#
# Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
##################################################################
###
lib_LIBRARIES = libciolib.a

libciolib_a_CXXFLAGS = @MPI_CFLAGS@ @TP_CFLAGS@ -I$(top_builddir)/include

libciolib_a_SOURCES = \
    cio_ActiveSubDomain.C \
    cio_DFI.C \
    cio_DFI_BOV.C \
    cio_DFI_SPH.C \
    cio_DFI_Write.C \
    cio_TextParser.C

ciolib_includedir = $(includedir)
ciolib_include_HEADERS = \
    $(top_builddir)/include/cio_ActiveSubDomain.h \
    $(top_builddir)/include/cio_DFI.h \
    $(top_builddir)/include/cio_DFI_BOV.h \
    $(top_builddir)/include/cio_DFI_SPH.h \
    $(top_builddir)/include/cio_Define.h \
    $(top_builddir)/include/cio_TextParser.h \
    $(top_builddir)/include/endianUtil.h

