###
#########################################################################
#
# CIOlib - Cartesian Input / Output library
#
# Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
#########################################################################
###

lib_LIBRARIES = libCIO.a

libCIO_a_CXXFLAGS = @MPI_CFLAGS@ @TP_CFLAGS@ -I$(top_builddir)/include

libCIO_a_SOURCES = \
    cio_ActiveSubDomain.C \
    cio_DFI.C \
    cio_DFI_BOV.C \
    cio_DFI_SPH.C \
    cio_DFI_Write.C \
    cio_TextParser.C

CIO_includedir = $(includedir)
CIO_include_HEADERS = \
    $(top_builddir)/include/cio_ActiveSubDomain.h \
    $(top_builddir)/include/cio_DFI.h \
    $(top_builddir)/include/cio_DFI_BOV.h \
    $(top_builddir)/include/cio_DFI_SPH.h \
    $(top_builddir)/include/cio_Define.h \
    $(top_builddir)/include/cio_TextParser.h \
    $(top_builddir)/include/cio_PathUtil.h \
    $(top_builddir)/include/cio_endianUtil.h

