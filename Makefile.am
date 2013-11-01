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

SUBDIRS = @MAKE_SUB_DIRS@

bin_SCRIPTS= cio-config
CLEANFILES = $(bin_SCRIPTS)
EXTRA_DIST = cio-config.in LICENSE #config_cio.sh
dist_data_DATA = AUTHORS COPYING ChangeLog README LICENSE

clean-local:
	find ./ -name *.obj | xargs rm -f
	find ./ -name *.exe | xargs rm -f
	find ./ -name *.lib | xargs rm -f

distclean-local: clean-local
	rm -f ./bin/*
	rm -f ./lib/*

