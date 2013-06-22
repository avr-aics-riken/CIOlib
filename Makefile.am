###
###################################################################
#
# CIO - Cartesian Input / Output library
#
# Copyright (c) RIKEN AICS, Japan. All right reserved. 2013
#
##################################################################
###

SUBDIRS = @MAKE_SUB_DIRS@

bin_SCRIPTS = cio-config
EXTRA_DIST = cio-config.in LICENSE #config_cio.sh
dist_data_DATA = AUTHORS COPYING ChangeLog README LICENSE

clean-local:
	find ./ -name *.obj | xargs rm -f
	find ./ -name *.exe | xargs rm -f
	find ./ -name *.lib | xargs rm -f

distclean-local: clean-local
	rm -f ./bin/*
	rm -f ./lib/*

