lib.name = dynsmooth~

class.sources = dynsmooth~.c

datafiles = \
	dynmooth~-help.pd \
	dynamicsmooth~.pd \
	dynamicsmooth-efficient~.pd

objectsdir = ./build
PDLIBBUILDER_DIR=./pd-lib-builder
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
