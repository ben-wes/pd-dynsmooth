lib.name = dynsmooth~

class.sources = dynsmooth~.c

datafiles = \
	dynsmooth~-help.pd \
	dynsmooth-abs~.pd \
	dynsmooth-abs-efficient~.pd

objectsdir = ./build
PDLIBBUILDER_DIR=./pd-lib-builder
include $(PDLIBBUILDER_DIR)/Makefile.pdlibbuilder
