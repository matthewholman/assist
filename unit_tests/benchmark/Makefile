ifndef REB_DIR
ifneq ($(wildcard ../../../rebound/.*),) # Check for REBOUND in default location
REB_DIR=../../../rebound
endif
ifneq ($(wildcard ../../../../rebound/.*),) # Check for ASSIST being inside REBOUND directory
REB_DIR=../../../
endif
endif
ifndef REB_DIR # REBOUND is not in default location and REB_DIR is not set
    $(error ASSIST not in the same directory as REBOUND.  To use a custom location, you Must set the REB_DIR environment variable for the path to your rebound directory, e.g., export REB_DIR=/Users/dtamayo/rebound.)
endif
PROBLEMDIR=$(shell basename `dirname \`pwd\``)"/"$(shell basename `pwd`)

include $(REB_DIR)/src/Makefile.defs

ASSIST_DIR=../../

all: librebound.so libassist.so
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) -I$(ASSIST_DIR)/src/ -I$(REB_DIR)/src/ -Wl,-rpath,./ $(OPT) $(PREDEF) problem.c -L. -lassist -lrebound $(LIB) -o rebound
	@echo ""
	@echo "Problem file compiled successfully."

librebound.so:
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C $(REB_DIR)/src/
	@echo "Creating link for shared library librebound.so ..."
	@-rm -f librebound.so
	@ln -s $(REB_DIR)/src/librebound.so .

libassist.so: librebound.so $(ASSIST_DIR)/src/*.h $(ASSIST_DIR)/src/*.c
	@echo "Compiling shared library libassist.so ..."
	$(MAKE) -C $(ASSIST_DIR)/src/
	@-rm -f libassist.so
	@ln -s $(ASSIST_DIR)/src/libassist.so .

clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	$(MAKE) -C $(REB_DIR)/src/ clean
	@echo "Cleaning up shared library libassist.so ..."
	@-rm -f libassist.so
	$(MAKE) -C $(ASSIST_DIR)/src/ clean
	@echo "Cleaning up local directory ..."
	@-rm -vf rebound
