ifndef REB_DIR
	
ifneq ($(wildcard ../../rebound/.*),) # Check for REBOUND in default location
REB_DIR=../../rebound
endif

ifneq ($(wildcard ../../../rebound/.*),) # Check for ASSIST being inside the REBOUND directory
REB_DIR=../..
endif

endif

# If a REBOUND source checkout is available, use its Makefile.defs.
# Otherwise fall back to the installed Python `rebound` package for `rebound.h`
# and the `librebound*` shared library.
REB_DEFS := $(wildcard $(REB_DIR)/src/Makefile.defs)

ifeq ($(strip $(REB_DEFS)),)
# --------------------------
# Python-package REBOUND mode
# --------------------------
PYTHON ?= ../.venv/bin/python
ifeq ($(wildcard $(PYTHON)),)
PYTHON := python3
endif

CC ?= cc
OPT ?= -std=c99 -Wpointer-arith -D_GNU_SOURCE -O3 -fPIC -Wall -g -Wno-unknown-pragmas -fPIC -DLIBASSIST
PREDEF ?=
LIB ?= -lm

REB_PY_DIR := $(shell $(PYTHON) -c "import rebound, pathlib; print(pathlib.Path(rebound.__file__).resolve().parent)" 2>/dev/null)
REB_LIB_PATH := $(shell $(PYTHON) -c "import rebound, pathlib; site=pathlib.Path(rebound.__file__).resolve().parent.parent; c=[]; [c.extend(sorted(site.glob(pat))) for pat in ('librebound*.so','librebound*.dylib','librebound*.dll')]; print(str(c[0]) if c else '')" 2>/dev/null)
REB_LIB_BASENAME := $(notdir $(REB_LIB_PATH))

# Ensure `@rpath/librebound*` can be resolved when `libassist.so` is loaded via
# Python (dlopen) by adding an rpath relative to this library.
RPATH_FLAGS :=
ifeq ($(shell uname),Darwin)
RPATH_FLAGS += -Wl,-rpath,@loader_path -Wl,-rpath,@loader_path/.
endif
ifeq ($(shell uname),Linux)
RPATH_FLAGS += -Wl,-rpath,$$ORIGIN -Wl,-rpath,$$ORIGIN/.
endif

ifeq ($(strip $(REB_PY_DIR)),)
    $(error Could not import the Python package 'rebound' using PYTHON='$(PYTHON)'. Install it (e.g. 'pip install rebound') or set REB_DIR to a REBOUND source checkout.)
endif
ifeq ($(strip $(REB_LIB_PATH)),)
    $(error Could not locate a 'librebound*' shared library next to the Python package 'rebound' (searched sibling dir of $(REB_PY_DIR)). Install the wheel that bundles the library or set REB_DIR.)
endif

REB_INC := -I$(REB_PY_DIR)

else
# -----------------------
# Source REBOUND checkout
# -----------------------
include $(REB_DIR)/src/Makefile.defs
OPT+= -fPIC -DLIBASSIST
REB_INC := -I$(REB_DIR)/src
REB_LIB_PATH := $(REB_DIR)/src/librebound.so
endif

ifndef ASSISTGITHASH
	ASSISTGITHASH = $(shell git rev-parse HEAD || echo '0000000000gitnotfound0000000000000000000')
	PREDEF+= -DASSISTGITHASH=$(ASSISTGITHASH)
endif

SOURCES=assist.c spk.c ascii_ephem.c forces.c tools.c
OBJECTS=$(SOURCES:.c=.o)
HEADERS=$(SOURCES:.c=.h)

all: $(SOURCES) librebound.so libassist.so
	
%.o: %.c $(HEADERS)
	@echo "Compiling source file $< ..."
	$(CC) -c $(OPT) $(PREDEF) $(REB_INC) -o $@ $<

librebound.so:
	@echo "Creating link for shared library librebound.so ..."
	@-rm -f librebound.so
ifeq ($(strip $(REB_DEFS)),)
	@ln -s $(REB_LIB_PATH) librebound.so
	@-rm -f $(REB_LIB_BASENAME)
	@ln -s $(REB_LIB_PATH) $(REB_LIB_BASENAME)
else
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C $(REB_DIR)/src/
	@ln -s $(REB_LIB_PATH) librebound.so
endif

libassist.so: $(OBJECTS)
	@echo ""        
	@echo "Linking shared library $@ ..."
ifeq ($(strip $(REB_DEFS)),)
	$(CC) $(OPT) -shared $(OBJECTS) $(LIB) -L. -lrebound $(RPATH_FLAGS) -o $@ 
else
	$(CC) $(OPT) -shared $(OBJECTS) $(LIB) -lrebound -L$(REB_DIR)/src -o $@ 
endif
	@echo ""        
	@echo "The shared library $@ has been created successfully."
	
clean:
	@echo "Cleaning up shared library librebound.so ..."
	@-rm -f librebound.so
	@if [ -n "$(REB_DIR)" ] && [ -f "$(REB_DIR)/src/Makefile" ]; then $(MAKE) -C $(REB_DIR)/src/ clean; fi
	@echo "Cleaning up shared library libassist.so ..."
	@-rm -f libassist.so
	@-rm -f *.o
	
