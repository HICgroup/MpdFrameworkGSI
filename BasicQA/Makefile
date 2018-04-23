# V.Kireyeu 05.04.2018

# Used compiler
CC=g++

ROOTCONFIG := root-config
ROOTCINT := rootcint

CFLAGS := $(shell $(ROOTCONFIG) --cflags)
CFLAGS += --std=c++11 -O2 

LDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
LDFLAGS += $(shell $(ROOTCONFIG) --glibs)


# Makefile will proceed all files with .cc suffix
# SOURCES=$(wildcard src/*.cc)
# Output executable has the same name without suffix
# OBJECTS=$(patsubst src/%.cc,%, $(SOURCES))

BIN_SOURCES = basicqa.cxx
BIN_OBJECTS = $(patsubst %.cxx,%, $(BIN_SOURCES))

# Library
LIB_SOURCES  = basicqa.cxx
LIB_OBJECTS  = $(patsubst %.cxx,%.so, $(LIB_SOURCES))
# Dict
DICT_SOURCES = $(patsubst %.cxx,%.h, $(LIB_SOURCES))
DICT_SOURCES += $(patsubst %.Cxx,%.LinkDef.h, $(LIB_SOURCES))
DICT_OBJECTS = $(patsubst %.cxx,%.Dict.cxx, $(LIB_SOURCES))
DICT_HEADER = $(patsubst %.cxx,%.Dict.h, $(LIB_SOURCES))

all: $(DICT_OBJECTS) $(LIB_OBJECTS) $(BIN_OBJECTS)

$(DICT_OBJECTS):
	@echo "Generating dictionary..."
	$(ROOTCINT) -f $@ -c -inlineInputHeader $(INCL) $(DICT_SOURCES)

$(LIB_OBJECTS):
	@echo "Generating library..."
	$(CC) -shared -o $@ -fPIC $(CFLAGS) $(INCL) $(LDFLAGS) $(DICT_OBJECTS)

$(BIN_OBJECTS):	$(BIN_SOURCES)
	@echo "Generating executable..."
	$(CC) $(CFLAGS) $(INCL) $(LDFLAGS) -Wl,-rpath=. $(LIB_OBJECTS) $< -o $@

clean:
	if [ -f $(DICT_HEADER) ]; then rm -vf $(DICT_OBJECTS) $(DICT_HEADER) ; fi
	if [ -f $(LIB_OBJECTS) ]; then rm -vf $(LIB_OBJECTS) ; fi
	if [ -f $(BIN_OBJECTS) ]; then rm -vf $(BIN_OBJECTS) ; fi

