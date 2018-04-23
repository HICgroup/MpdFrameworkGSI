# V.Kireyeu 05.04.2018

# Used compiler
CC=g++

ROOTCONFIG := root-config
ROOTCINT := rootcint
# Path to the MpdRoot installation
MPDROOTDIR := ${VMCWORKDIR}
FAIRINCDIR := $(shell $(ROOTCONFIG) --prefix)/include/


CFLAGS := $(shell $(ROOTCONFIG) --cflags)
CFLAGS += --std=c++11 -O2 
INCL := -I$(MPDROOTDIR)/base/event 
INCL += -I$(MPDROOTDIR)/base/source 
INCL += -I$(MPDROOTDIR)/base/steer 
INCL += -I$(MPDROOTDIR)/base/field 
INCL += -I$(MPDROOTDIR)/base/sim 
INCL += -I$(MPDROOTDIR)/geobase 
INCL += -I$(MPDROOTDIR)/tpc 
INCL += -I$(MPDROOTDIR)/kalman 
INCL += -I$(MPDROOTDIR)/lhetrack 
INCL += -I$(MPDROOTDIR)/mcstack 
INCL += -I$(MPDROOTDIR)/strawendcap 
INCL += -I$(MPDROOTDIR)/etof 
INCL += -I$(MPDROOTDIR)/tof 
INCL += -I$(MPDROOTDIR)/sft 
INCL += -I$(MPDROOTDIR)/sts 
INCL += -I$(MPDROOTDIR)/parbase 
INCL += -I$(MPDROOTDIR)/mpddata 
INCL += -I$(MPDROOTDIR)/mpdbase 
INCL += -I$(MPDROOTDIR)/fairtools 
INCL += -I$(MPDROOTDIR)/clustering 
INCL += -I$(MPDROOTDIR)/geometry 
INCL += -I$(MPDROOTDIR)/mpdpid 
INCL += -I$(MPDROOTDIR)/passive
INCL += -I$(MPDROOTDIR)/cpc
INCL += -I$(MPDROOTDIR)/generators
INCL += -I$(MPDROOTDIR)/mpdfield
INCL += -I$(MPDROOTDIR)/geane
INCL += -I$(MPDROOTDIR)/zdc
INCL += -I$(FAIRINCDIR)

LDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
LDFLAGS += $(shell $(ROOTCONFIG) --glibs)
LDFLAGS += -L$(MPDROOTDIR)/build/lib 
LDFLAGS += -lMpdField -lFairTools -lGeoBase  -lBase -lMCStack -lPassive -lTrkBase  
LDFLAGS += -lMpdBase -lMpdData -lMpdGeneralGenerator  -lKalman -lGeane -lCluster 
LDFLAGS += -ltpc -lTof -lEtof -lStrawendcap  -lZdc -lFfd -lCpc -lStrawECT 
LDFLAGS += -lLHETrack -lSts -lMpdPid 


# Makefile will proceed all files with .cc suffix
# SOURCES=$(wildcard src/*.cc)
# Output executable has the same name without suffix
# OBJECTS=$(patsubst src/%.cc,%, $(SOURCES))

BIN_SOURCES = reducedTreeCreator.C
BIN_OBJECTS = $(patsubst %.C,%, $(BIN_SOURCES))

# Library
LIB_SOURCES  = reducedTreeCreator.C
LIB_OBJECTS  = $(patsubst %.C,%.so, $(LIB_SOURCES))
# Dict
DICT_SOURCES = $(patsubst %.C,%.h, $(LIB_SOURCES))
DICT_SOURCES += $(patsubst %.C,%.LinkDef.h, $(LIB_SOURCES))
DICT_OBJECTS = $(patsubst %.C,%.Dict.cxx, $(LIB_SOURCES))
DICT_HEADER = $(patsubst %.C,%.Dict.h, $(LIB_SOURCES))

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

