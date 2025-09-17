# Makefile (drop-in replacement)

CXX      ?= g++
CXXSTD    = -std=c++17
OPTFLAGS ?= -O2
WARN     = -Wall -Wextra -Wshadow -pedantic

# Pythia8 (required)
PY8CXX  := $(shell pythia8-config --cxxflags 2>/dev/null)
PY8LIBS := $(shell pythia8-config --libs --ldflags 2>/dev/null)

ifeq ($(strip $(PY8CXX)),)
  $(error "pythia8-config not found or not usable. Ensure Pythia8 is installed and in PATH.")
endif

# FastJet (optional)
FJCXX  := $(shell which fastjet-config >/dev/null 2>&1 && fastjet-config --cxxflags)
FJLIBS := $(shell which fastjet-config >/dev/null 2>&1 && fastjet-config --libs)

# ROOT (optional)
ROOTCXX  := $(shell which root-config >/dev/null 2>&1 && root-config --cflags)
ROOTLIBS := $(shell which root-config >/dev/null 2>&1 && root-config --libs)


CXXFLAGS :=  $(CXXSTD) $(WARN) -fPIC $(PY8CXX) $(OPTFLAGS) $(FJCXX) $(ROOTCXX)
LDFLAGS  := $(PY8LIBS) $(FJLIBS) $(ROOTLIBS)


all: makeTree

makeTree: makeTree.cc
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) -o $@

clean:
	rm -f makeTree *.o *.root
