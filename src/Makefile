# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

CFLAGS = -DSFMT_MEXP=19937 -DHAVE_SSE2 -msse2
FLAGS = -std=c++11 -Wall
ifdef DEBUG
FLAGS += -O0 -g -DDEBUG
ifdef VALGRIND
CFLAGS += -DVALGRIND_CHECKS
endif
else
ifdef PROFILE
FLAGS += -O1 -pg
else
FLAGS += -O3
endif
endif
ifdef FAST_ERFC
CFLAGS += -DFAST_ERFC
endif
ifdef LATTICE_DYNAMICS
CFLAGS += -DLATTICE_DYNAMICS -I/usr/include/eigen3
endif

OBJDIR = build

SRC = $(filter-out main%cpp, $(wildcard *.cpp))
MAINS := mc md dimer mix virial virial-overlap eam eam-mc water-mc water-min
MAINSRC = $(patsubst %,main-%.cpp,$(MAINS))
INC = $(wildcard *.h)
OBJECTS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(SRC))
MAINOBJ := $(patsubst %.cpp,$(OBJDIR)/%.o,$(MAINSRC))
MAINEXE := $(patsubst main-%.cpp,$(OBJDIR)/%,$(MAINSRC))

default: directories $(MAINEXE)

directories:
	@mkdir -p $(OBJDIR)

$(MAINEXE): $(OBJDIR)/% : $(OBJDIR)/main-%.o $(OBJECTS) 
	g++ $(FLAGS) -o $@ $(OBJECTS) $< -lm

$(OBJDIR)/%.o: %.cpp *.h
	g++ -fPIC -c $(FLAGS) $(CFLAGS) -o $@ $<
