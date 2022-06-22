# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.

-include config.mk

# DEBUG flags
#FLAGS = -O0 -std=c++11 -Wall
FLAGS = -std=c++11 -Wall
ifdef DEBUG
FLAGS += -O0 -g
else
FLAGS += -Os
endif
CFLAGS = $(FLAGS) -I../src/ -I. -Ibuild -DSFMT_MEXP=19937
ifdef LATTICE_DYNAMICS
CFLAGS += -DLATTICE_DYNAMICS -I/usr/include/eigen3
ALLIDL = mc.idl ld.idl
else
ALLIDL = mc.idl
endif
LFLAGS = $(FLAGS) -s WASM=1
#LFLAGS = $(FLAGS) -s WASM=1 -s ASSERTIONS=2 -s SAFE_HEAP=1 -s STACK_OVERFLOW_CHECK=2 -s DEMANGLE_SUPPORT=1
# define WEBIDL or create webidl.mk
ifndef WEBIDL
  # webidl.mk needs to include a line like:
  #WEBIDL = python /usr/data/build/emsdk/emscripten/incoming/tools/webidl_binder.py 
	include webidl.mk
endif

OBJDIR = build

WEBSTUFF = $(wildcard *.php) $(wildcard *.html) $(wildcard *.js) $(wildcard *.css)
EMSRC = $(wildcard *.cpp)
SRC = $(filter-out ../src/main%.cpp, $(wildcard ../src/*.cpp))
INC = $(wildcard *.h) $(wildcard ../src/*.h)
IDL = mc.idl
GLUECPP := $(patsubst %.idl,$(OBJDIR)/glue-%.cpp,$(IDL))
GLUEJS := $(patsubst %.idl,$(OBJDIR)/glue-%.js,$(IDL))
POSTJS := $(patsubst %,--post-js %,$(GLUEJS))
EMOBJECTS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(EMSRC))
OBJECTS := $(patsubst ../src/%.cpp,$(OBJDIR)/%.o,$(SRC))

.PHONY: $(OBJDIR)

default: $(OBJDIR) $(OBJDIR)/mc.js

$(OBJDIR)/mc.js: $(EMOBJECTS) $(OBJECTS) $(GLUEOBJECTS) $(GLUEJS)
	em++ $(LFLAGS) $(POSTJS) -o $@ $(EMOBJECTS) $(OBJECTS) $(GLUEOBJECTS)  -lm

$(OBJDIR)/mc.idl: $(ALLIDL)
	cat $(ALLIDL) > $(OBJDIR)/mc.idl

$(OBJDIR):
	@mkdir -p $(OBJDIR)
	@for f in $(WEBSTUFF) $(EMSRC) $(SRC) $(INC) Makefile; do \
		ln -s -f ../$$f $(OBJDIR)/; \
	done

$(GLUEOBJECTS): %.o: %.cpp *.h ../src/*.h
	em++ $(CFLAGS) -c -o $@ $<

$(EMOBJECTS): $(OBJDIR)/%.o: %.cpp *.h ../src/*.h $(GLUECPP)
	em++ $(CFLAGS) -c -o $@ $<

$(OBJECTS): $(OBJDIR)/%.o: ../src/%.cpp *.h ../src/*.h
	em++ $(CFLAGS) -c -o $@ $<

$(GLUECPP): $(OBJDIR)/glue-%.cpp: $(OBJDIR)/%.idl
	$(WEBIDL) $< $(patsubst %.cpp,%,$@)

$(GLUEJS): $(OBJDIR)/glue-%.js: $(OBJDIR)/%.idl
	$(WEBIDL) $< $(patsubst %.js,%,$@)

