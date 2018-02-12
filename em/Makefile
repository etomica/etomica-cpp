CFLAGS = -Os -Wall -I../src/ -I. -Ibuild -DSFMT_MEXP=19937
LFLAGS = -Os -Wall -s WASM=1
WEBIDL = python /usr/data/build/emsdk/emscripten/incoming/tools/webidl_binder.py 
# DEBUG flags
#FLAGS = -O0 -Wall -s ASSERTIONS=2 -s SAFE_HEAP=1 -s STACK_OVERFLOW_CHECK=2 -s DEMANGLE_SUPPORT=1

OBJDIR = build

WEBSTUFF = $(wildcard *.html) $(wildcard *.js) $(wildcard *.css)
EMSRC = $(wildcard *.cpp)
SRC = $(filter-out ../src/main.cpp, $(wildcard ../src/*.cpp))
INC = $(wildcard *.h) $(wildcard ../src/*.h)
IDL = $(wildcard *.idl)
GLUECPP := $(patsubst %.idl,$(OBJDIR)/glue-%.cpp,$(IDL))
GLUEJS := $(patsubst %.idl,$(OBJDIR)/glue-%.js,$(IDL))
POSTJS := $(patsubst %,--post-js %,$(GLUEJS))
EMOBJECTS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(EMSRC))
OBJECTS := $(patsubst ../src/%.cpp,$(OBJDIR)/%.o,$(SRC))

.PHONY: $(OBJDIR)

default: $(OBJDIR) $(OBJDIR)/mc.js

$(OBJDIR)/mc.js: $(EMOBJECTS) $(OBJECTS) $(GLUEOBJECTS) $(GLUEJS)
	em++ $(LFLAGS) $(POSTJS) -o $@ $(EMOBJECTS) $(OBJECTS) $(GLUEOBJECTS)  -lm

$(OBJDIR):
	@mkdir -p $(OBJDIR)
	@for f in $(WEBSTUFF) $(EMSRC) $(SRC) $(INC) $(IDL) Makefile; do \
		ln -s -f ../$$f $(OBJDIR)/; \
	done

$(GLUEOBJECTS): %.o: %.cpp *.h ../src/*.h
	em++ $(CFLAGS) -c -o $@ $<

$(EMOBJECTS): $(OBJDIR)/%.o: %.cpp *.h ../src/*.h $(GLUECPP)
	em++ $(CFLAGS) -c -o $@ $<

$(OBJECTS): $(OBJDIR)/%.o: ../src/%.cpp *.h ../src/*.h
	em++ $(CFLAGS) -c -o $@ $<

$(GLUECPP): $(OBJDIR)/glue-%.cpp: %.idl
	$(WEBIDL) $< $(patsubst %.cpp,%,$@)

$(GLUEJS): $(OBJDIR)/glue-%.js: %.idl
	$(WEBIDL) $< $(patsubst %.js,%,$@)

