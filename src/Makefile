
FLAGS = -O3 -std=c++11 -Wall -msse2 -DHAVE_SSE2

OBJDIR = build

SRC = $(filter-out main%cpp, $(wildcard *.cpp))
MAINSRC = main.cpp
MAINMDSRC = mainMD.cpp
INC = $(wildcard *.h)
OBJECTS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(SRC))
MAINOBJ := $(patsubst %.cpp,$(OBJDIR)/%.o,$(MAINSRC))
MAINMDOBJ := $(patsubst %.cpp,$(OBJDIR)/%.o,$(MAINMDSRC))

default: directories $(OBJDIR)/mc $(OBJDIR)/md

directories:
	@mkdir -p $(OBJDIR)

$(OBJDIR)/mc: $(OBJECTS) $(MAINOBJ)
	g++ $(FLAGS) -o $@ $(OBJECTS) $(MAINOBJ) -lm

$(OBJDIR)/md: $(OBJECTS) $(MAINMDOBJ)
	g++ $(FLAGS) -o $@ $(OBJECTS) $(MAINMDOBJ) -lm

$(OBJECTS) $(MAINOBJ) $(MAINMDOBJ): $(OBJDIR)/%.o: %.cpp *.h
	g++ -fPIC -c $(FLAGS) -DSFMT_MEXP=19937 -o $@ $<
