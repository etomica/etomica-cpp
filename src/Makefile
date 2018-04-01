
CFLAGS = -DSFMT_MEXP=19937 -DHAVE_SSE2 -msse2
FLAGS = -std=c++11 -Wall
ifdef DEBUG
FLAGS += -O0 -g -DDEBUG
else
FLAGS += -O3
endif

OBJDIR = build

SRC = $(filter-out main%cpp, $(wildcard *.cpp))
MAINS := mc md dimer mix virial virial-overlap eam eam-mc
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
