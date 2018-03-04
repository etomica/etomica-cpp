
CFLAGS = -DSFMT_MEXP=19937 -DHAVE_SSE2 -msse2
FLAGS = -std=c++11 -Wall
ifdef DEBUG
FLAGS += -O0 -g -DDEBUG
else
FLAGS += -O3
endif

OBJDIR = build

SRC = $(filter-out main%cpp, $(wildcard *.cpp))
MAINSRC = main.cpp
MAINMDSRC = mainMD.cpp
MAINDIMERSRC = mainDimer.cpp
INC = $(wildcard *.h)
OBJECTS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(SRC))
MAINOBJ := $(patsubst %.cpp,$(OBJDIR)/%.o,$(MAINSRC))
MAINMDOBJ := $(patsubst %.cpp,$(OBJDIR)/%.o,$(MAINMDSRC))
MAINDIMEROBJ := $(patsubst %.cpp,$(OBJDIR)/%.o,$(MAINDIMERSRC))

default: directories $(OBJDIR)/mc $(OBJDIR)/md $(OBJDIR)/mcDimer

directories:
	@mkdir -p $(OBJDIR)

$(OBJDIR)/mc: $(OBJECTS) $(MAINOBJ)
	g++ $(FLAGS) -o $@ $(OBJECTS) $(MAINOBJ) -lm

$(OBJDIR)/md: $(OBJECTS) $(MAINMDOBJ)
	g++ $(FLAGS) -o $@ $(OBJECTS) $(MAINMDOBJ) -lm

$(OBJDIR)/mcDimer: $(OBJECTS) $(MAINDIMEROBJ)
	g++ $(FLAGS) -o $@ $(OBJECTS) $(MAINDIMEROBJ) -lm

$(OBJECTS) $(MAINOBJ) $(MAINMDOBJ) $(MAINDIMEROBJ): $(OBJDIR)/%.o: %.cpp *.h
	g++ -fPIC -c $(FLAGS) $(CFLAGS) -o $@ $<
