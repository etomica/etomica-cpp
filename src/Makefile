
FLAGS = -O3 -std=c++11 -Wall -msse2 -DHAVE_SSE2

OBJDIR = build

SRC = $(wildcard *.cpp)
CSRC = $(wildcard *.c)
INC = $(wildcard *.h)
OBJECTS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(SRC))
COBJECTS := $(patsubst %.c,$(OBJDIR)/%.o,$(CSRC))

default: directories $(OBJDIR)/mc

directories:
	@mkdir -p $(OBJDIR)

$(OBJDIR)/mc: directories $(OBJECTS) $(COBJECTS)
	g++ $(FLAGS) -o $@ $(OBJECTS) $(COBJECTS) -lm

$(COBJECTS): $(OBJDIR)/%.o: %.c *.h
	g++ -fPIC -c $(FLAGS) -DSFMT_MEXP=19937 -o $@ $<

$(OBJECTS): $(OBJDIR)/%.o: %.cpp *.h
	g++ -fPIC -c $(FLAGS) -DSFMT_MEXP=19937 -o $@ $<
