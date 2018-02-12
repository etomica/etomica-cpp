default: mc

FLAGS = -O0 -g -std=c++11 -Wall -msse2 -DHAVE_SSE2

OBJ = box.o potential.o move.o potential-master.o integrator.o SFMT.o main.o \
			potential-master-cell.o average.o meter-potential-energy.o data-pump.o \
			potential-callback-pressure.o meter-full-compute.o random.o \
			potential-callback-hma.o meter-num-atoms.o meter-density.o

SFMT.o: SFMT.c SFMT.h SFMT-common.h SFMT-params.h SFMT-params19937.h
	g++ -fPIC -c $(FLAGS) -DSFMT_MEXP=19937 -o $@ $<

#	g++ -fPIC -c -O0 -msse2 -g -DHAVE_SSE2 -DSFMT_MEXP=19937 -o $@ $<

%.o: %.cpp *.h
	g++ -fPIC -c $(FLAGS) -DSFMT_MEXP=19937 -o $@ $<

mc: $(OBJ)
	g++ $(FLAGS) -o $@ $(OBJ) -lm

test: $(OBJ) test.o
	g++ $(FLAGS) -o $@ test.o potential.o box.o potential-master.o -lm
