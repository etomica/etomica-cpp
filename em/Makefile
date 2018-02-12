default: mc.js

# FLAGS = -O3 -std=c++11
FLAGS = -Os -s WASM=1 -Wall
#FLAGS = -O0 -Wall -s ASSERTIONS=2 -s SAFE_HEAP=1 -s STACK_OVERFLOW_CHECK=2 -s DEMANGLE_SUPPORT=1

OBJ = box.o potential.o move.o potential-master.o integrator.o SFMT.o \
			potential-master-cell.o average.o meter-potential-energy.o data-pump.o \
			potential-callback-pressure.o meter-full-compute.o \
			potential-callback-hma.o meter-num-atoms.o meter-density.o \
			random.o glue-wrapper.o

SFMT.o: SFMT.c SFMT.h SFMT-common.h SFMT-params.h SFMT-params19937.h
	emcc -c $(FLAGS) -DSFMT_MEXP=19937 -o $@ $<

#	g++ -fPIC -c -O0 -msse2 -g -DHAVE_SSE2 -DSFMT_MEXP=19937 -o $@ $<

glue-wrapper.o: glue.cpp

%.o: %.cpp *.h
	em++ -c $(FLAGS) -DSFMT_MEXP=19937 -o $@ $<

mc.js: $(OBJ) glue.js
	em++ $(FLAGS) --post-js glue.js -o $@ $(OBJ) -lm

glue.cpp: mc.idl
	python /usr/data/build/emsdk/emscripten/incoming/tools/webidl_binder.py mc.idl glue

glue.js: mc.idl
	python /usr/data/build/emsdk/emscripten/incoming/tools/webidl_binder.py mc.idl glue

