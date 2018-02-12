#include "random.h"

Random::Random() {
  seed = makeSeed();
  sfmt_init_gen_rand(&sfmt, seed);
}

Random::Random(int s) {
  seed = s;
  sfmt_init_gen_rand(&sfmt, seed);
}

int Random::getSeed() {
  return seed;
}
