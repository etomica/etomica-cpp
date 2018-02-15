#include <math.h>
#include "random.h"

#include "SFMT.h"
#ifdef HAVE_GETRANDOM
#include <sys/random.h>
#else
#include <fcntl.h>
#include <sys/time.h>
#include <unistd.h>
#endif

#define MAXUINT ((uint32_t)-2)

Random::Random() : hasNextGaussian(false) {
  seed = makeSeed();
  sfmt_init_gen_rand(&sfmt, seed);
}

Random::Random(int s) : hasNextGaussian(false) {
  seed = s;
  sfmt_init_gen_rand(&sfmt, seed);
}

int Random::makeSeed() {
  uint32_t seed;
#ifdef HAVE_GETRANDOM
  getrandom(&seed, 4, 0);
#else
  int fd = open("/dev/urandom", O_RDONLY);
  if (fd<0) {
    struct timeval t;
    gettimeofday(&t, NULL);
    unsigned long long x = (unsigned long long)(t.tv_sec) * 1000000 + (unsigned long long)(t.tv_usec);
    seed = (uint32_t)x;
  }
  else {
    ssize_t b = read(fd, &seed, 4);
    if (b!=4) {
      fprintf(stderr, "unable to read from urandom\n");
      exit(1);
    }
    close(fd);
  }
#endif
  return seed;
}

int Random::getSeed() {
  return seed;
}

int Random::nextInt(int max) {
  uint32_t maxRand = MAXUINT - ((MAXUINT+1) % max);
  uint32_t s;
  do {
    s = sfmt_genrand_uint32(&sfmt);
  } while (s > maxRand);
  return s % max;
}

double Random::nextDouble() {
  return sfmt_genrand_res53(&sfmt);
}

double Random::nextDouble32() {
  return sfmt_genrand_real1(&sfmt);
}

double Random::nextGaussian() {
  if (hasNextGaussian) {
    hasNextGaussian = false;
    return nextG;
  }

  double x1, x2, w;
  do {
    x1 = nextDouble();
    x2 = nextDouble();
    w = x1 * x1 + x2 * x2;
  } while (w >= 1);
  w = sqrt(-2 * log(w) / w);
  int signs = nextInt(4);
  if ((signs & 0x00000001) == 1) x1 = -x1;
  if ((signs & 0x00000002) == 2) x2 = -x2;
  nextG = x2 * w;
  hasNextGaussian = true;
  return x1 * w;
}
