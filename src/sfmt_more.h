#pragma once
#include <stdlib.h>

#if defined(__cplusplus)
extern "C" {
#endif
#include "SFMT.h"
#ifdef HAVE_GETRANDOM
#include <sys/random.h>
#else
#include <fcntl.h>
#include <sys/time.h>
#include <unistd.h>
#endif

#define MAXUINT ((uint32_t)-2)

inline static uint32_t sfmt_genrand_int_max(sfmt_t *sfmt, int max) {
  uint32_t maxRand = MAXUINT - ((MAXUINT+1) % max);
  uint32_t s;
  do {
    s = sfmt_genrand_uint32(sfmt);
  } while (s > maxRand);
  return s % max;
}

inline static uint32_t makeSeed() {
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
#if defined(__cplusplus)
}
#endif
