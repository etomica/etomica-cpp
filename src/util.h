/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include <sys/time.h>
#include <string.h>
#include <ctype.h>

inline double getTime() {
  struct timeval t;
  gettimeofday(&t, NULL);
  unsigned long long tmus = (unsigned long long)(t.tv_sec) * 1000000 + (unsigned long long)(t.tv_usec);
  return tmus*1e-6;
}

#ifdef FAST_ERFC
// we are overriding erfc from libc.  this is ~6x faster and has adequate accuracy
inline double erfc(double x) {
  double t = 1.0 / (1.0 + 0.3275911 * x);
  return exp(-x * x) * (t * (
        0.254829592 + t * (
          -0.284496736 + t * (
            1.421413741 + t * (
              -1.453152027 + 1.061405429 * t)))));
}
#endif

inline long factorial(int i) {
  long f = 1;
  for (int j=2; j<=i; j++) f *= j;
  return f;
}

inline char* trim(char* s) {
  while (*s && isspace(*s)) s++;
  char* end = s + strlen(s) - 1;
  while (end > s && isspace(*end)) end--;
  *(end+1) = '\0';
  return s;
}

