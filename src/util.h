#include <sys/time.h>

double getTime() {
  struct timeval t;
  gettimeofday(&t, NULL);
  unsigned long long tmus = (unsigned long long)(t.tv_sec) * 1000000 + (unsigned long long)(t.tv_usec);
  return tmus*1e-6;
}
