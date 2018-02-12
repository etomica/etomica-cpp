#pragma once

#include "SFMT.h"
#include "sfmt_more.h"

class Random {
  private:
    int seed;
  public:
    Random();
    Random(int seed);
    sfmt_t sfmt;
    int getSeed();
};
