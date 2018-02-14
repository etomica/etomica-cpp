#pragma once

#include "SFMT.h"
#include "sfmt_more.h"

class Random {
  private:
    int seed;
    static int makeSeed();
    bool hasNextGaussian;
    double nextG;
  public:
    Random();
    Random(int seed);
    sfmt_t sfmt;
    int getSeed();

    int nextInt(int max);
    double nextDouble();
    double nextDouble32();
    double nextGaussian();
};
