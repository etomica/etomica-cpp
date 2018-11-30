/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include "SFMT.h"

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
    void setSeed(int seed);

    int nextInt(int max);
    double nextDouble();
    double nextDouble32();
    double nextGaussian();
    void onSphere(double *v);
    void inSphere(double *v);
};
