/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "alloc2d.h"

class ArrayUtil {
  private:
    void* p;
  public:
    ArrayUtil(void* pointer) : p(pointer) {}
    ~ArrayUtil() {}
    double x(int i) { return ((double*)p)[i]; }
    double x2d(int i, int j) { return ((double**)p)[i][j]; }
    int ix(int i) { return ((int*)p)[i]; }
    double ix2d(int i, int j) { return ((int**)p)[i][j]; }
    // free a arrays returned by native code
    void freeau() { free(p); }
    void free2Dau() { free2D((void**)p); }
};
