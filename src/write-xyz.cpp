/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdlib.h>
#include "action.h"
#include "box.h"

void WriteXYZ::go(const char* filename, Box& box, char** symbols, bool append) {
  FILE* f = append ? fopen(filename, "a") : fopen(filename, "w");
  if (!f) {
    fprintf(stderr, "Could not open %s for writing\n", filename);
    abort();
  }
  fprintf(f, "%d\n#\n", box.getNumAtoms());
  char s1[2] = {'A', '\0'};
  for (int iAtom=0; iAtom<box.getNumAtoms(); iAtom++) {
    int iType = box.getAtomType(iAtom);
    const char* sym;
    if (symbols) sym = symbols[iType];
    else {
      s1[0] = 'A' + iType;
      sym = s1;
    }
    double* ri = box.getAtomPosition(iAtom);
    fprintf(f, "%s %f %f %f\n", sym, ri[0], ri[1], ri[2]);
  }
  fclose(f);
}
