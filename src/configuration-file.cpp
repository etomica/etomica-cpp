/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


#include "action.h"
#include "box.h"

ConfigurationFile::ConfigurationFile(Box& b, const char* f) : box(b), filename(f) {
  setReplication(1,1,1);
}

ConfigurationFile::~ConfigurationFile() { }

void ConfigurationFile::setReplication(int x, int y, int z) {
  replicates[0] = x;
  replicates[1] = y;
  replicates[2] = z;
}

void ConfigurationFile::go() {
  char buf[512];

  FILE* inp;
  if (!(inp=fopen (filename, "r"))) {
    fprintf (stderr, "unable to open configuration file '%s'\n", filename);
    abort();
  }

  int numAtoms = box.getNumAtoms();
  int nRep = replicates[0]*replicates[1]*replicates[2];
  if (numAtoms%nRep != 0) {
    fprintf(stderr, "Number of atoms must be an even multiple of the number of replicates\n");
    abort();
  }
  int nBasis = numAtoms/nRep;
  const double* bs = box.getBoxSize();
  // This scheme assumes a single species
  for (int iAtom=0; iAtom<nBasis; iAtom++) {
    if (!fgets(buf, 510, inp)) {
      fprintf(stderr, "Not enough atoms in config file.  I only found %d\n", iAtom);
      abort();
    }
    double r[3];
    char* c = buf;
    for (int j=0; j<3; j++) {
      char* c2;
      r[j] = strtod(c, &c2);
      if (c == c2) {
        fprintf(stderr, "Could not parse line %d from config file\n", iAtom+1);
        abort();
      }
      c = c2;
    }
    int irep=0;
    double o[3];
    for (int ix=0; ix<replicates[0]; ix++) {
      // -0.5*bs + ix*bs/nrep  -  0.5*bs/nrep
      o[0] = (-0.5 + (ix-0.5)/replicates[0]) * bs[0];
      for (int iy=0; iy<replicates[0]; iy++) {
        o[1] = (-0.5 + (iy-0.5)/replicates[1]) * bs[1];
        for (int iz=0; iz<replicates[0]; iz++) {
          o[2] = (-0.5 + (iz-0.5)/replicates[2]) * bs[2];
          double *ri = box.getAtomPosition(nBasis*irep + iAtom);
          for (int k=0; k<3; k++) ri[k] = r[k] + o[k];
          box.nearestImage(ri);
          irep++;
        }
      }
    }
  }
  fclose(inp);
}
