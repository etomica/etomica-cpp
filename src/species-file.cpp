#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "species.h"

SpeciesFile::SpeciesFile(const char* filename) : Species(0,0), typeOffset(0) {
  FILE* f;
  if (!(f = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open species file '%s'\n", filename);
    abort();
  }
  vector<double*> tmpPositions;
  char buf[512];
  while (true) {
    if (!fgets(buf, 510, f)) {
      fprintf(stderr, "Species file %s ended without any atoms!\n", filename);
      abort();
    }
    if (buf[0] == '#') continue; // comment
    char* c = buf;
    char* symbol = strsep(&c, " ");
    if (!symbol || symbol[0] == '#') continue;
    if (strncmp(symbol, "atoms:", 6) == 0) break;
    if (strncmp(symbol, "types:", 6) == 0) continue;
    if (!c) {
      fprintf(stderr, "found atom type symbol and no mass in file %s\n", filename);
      abort();
    }
    if (strlen(symbol) > 4) {
      fprintf(stderr, "Atom type symbols cannot exceed 4 characters from file %s\n", filename);
      abort();
    }
    for (int i=0; i<(int)strlen(symbol); i++) {
      if (!(symbol[i] >= 'A' && symbol[i] <= 'Z') &&
          !(symbol[i] >= 'a' && symbol[i] <= 'z') &&
          !(symbol[i] >= '0' && symbol[i] <= '9')) {
        fprintf(stderr, "Atom type symbol must be alphnumeric\n");
        abort();
      }
    }
    for (int i=0; i<(int)typeSymbols.size(); i++) {
      if (strcmp(typeSymbols[i],symbol) == 0) {
        fprintf(stderr, "Symbol %s listed twice! in %s\n", symbol, filename);
        abort();
      }
    }
    typeSymbols.push_back(strdup(symbol));
    double m = atof(c);
    typeMass.push_back(m);
  }
  if (typeSymbols.size() == 0) {
    fprintf(stderr, "Did not find any atom types in file %s\n", filename);
    abort();
  }
  // we read the types and encountered an "atoms:" line
  while (true) {
    if (!fgets(buf, 510, f)) {
      break;
    }
    if (buf[0] == '#') continue; // comment
    char* c = buf;
    char* symbol = strsep(&c, " ");
    if (!symbol || symbol[0]=='#') continue;
    int iType = -1;
    for (int i=0; i<(int)typeSymbols.size(); i++) {
      if (strcmp(typeSymbols[i], symbol) == 0) {
        iType = i;
        break;
      }
    }
    if (iType == -1) {
      fprintf(stderr, "Type symbol %s not matched in file %s\n", symbol, filename);
      abort();
    }
    types.push_back(iType);
    double* ri = new double[3];
    for (int j=0; j<3; j++) {
      char* c2;
      ri[j] = strtod(c, &c2);
      if (c==c2) {
        fprintf(stderr, "Could not find coordinates from species file %s\n", filename);
        abort();
      }
      c = c2;
    }
    tmpPositions.push_back(ri);
  }
  if (types.size() == 0) {
    fprintf(stderr, "Did not find any atoms in file %s\n", filename);
    abort();
  }
  fclose(f);

  setup(types.size(), typeSymbols.size());
  // now copy over positions
  for (int i=0; i<numAtoms; i++) {
    double* ri = tmpPositions[i];
    for (int j=0; j<3; j++) positions[i][j] = ri[j];
    delete[] ri;
  }
}

SpeciesFile::~SpeciesFile() {
  for (int i=0; i<numAtomTypes; i++) {
    char* s = typeSymbols[i];
    free(s);
  }
}

void SpeciesFile::init(AtomInfo& ai) {
  Species::init(ai);
  typeOffset = 0;
  for (int i=0; i<numAtomTypes; i++) {
    int myType = ai.addAtomType(typeMass[i]);
    typeOffset = myType-i;
  }
  for (int i=0; i<numAtoms; i++) {
    atomTypes[i] = types[i] + typeOffset;
  }
}


int SpeciesFile::getTypeForSymbol(const char* symbol) {
  for (int i=0; i<numAtomTypes; i++) {
    if (strcmp(typeSymbols[i], symbol) == 0) return i + typeOffset;
  }
  fprintf(stderr, "Unknown symbol %s\n", symbol);
  abort();
}
