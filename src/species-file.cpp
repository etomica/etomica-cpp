#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "species.h"
#include "rigid-constraint.h"

SpeciesFile::SpeciesFile(const char* filename) : Species(0,0), typeOffset(0) {
  FILE* f;
  if (!(f = fopen(filename, "r"))) {
    fprintf(stderr, "Unable to open species file '%s'\n", filename);
    abort();
  }
  readAtomTypes(f, filename);
  vector<double*> tmpPositions;
  char* nextSection = readAtoms(f, filename, tmpPositions);
  while (nextSection) {
    if (!strcmp(nextSection, "orientation")) {
      free(nextSection);
      nextSection = readOrientations(f, filename);
    }
    else if (!strcmp(nextSection, "constraints")) {
      free(nextSection);
      nextSection = readConstraints(f, filename);
    }
  }
  fclose(f);

  // allocate positions, atomTypes
  setup(types.size(), typeSymbols.size());
  // now copy over positions
  for (int i=0; i<numAtoms; i++) {
    double* ri = tmpPositions[i];
    for (int j=0; j<3; j++) positions[i][j] = ri[j];
    delete[] ri;
  }
}

void SpeciesFile::readAtomTypes(FILE* f, const char* filename) {
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
}

char* SpeciesFile::trim(char* s) {
  while (*s && isspace(*s)) s++;
  char* end = s + strlen(s) - 1;
  while (end > s && isspace(*end)) end--;
  *(end+1) = '\0';
  return s;
}

char* SpeciesFile::readAtoms(FILE* f, const char* filename, vector<double*> &tmpPositions) {
  char buf[512];
  char* nextSection = nullptr;
  while (true) {
    if (!fgets(buf, 511, f)) {
      break;
    }
    if (buf[0] == '#') continue; // comment
    char* c = trim(buf);
    if (c[strlen(c)-1] == ':') {
      c[strlen(c)-1] = '\0';
      nextSection = strdup(c);
      break;
    }
    char* atomNumStr = strsep(&c, " ");
    if (!atomNumStr || atomNumStr[0]=='#') continue;
    int atomNum = atoi(atomNumStr);
    if (atomNum != (int)types.size()+1) {
      fprintf(stderr, "atom %d on line %lu!\n", atomNum, types.size()+1);
      abort();
    }
    if (!c) {
      fprintf(stderr, "Found atom #, but no symbol in file %s\n", filename);
      abort();
    }
    char* symbol = strsep(&c, " ");
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
  return nextSection;
}

char* SpeciesFile::readOrientations(FILE* f, const char* filename) {
  int oCount = 0;
  char* nextSection = nullptr;
  char buf[512];
  while (true) {
    if (!fgets(buf, 510, f)) {
      break;
    }
    char* c = trim(buf);
    if (c[strlen(c)-1] == ':') {
      c[strlen(c)-1] = '\0';
      nextSection = strdup(c);
      break;
    }
    char* oNumStr = strsep(&c, " ");
    if (!oNumStr || oNumStr[0]=='#') continue;
    int oNum = atoi(oNumStr);
    if (oCount+1 != oNum) {
      fprintf(stderr, "orientation %d on line %d in file %s\n", oNum, oCount+1, filename);
      abort();
    }
    if (!c) {
      fprintf(stderr, "Found orientation line, but no atoms in file %s\n", filename);
      abort();
    }
    char* atomStr = strsep(&c, " ");
    if (!atomStr || !c) {
      fprintf(stderr, "Could not find first orienation atom in file %s\n", filename);
      abort();
    }
    int atom1 = atoi(atomStr);
    if (atom1 <= 0 || atom1 > (int)types.size() + 1) {
      fprintf(stderr, "Invalid atom %d for orientation in file %s\n", atom1, filename);
      abort();
    }
    atomStr = strsep(&c, " ");
    if (!atomStr) {
      fprintf(stderr, "Could not find first orienation atom in file %s\n", filename);
      abort();
    }
    int atom2 = atoi(atomStr);
    if (atom2 <= 0 || atom2 > (int)types.size() + 1) {
      fprintf(stderr, "Invalid atom %d for orientation in file %s\n", atom2, filename);
      abort();
    }
    if (atom2 == atom1) {
      fprintf(stderr, "Orientation atoms must be different from file %s\n", filename);
      abort();
    }
    axisAtoms[oCount][0] = atom1-1;
    axisAtoms[oCount][1] = atom2-1;
    oCount++;
  }
  return nextSection;
}

char* SpeciesFile::readConstraints(FILE* f, const char* filename) {
  return nullptr;
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
