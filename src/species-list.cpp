#include <cstddef>
#include <stdlib.h>
#include <stdio.h>
#include "species.h"

SpeciesList::SpeciesList() : nSpecies(0), allSpecies(nullptr), fixed(false) {}

SpeciesList::~SpeciesList() {
  if (allSpecies) free(allSpecies);
}

int SpeciesList::size() {
  return nSpecies;
}

void SpeciesList::add(Species* s) {
  if (fixed) {
    fprintf(stderr, "SpeciesList is already in use and cannot be added to\n");
    abort();
  }
  allSpecies = (Species**)realloc(allSpecies, (nSpecies+1)*sizeof(Species*));
  allSpecies[nSpecies] = s;
  nSpecies++;
  s->init(atomInfo);
}

Species* SpeciesList::get(int i) {
  fixed = true;
  return allSpecies[i];
}

AtomInfo& SpeciesList::getAtomInfo() {
  return atomInfo;
}

bool SpeciesList::isPurelyAtomic() {
  fixed = true;
  for (int i=0; i<nSpecies; i++) {
    if (allSpecies[i]->getNumAtoms() > 1) return false;
  }
  return true;
}
