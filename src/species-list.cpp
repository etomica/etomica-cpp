#include <cstddef>
#include <stdlib.h>
#include <stdio.h>
#include "species.h"

SpeciesList::SpeciesList() : nSpecies(0), allSpecies(nullptr) {}

SpeciesList::~SpeciesList() {
  if (allSpecies) free(allSpecies);
}

int SpeciesList::size() const {
  return nSpecies;
}

void SpeciesList::add(Species* s) {
  allSpecies = (Species**)realloc(allSpecies, (nSpecies+1)*sizeof(Species*));
  allSpecies[nSpecies] = s;
  nSpecies++;
  s->init(atomInfo);
}

Species* SpeciesList::get(int i) const {
  return allSpecies[i];
}

AtomInfo& SpeciesList::getAtomInfo() {
  return atomInfo;
}

int SpeciesList::getNumAtomTypes() const {
  return atomInfo.getNumTypes();
}

bool SpeciesList::isPurelyAtomic() const {
  for (int i=0; i<nSpecies; i++) {
    if (allSpecies[i]->getNumAtoms() > 1) return false;
  }
  return true;
}
