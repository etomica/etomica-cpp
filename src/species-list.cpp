#include <cstddef>
#include <stdlib.h>
#include <stdio.h>
#include "species.h"

SpeciesList::SpeciesList() : nSpecies(0), allSpecies(nullptr) {}

SpeciesList::~SpeciesList() {
  for (int i=0; i<nSpecies; i++) {
    delete allSpecies[i];
  }
  if (allSpecies) free(allSpecies);
}

int SpeciesList::size() const {
  return nSpecies;
}

/**
 * Add Species s to the list.  Species will be destroyed along with the SpeciesList.
 */
int SpeciesList::add(Species* s) {
  allSpecies = (Species**)realloc(allSpecies, (nSpecies+1)*sizeof(Species*));
  allSpecies[nSpecies] = s;
  nSpecies++;
  s->init(atomInfo);
  return nSpecies-1;
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
