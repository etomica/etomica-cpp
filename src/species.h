#pragma once

#include "atom-info.h"

class Species {
  protected:
    int numAtoms, numAtomTypes;
    int* atomTypes;

  public:
    Species(int numAtoms, int numAtomTypes);
    virtual ~Species();
    // to be called by SpeciesList when Species is added
    virtual void init(AtomInfo& atomInfo) = 0;
    int* getAtomTypes();
    int getNumAtoms();
};

class SpeciesSimple : public Species {
  private:
    double mass;

  public:
    SpeciesSimple(int numAtoms, double mass);
    virtual ~SpeciesSimple() {}
    virtual void init(AtomInfo& atomInfo);
};

class SpeciesList {
  private:
    int nSpecies;
    Species** allSpecies;
    AtomInfo atomInfo;
    bool fixed;

  public:
    SpeciesList();
    ~SpeciesList();
    int size();
    void add(Species* species);
    Species* get(int iSpecies);
    AtomInfo& getAtomInfo();
    bool isPurelyAtomic();
};
