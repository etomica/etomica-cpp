#pragma once

#include "atom-info.h"

class Species {
  protected:
    int numAtoms, numAtomTypes;
    int* atomTypes;
    double** positions;

  public:
    Species(int numAtoms, int numAtomTypes);
    virtual ~Species();
    // to be called by SpeciesList when Species is added
    virtual void init(AtomInfo& atomInfo) = 0;
    int* getAtomTypes();
    int getNumAtoms();
    double* getAtomPosition(int iAtom);
};

class SpeciesSimple : public Species {
  private:
    double mass;

  public:
    SpeciesSimple(int numAtoms, double mass);
    virtual ~SpeciesSimple() {}
    virtual void init(AtomInfo& atomInfo);
    void setAtomPosition(int iAtom, double* iPosition);
};

class SpeciesList {
  private:
    int nSpecies;
    Species** allSpecies;
    AtomInfo atomInfo;

  public:
    SpeciesList();
    ~SpeciesList();
    int size() const;
    void add(Species* species);
    Species* get(int iSpecies) const;
    AtomInfo& getAtomInfo();
    int getNumAtomTypes() const;
    bool isPurelyAtomic() const;
};
