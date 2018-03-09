#pragma once

#include <stdio.h>
#include <stdlib.h>
#include "species.h"

class Box {
  protected:
    double boxHalf[3];
    double boxSize[3];

    double ***positions, ***velocities;

    const int knownNumSpecies;
    int *numAtomsBySpecies, *maxNumAtomsBySpecies, *speciesNumAtoms;
    int *numMoleculesBySpecies, *maxNumMoleculesBySpecies;
    int **firstAtom, **moleculeIdx;
    int **atomTypes;

    SpeciesList &speciesList;

  public:
    Box(SpeciesList &speciesList);
    virtual ~Box();

    SpeciesList& getSpeciesList() {return speciesList;}
    double* getBoxSize() {return boxSize;}
    void boxSizeUpdated();
    int getTotalNumMolecules() {
      int s = 0;
      for (int i=0; i<knownNumSpecies; i++) s += numMoleculesBySpecies[i];
      return s;
    }
    inline int getNumMolecules(int iSpecies) { return numMoleculesBySpecies[iSpecies]; }
    int getNumAtoms() {
      int s = 0;
      for (int i=0; i<knownNumSpecies; i++) s += numAtomsBySpecies[i];
      return s;
    }
    double* getAtomPosition(int i) {
      int idx = i, iSpecies = 0;
      for ( ; iSpecies<knownNumSpecies-1 && idx >= maxNumAtomsBySpecies[iSpecies]; iSpecies++) {
        idx -= maxNumAtomsBySpecies[iSpecies];
      }
#ifdef DEBUG
      if (idx>=numAtomsBySpecies[iSpecies]) {
        printf("getAtomPosition oops i %d is more atoms than I have\n", i);
        abort();
      }
#endif
      return positions[iSpecies][idx];
    }
    int getAtomType(int i) {
      int idx = i, iSpecies = 0;
      for ( ; iSpecies<knownNumSpecies-1 && idx >= maxNumAtomsBySpecies[iSpecies]; iSpecies++) {
        idx -= maxNumAtomsBySpecies[iSpecies];
      }
#ifdef DEBUG
      if (idx>=numAtomsBySpecies[iSpecies]) {
        printf("gAT oops i %d is more atoms than I have\n", i);
        abort();
      }
#endif
      return atomTypes[iSpecies][idx];
    }
    int getFirstAtom(int iSpecies, int iMoleculeInSpecies) { return firstAtom[iSpecies][iMoleculeInSpecies]; }
    int getGlobalMoleculeIndex(int iSpecies, int iMoleculeInSpecies);
    void nearestImage(double *dr);
    void initCoordinates();
    void setBoxSize(double x, double y, double z);
    void setNumMolecules(int iSpecies, int numMolecules);
    double* getAtomVelocity(int iAtom);
    void enableVelocities();
    void getAtomInfo(int iMolecule, int iSpecies);
    // gives species index, first and last atom indicies for molecule with global index iMolecule
    void getMoleculeInfo(int iMolecule, int &iSpecies, int &firstAtom, int &lastAtom);
    void getMoleculeInfoAtom(int iAtom, int &iMolecule, int &iSpecies, int &firstAtom);
};
