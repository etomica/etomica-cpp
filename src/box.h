#pragma once

#include "species.h"

class Box {
  protected:
    double boxHalf[3];
    double boxSize[3];

    double ***positions, ***velocities;

    int knownNumSpecies;
    int *numAtomsBySpecies, *numMoleculesBySpecies, *maxNumMoleculesBySpecies;
    int **firstAtom, **moleculeIdx;

    SpeciesList &speciesList;

  public:
    Box(SpeciesList &speciesList);
    virtual ~Box();

    double* getBoxSize() {return boxSize;}
    void boxSizeUpdated();
    int getNumAtoms();
    int getNumMolecules();
    void nearestImage(double *dr);
    void initCoordinates();
    void setBoxSize(double x, double y, double z);
    void setNumMolecules(int iSpecies, int numMolecules);
    double* getAtomPosition(int iAtom);
    double* getAtomVelocity(int iAtom);
    void enableVelocities();
};
