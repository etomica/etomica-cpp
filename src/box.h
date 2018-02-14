#pragma once

class Box {
  private:
    double** positions;
    double** velocities;
    double boxHalf[3];
    double boxSize[3];

    int numAtoms, maxNumAtoms;

  public:
    Box();
    ~Box();

    double* getBoxSize() {return boxSize;}
    int getNumAtoms();
    void setNumAtoms(int numAtoms);
    void boxSizeUpdated();
    double* getAtomPosition(int iAtom);
    double* getAtomVelocity(int iAtom);
    void nearestImage(double *dr);
    void initCoordinates();
    void setBoxSize(double x, double y, double z);
    void enableVelocities();
};
