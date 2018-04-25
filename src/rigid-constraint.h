#pragma once

#include <vector>

using namespace std;

class Box;
class Species;
class RotationMatrix;

class RigidConstraint {
  protected:
    Species& species;
    const vector<int>& rigidAtoms;
    const vector<double>& bondLengths;
    const vector<int>& extraAtoms;
    vector<double*> extraAtomCoeff;
    bool fullRigid;

    void sumToCOM(Box& box, int iAtom, double* r0, double com[3], double mass, double &totMass);
    void sumToMoment(Box& box, int iAtom, double com[3], RotationMatrix moment, double mass);

  public:
    //                     rigidAtoms      bondLengths
    // simple pair             2                1
    // 3-atom chain            3                2
    // 3-atom ring             3                3
    // 4-atom chain            4                3
    // 4-atom ring             4                4
    // tetrahedron             4                6
    // n-atom chain            n               n-1
    // n-atom ring             n                n
    RigidConstraint(Species& species, const vector<int>& rigidAtoms, const vector<double>& bondLengths, const vector<int>& atoms);
    ~RigidConstraint();

    // redistribute forces from implicitly constrained atoms
    virtual void redistributeForces(Box& box, int iFirstAtom, double** force);
    // place implicitly constrained atoms based on constrained ones
    virtual void relaxMolecules(Box& box, int iFirstAtom);
}; 

