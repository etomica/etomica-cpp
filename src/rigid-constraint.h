#pragma once

#include <vector>

using namespace std;

class Box;

class RigidConstraint {
  public:
    vector<int*> pairs;
    vector<int> atoms;

    RigidConstraint(vector<int*>& pairs, vector<int>& atoms);
    ~RigidConstraint();

    // redistribute forces from implicitly constrained atoms
    virtual void redistributeForces(Box& box, int iFirstAtom, double** force);
    // place implicitly constrained atoms based on constrained ones
    virtual void relaxMolecules(Box& box, int iFirstAtom);
}; 

