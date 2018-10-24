/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include <vector>

using namespace std;

class Box;
class Species;

class RigidConstraint {
  protected:
    Species& species;
    const vector<int>& rigidAtoms;
    const vector<double>& bondLengths;
    const vector<int>& extraAtoms;
    vector<double*> extraAtomCoeff;
    bool fullRigid;

    void sumToCOM(Box& box, int iAtom, double* r0, double com[3], double mass, double &totMass);

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
    virtual ~RigidConstraint();

    // redistribute forces from implicitly constrained atoms
    void redistributeForces(Box& box, int iFirstAtom, double** force);
    // place implicitly constrained atoms based on constrained ones
    void relaxMolecule(Box& box, int iFirstAtom);

    virtual void adjustPositions(Box& box, int iFirstAtom, double dt) {}

    bool getFullRigid();
    const vector<int>& getRigidAtoms();
    const vector<double>& getBondLengths();
}; 

