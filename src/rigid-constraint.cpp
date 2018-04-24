#include "rigid-constraint.h"

RigidConstraint::RigidConstraint(vector<int*>& p, vector<int>& a) : pairs(p), atoms(a) {
}

RigidConstraint::~RigidConstraint() {
}

// redistribute forces from implicitly constrained atoms
void RigidConstraint::redistributeForces(Box& box, int iFirstAtom, double** force) {
}

// place implicitly constrained atoms based on constrained ones
void RigidConstraint::relaxMolecules(Box& box, int iFirstAtom) {
}
