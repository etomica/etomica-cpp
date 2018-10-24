/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rigid-constraint.h"
#include "species.h"
#include "rotation-matrix.h"
#include "matrix.h"
#include "vector.h"
#include "box.h"
#include "alloc2d.h"

RigidConstraint::RigidConstraint(Species& s, const vector<int>& ra, const vector<double>& bl, const vector<int>& ea) : species(s), rigidAtoms(ra), bondLengths(bl), extraAtoms(ea) {
  if (rigidAtoms.size() < 2) {
   fprintf(stderr, "I need at least 2 rigid atoms\n"); 
   abort();
  }
  fullRigid = false;
  int bls = bondLengths.size();
  if (rigidAtoms.size() == 2) {
    if (bls != 1) {
      fprintf(stderr, "I need 1 bond length for 2 rigid atoms\n");
      abort();
    }
    fullRigid = true;
  }
  if (rigidAtoms.size() == 3) {
    if (bls < 2 || bls > 3) {
      fprintf(stderr, "I need 2 (chain) or 3 (ring) bond lengths for 3 rigid atoms\n");
      abort();
    }
    fullRigid = bondLengths.size() == 3;
  }
  else if (rigidAtoms.size() == 4) {
    if (bls < 3 || bls == 5 || bls > 6) {
      fprintf(stderr, "I need 3 (chain), 4 (ring) or 6 (tetrahedron) for 4 rigid atoms\n");
      abort();
    }
    fullRigid = bls == 6;
  }
  else {
    if (bls < (int)rigidAtoms.size()-1 || bls > (int)rigidAtoms.size()) {
      fprintf(stderr, "For %d rigid atoms, I need %d (chain) or %d (ring) atoms\n", (int)rigidAtoms.size(), (int)rigidAtoms.size()-1, (int)rigidAtoms.size());
      abort();
    }
    fullRigid = false;
  }
  if (extraAtoms.size() > 0) {
    if (!fullRigid) {
      fprintf(stderr, "Extra atoms are only possible if rigid atoms are fully rigid\n");
      abort();
    }
    for (int i=0; i<(int)extraAtoms.size(); i++) {
      if (species.getMass(extraAtoms[i]) != 0) {
        fprintf(stderr, "Extra atoms in rigid constraint must have no mass!\n");
        abort();
      }
    }
      
    double dr[rigidAtoms.size()-1][3];
    double* r0 = species.getAtomPosition(rigidAtoms[0]);
    for (int i=1; i<(int)rigidAtoms.size(); i++) {
      double* ri = species.getAtomPosition(rigidAtoms[i]);
      for (int k=0; k<3; k++) dr[i-1][k] = ri[k] - r0[k];
    }
    if (rigidAtoms.size() == 2) {
      // linear
      for (int i=0; i<(int)extraAtoms.size(); i++) {
        double* ri = species.getAtomPosition(extraAtoms[i]);
        double c = 0;
        double drMax = 0;
        for (int k=0; k<3; k++) {
          double drik = ri[k] - r0[k];
          if (fabs(dr[0][k]) > 1e-9) {
            double cNew = drik / dr[0][k];
            if (c && fabs(cNew-c) > 1e-9) {
              fprintf(stderr, "Oops, extra atom doesn't actually line up with rigid pair\n");
              abort();
            }
            if (dr[0][k] > drMax) {
              c = cNew;
            }
          }
          else if (fabs(drik) > 1e-6) {
            fprintf(stderr, "Oops, extra atom doesn't actually line up with rigid pair\n");
            abort();
          }
        }
        double* coeff = new double[2];
        coeff[0] = 1-c;  coeff[1] = c;
        extraAtomCoeff.push_back(coeff);
      }
    }
    else if (rigidAtoms.size() == 3) {
      // planar
      RotationMatrix mat;
      mat.setRows(r0, dr[0], dr[1]);
      mat.invert();
      for (int i=0; i<(int)extraAtoms.size(); i++) {
        double* ri = species.getAtomPosition(extraAtoms[i]);
        double dri[3];
        for (int k=0; k<3; k++) {
          dri[k] = ri[k] - r0[k];
        }
        mat.transform(dri);
        if (fabs(dri[0]/dri[1]) > 1e-9 || fabs(dri[0]/dri[2]) > 1e-9) {
          fprintf(stderr, "Extra atom does not seem to be coplanar with rigid atoms\n");
          abort();
        }
        double* coeff = new double[3];
        coeff[0] = 1-dri[1]-dri[2];  coeff[1] = dri[1];  coeff[2] = dri[2];
        extraAtomCoeff.push_back(coeff);
      }
    }
    else {
      // full 3-D
      RotationMatrix mat;
      mat.setRows(dr[0], dr[1], dr[2]);
      mat.invert();
      for (int i=0; i<(int)extraAtoms.size(); i++) {
        double* ri = species.getAtomPosition(extraAtoms[i]);
        double dri[3];
        for (int k=0; k<3; k++) {
          dri[k] = ri[k] - r0[k];
        }
        mat.transform(dri);
        double* coeff = new double[4];
        coeff[0] = 1-dri[0]-dri[1]-dri[2];  coeff[1] = dri[0];  coeff[2] = dri[1];  coeff[3] = dri[2];
        extraAtomCoeff.push_back(coeff);
      }
    }
  }
}

RigidConstraint::~RigidConstraint() {
}

void RigidConstraint::sumToCOM(Box& box, int iAtom, double* r0, double com[3], double mass, double &totMass) {
  if (mass==0) return;
  double *ri = box.getAtomPosition(iAtom);
  totMass += mass;
  if (totMass == mass) {
    for (int k=0; k<3; k++) com[k] = mass*ri[k];
  }
  else {
    double dr[3];
    for (int k=0; k<3; k++) dr[k] = ri[k] - r0[k];
    box.nearestImage(dr);
    for (int k=0; k<3; k++) com[k] += mass*(r0[k] + dr[k]);
  }
}

#define sgn(x) ((x>0) - (x<0))

// redistribute forces from implicitly constrained atoms
void RigidConstraint::redistributeForces(Box& box, int iFirstAtom, double** force) {
  if (extraAtoms.size()==0) return;
  double com[3] = {0.0,0.0,0.0};
  double totMass = 0;
  double *r0 = box.getAtomPosition(iFirstAtom);
  for (int i=0; i<(int)rigidAtoms.size(); i++) {
    sumToCOM(box, iFirstAtom+rigidAtoms[i], r0, com, species.getMass(rigidAtoms[i]), totMass);
  }
  for (int i=0; i<(int)extraAtoms.size(); i++) {
    sumToCOM(box, iFirstAtom+extraAtoms[i], r0, com, species.getMass(extraAtoms[i]), totMass);
  }
  for (int k=0; k<3; k++) com[k] /= totMass;
  double torque[3] = {0.0,0.0,0.0};
  double fnet[3] = {0.0,0.0,0.0};
  double t[3];
  for (int i=0; i<(int)extraAtoms.size(); i++) {
    int iAtom = iFirstAtom + extraAtoms[i];
    for (int k=0; k<3; k++) fnet[k] += force[iAtom][k];
    double* ri = box.getAtomPosition(iAtom);
    double dr[3];
    for (int k=0; k<3; k++) dr[k] = ri[k] - com[k];
    box.nearestImage(dr);
    Vector::cross(dr, force[iAtom], t);
    for (int k=0; k<3; k++) {
      torque[k] += t[k];
      force[iAtom][k] = 0;
    }
  }
  double dr[rigidAtoms.size()][3];
  for (int i=0; i<(int)rigidAtoms.size(); i++) {
    int iAtom = iFirstAtom + rigidAtoms[i];
    double* ri = box.getAtomPosition(iAtom);
    for (int k=0; k<3; k++) {
      force[iAtom][k] += fnet[k] * species.getMass(rigidAtoms[i]) / totMass;
      dr[i][k] = ri[k] - com[k];
    }
    box.nearestImage(dr[i]);
  }
  if (rigidAtoms.size() == 2) {
    double ri0 = sqrt(Vector::dot(dr[0], dr[0]));
    double ri1 = sqrt(Vector::dot(dr[1], dr[1]));
    double fd0[3]; // direction force on 0 needs to be to positively contribute to torque 
    Vector::cross(torque, dr[0], fd0);
    double tMag = 0;
    double fd02 = 0;
    for (int k=0; k<3; k++) {
      tMag += torque[k]*torque[k];
      fd02 += fd0[k]*fd0[k];
    }
    tMag = sqrt(tMag);
    double fd02sqrt = sqrt(fd02);
    for (int k=0; k<3; k++) fd0[k] /= fd02sqrt;
    // take f to be in direction of fd0
    // f ri0 - f sgn(dr[0] dot dr[1]) ri1 = torque
    // f = torque / (ri0 - sgn(dr[0] dot dr[1]) ri1)
    // f0 = f, f1 = - sgn(dr[0] dot dr[1]) f
    double s = sgn(Vector::dot(dr[0], dr[1]));
    double x = ri0 - s*ri1;
    for (int k=0; k<3; k++) {
      force[iFirstAtom+rigidAtoms[0]][k] += fd0[k]*tMag/x;
      force[iFirstAtom+rigidAtoms[1]][k] -= fd0[k]*tMag/x;
    }
  }
  else if (rigidAtoms.size() == 3 || rigidAtoms.size() == 4) {
    // even if we have 4 rigid atoms, redistribute the forces to the first 3
    /**
     * tx = sum[mi (ry Fz - Fy rz)]
     * ty = sum[mi (rz Fx - Fz rx)]
     * tz = sum[mi (rx Fy - Fx ry)]
     *
     * sum[Fx] = 0
     * sum[Fy] = 0
     * sum[Fz] = 0
     *
     * ri dot Fi = 0
     */
    // rows of m correspond to F1x, F1y, F1z, F2x, F2y, F2z, etc.
    Matrix m(9,9);
    double b[9] = {torque[0],torque[1],torque[2],0,0,0,0,0,0};
    for (int i=0; i<3; i++) {
      double* ri = box.getAtomPosition(rigidAtoms[i]);
      m.matrix[0][3*i+1] -= ri[2];
      m.matrix[0][3*i+2] += ri[1];
      m.matrix[1][3*i] += ri[2];
      m.matrix[1][3*i+2] -= ri[0];
      m.matrix[2][3*i] -= ri[1];
      m.matrix[2][3*i+1] += ri[0];
      m.matrix[3][3*i] = 1;
      m.matrix[4][3*i+1] = 1;
      m.matrix[5][3*i+2] = 1;
      m.matrix[6+i][3*i] = ri[0];
      m.matrix[6+i][3*i+1] = ri[1];
      m.matrix[6+i][3*i+2] = ri[2];
    }
    m.transform(b);
    for (int i=0; i<3; i++) {
      double *f = force[rigidAtoms[i]];
      for (int k=0; k<3; k++) f[k] += b[3*i+k];
    }
  }
}

// place implicitly constrained atoms based on constrained ones
void RigidConstraint::relaxMolecule(Box& box, int iFirstAtom) {
  for (int i=0; i<(int)extraAtoms.size(); i++) {
    double* ri = box.getAtomPosition(iFirstAtom + extraAtoms[i]);
    ri[0] = ri[1] = ri[2] = 0;
    double* c = extraAtomCoeff[i];
    for (int j=0; j<(int)rigidAtoms.size(); j++) {
      double* rj = box.getAtomPosition(iFirstAtom + rigidAtoms[j]);
      for (int k=0; k<3; k++) ri[k] += c[j]*rj[k];
    }
  }
}

const vector<int>& RigidConstraint::getRigidAtoms() {
  return rigidAtoms;
}

const vector<double>& RigidConstraint::getBondLengths() {
  return bondLengths;
}
