/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include <math.h>
#include "shake-listener.h"
#include "rigid-constraint.h"
#include "species.h"
#include "box.h"
#include "species.h"
#include "alloc2d.h"

ShakeListener::ShakeListener(SpeciesList& sl, Box& b, IntegratorMD& i) : speciesList(sl), box(b), integrator(i), tol(1e-9), maxIterations(500), maxNumAtoms(0), maxNumBonds(0) {
  rOld = rNew = rijOld = rijNew = nullptr;
  rhs = rhsOld = dl = dd = du = dl0 = dd0 = du0 = nullptr;
  lambda = nullptr;
  callStepStarted = callPreForce = callStepFinished = true;

  constraints = new vector<RigidConstraint*>[speciesList.size()];
}

ShakeListener::~ShakeListener() {
  free2D((void**)rOld);
  free2D((void**)rNew);
  free2D((void**)rijOld);
  free2D((void**)rijNew);
  free(dl);
  free(dd);
  free(du);
  free(dl0);
  free(dd0);
  free(du0);
  free(rhs);
  free(rhsOld);
  delete[] constraints;
}

void ShakeListener::init() {
  int nAtoms = 0;
  int maxBonds = 0;
  for (int iSpecies=0; iSpecies<speciesList.size(); iSpecies++) {
    vector<RigidConstraint*> &iConstraints = constraints[iSpecies];
    if (iConstraints.size() == 0) continue;
    int iBonds = 0;
    int iAtoms = 0;
    for (int j=0; j<(int)iConstraints.size(); j++) {
      iBonds += iConstraints[j]->getBondLengths().size();
      int n = iConstraints[j]->getRigidAtoms().size();
      iAtoms += n;
    }
    int nm = box.getNumMolecules(iSpecies);
    nAtoms += nm*iAtoms;
    if (iBonds > maxBonds) maxBonds = iBonds;
  }
  if (nAtoms > maxNumAtoms) {
    // reallocate everything
    rOld = (double**)realloc2D((void**)rOld, nAtoms, 3, sizeof(double));
    maxNumAtoms = nAtoms;
  }
  if (maxBonds > maxNumBonds) {
    // reallocate everything, should only be needed first time
    rNew = (double**)realloc2D((void**)rNew, maxBonds+1, 3, sizeof(double));
    rijOld = (double**)realloc2D((void**)rijOld, maxBonds, 3, sizeof(double));
    rijNew = (double**)realloc2D((void**)rijNew, maxBonds, 3, sizeof(double));
    dl = (double*)realloc(dl, maxBonds*sizeof(double));
    dd = (double*)realloc(dd, maxBonds*sizeof(double));
    du = (double*)realloc(du, maxBonds*sizeof(double));
    rhs = (double*)realloc(rhs, maxBonds*sizeof(double));
    rhsOld = (double*)realloc(rhsOld, maxBonds*sizeof(double));
    lambda = (double*)realloc(lambda, maxBonds*sizeof(double));
    maxNumBonds = maxBonds;
  }
}

void ShakeListener::stepStarted() {
  int oldOffset = 0;
  int iFirstAtom = 0;
  for (int iSpecies=0; iSpecies<speciesList.size(); iSpecies++) {
    vector<RigidConstraint*> &iConstraints = constraints[iSpecies];
    int nm = box.getNumMolecules(iSpecies);
    Species* s = speciesList.get(iSpecies);
    int na = s->getNumAtoms();
    if (iConstraints.size() == 0) {
      iFirstAtom += nm*na;
      continue;
    }
    for (int iMolecule=0; iMolecule<nm; iMolecule++) {
      for (vector<RigidConstraint*>::iterator it = iConstraints.begin(); it!=iConstraints.end(); it++) {
        const vector<int> & rigidAtoms = (*it)->getRigidAtoms();
        for (int i=0; i<(int)rigidAtoms.size(); i++) {
          double* ri = box.getAtomPosition(iFirstAtom+rigidAtoms[i]);
          for (int k=0; k<3; k++) {
            rOld[oldOffset+rigidAtoms[i]][k] = ri[k];
          }
        }
      }
    }
  }
}

void ShakeListener::preForce() {
  double dt = integrator.getTimeStep();
  int iFirstAtom = 0;
  int oldOffset = 0;
  for (int iSpecies=0; iSpecies<speciesList.size(); iSpecies++) {
    vector<RigidConstraint*> &iConstraints = constraints[iSpecies];
    int nm = box.getNumMolecules(iSpecies);
    Species* s = speciesList.get(iSpecies);
    int na = s->getNumAtoms();
    if (iConstraints.size() == 0) {
      iFirstAtom += nm*na;
      continue;
    }
    for (int iMolecule=0; iMolecule<nm; iMolecule++) {
      for (vector<RigidConstraint*>::iterator it = iConstraints.begin(); it!=iConstraints.end(); it++) {
        const vector<double> & bl = (*it)->getBondLengths();
        const vector<int> & rigidAtoms = (*it)->getRigidAtoms();
        for (int i=0; i<(int)rigidAtoms.size(); i++) {
          double* ri = box.getAtomPosition(iFirstAtom+rigidAtoms[i]);
          for (int k=0; k<3; k++) rNew[i][k] = ri[k];
        }
        double* rNewPrev = rNew[iFirstAtom+rigidAtoms[0]];
        double* rOldPrev = rOld[oldOffset+rigidAtoms[0]];
        double maxError = 0;
        for (int i=1; i<(int)rigidAtoms.size(); i++) {
          int b = i-1;
          double* riNew = rNew[i];
          double* riOld = rOld[oldOffset+rigidAtoms[i]];
          if (b>0) dl[b-1] = 0;
          dd[b] = du[b] = 0;
          double rij2 = 0;
          for (int k=0; k<3; k++) {
            rijNew[b][k] = rNewPrev[k] - riNew[k];
            rij2 += rijNew[b][k]*rijNew[b][k];
            rijOld[b][k] = rOldPrev[k] - riOld[k];
            if (b>0) dl[b-1] += -2*rijOld[b-1][k]*rijNew[b][k];
            dd[b] += 2*rijOld[b][k]*rijNew[b][k];
            if (b>1) du[b] += -2*rijOld[b][k]*rijNew[b-1][k];
          }
          rhsOld[b] = rhs[b] = bl[b]*bl[b] - rij2;
          double e = fabs(rhs[b]) / (2*bl[b]*bl[b]);
          if (e>maxError) maxError = e;
          rNewPrev = riNew;
          rOldPrev = riOld;
        }
        const bool ring = rigidAtoms.size() == bl.size() && bl.size() > 2;
        int b = bl.size()-1;
        dl[b-1] = dd[b] = du[b] = dl[b] = du[0] = 0;
        if (ring) {
          // ring
          double* riNew = rNew[0];
          double* riOld = rOld[oldOffset+rigidAtoms[0]];
          double rij2 = 0;
          for (int k=0; k<3; k++) {
            rijNew[b][k] = rNewPrev[k] - riNew[k];
            rij2 += rijNew[b][k]*rijNew[b][k];
            rijOld[b][k] = rOldPrev[k] - riOld[k];
            dl[b-1] += -2*rijOld[b-1][k]*rijNew[b][k];
            dd[b] += 2*rijOld[b][k]*rijNew[b][k];
            du[b] += -2*rijOld[b][k]*rijNew[b-1][k];

            dl[b] += -2*rijOld[b][k]*rijNew[0][k];
            du[0] += -2*rijOld[0][k]*rijNew[b][k];
          }
          rhsOld[b] = rhs[b] = bl[b]*bl[b] - rij2;
          double e = fabs(rhs[b]) / (2*bl[b]*bl[b]);
          if (e>maxError) maxError = e;
        }
        int iter = 0;
        while (maxError >= tol && iter < maxIterations) {
          for (int i=0; i<(int)rigidAtoms.size(); i++) {
            du0[i] = du[i]; dd0[i] = dd[i]; du0[i] = du[i];
          }
          tridiagSolve(bl.size(), dl0, dd0, du0, rhs);

          for (int i=0; i<(int)bl.size(); i++) {
            lambda[i] = rhs[i];
          }
          double* rPrev = nullptr;
          maxError = 0;
          for (int i=0; i<(int)rigidAtoms.size(); i++) {
            int b = i-1;
            double* ri = box.getAtomPosition(iFirstAtom+rigidAtoms[i]);
            double rij2 = 0;
            for (int k=0; k<3; k++) {
              ri[k] = rNew[i][k];
              if (i<(int)bl.size()) ri[k] += lambda[i] * rijOld[i][k];
              if (i>0) {
                ri[k] -= lambda[b]*rijOld[b][k];
                double rij = rPrev[k] - ri[k];
                rij2 += rij*rij;
              }
              else if (ring) ri[k] -= lambda[bl.size()-1]*rijOld[bl.size()-1][k];
            }
            if (i>0) {
              rhsOld[b] = rhs[b] = bl[b]*bl[b] - rij2;
              double e = fabs(rhs[b]) / (2*bl[b]*bl[b]);
              if (e>maxError) maxError = e;
            }
            rPrev = ri;
          }
          if (ring) {
            int b = bl.size()-1;
            double* ri = box.getAtomPosition(iFirstAtom+rigidAtoms[0]);
            double rij2 = 0;
            for (int k=0; k<3; k++) {
              double rij = rPrev[k] - ri[k];
              rij2 += rij*rij;
            }
            rhsOld[b] = rhs[b] = bl[b]*bl[b] - rij2;
            double e = fabs(rhs[b]) / (2*bl[b]*bl[b]);
            if (e>maxError) maxError = e;
          }
          iter++;
        }
        if (maxError > tol) {
          fprintf(stderr, "Could not converge to solution in MILC SHAKE (%d iterations, %e > %e)\n", maxIterations, maxError, tol);
          abort();
        }

        // update velocities
        for (int i=0; i<(int)rigidAtoms.size(); i++) {
          int b = i-1;
          double* vi = box.getAtomVelocity(iFirstAtom+rigidAtoms[i]);
          for (int k=0; k<3; k++) {
            if (i<(int)bl.size()) vi[k] += lambda[i] * rijOld[i][k] / dt;
            if (i>0) vi[k] -= lambda[b]*rijOld[b][k]/dt;
            else if (ring) vi[k] -= lambda[bl.size()-1]*rijOld[bl.size()-1][k]/dt;
          }
        }
      }
      iFirstAtom += na;
      oldOffset += na;
    }
  }
}

void ShakeListener::stepFinished() {
  int iFirstAtom = 0;
  double dt = integrator.getTimeStep();
  for (int iSpecies=0; iSpecies<speciesList.size(); iSpecies++) {
    vector<RigidConstraint*> &iConstraints = constraints[iSpecies];
    int nm = box.getNumMolecules(iSpecies);
    Species* s = speciesList.get(iSpecies);
    int na = s->getNumAtoms();
    if (iConstraints.size() == 0) {
      iFirstAtom += nm*na;
      continue;
    }
    for (int iMolecule=0; iMolecule<nm; iMolecule++) {
      for (vector<RigidConstraint*>::iterator it = iConstraints.begin(); it!=iConstraints.end(); it++) {
        const vector<double> & bl = (*it)->getBondLengths();
        const vector<int> & rigidAtoms = (*it)->getRigidAtoms();
        double *rPrev = box.getAtomPosition(iFirstAtom+rigidAtoms[0]);
        double *vPrev = box.getAtomVelocity(iFirstAtom+rigidAtoms[0]);
        for (int b=0; b<(int)bl.size(); b++) {
          int j = b+1;
          if (j==(int)rigidAtoms.size()) j=0;
          double* ri = box.getAtomPosition(iFirstAtom+rigidAtoms[j]);
          double* vi = box.getAtomVelocity(iFirstAtom+rigidAtoms[j]);
          rhs[b] = 0;
          if (b>0) dl[b-1] = 0;
          dd[b] = du[b] = 0;
          for (int k=0; k<3; k++) {
            rijNew[b][k] = rPrev[k] - ri[k];
            double vijk = vPrev[k] - vi[k];
            rhs[b] -= vijk*rijNew[b][k];

            if (b>0) dl[b-1] += -2*rijNew[b-1][k]*rijNew[b][k];
            dd[b] += 2*rijNew[b][k]*rijNew[b][k];
          }
          du[b] = dl[b];
          rPrev = ri;
          vPrev = vi;

          if (b>j) {
            // last ring bond.
            dl[b] = 0;
            for (int k=0; k<3; k++) {
              dl[b] += -2*rijNew[b][k]*rijNew[0][k];
            }
            du[b] = dl[b];
          }
        }
        tridiagSolve(bl.size(), dl, dd, du, rhs);
        const bool ring = rigidAtoms.size() == bl.size() && bl.size() > 2;

        // update velocities
        for (int i=0; i<(int)rigidAtoms.size(); i++) {
          int b = i-1;
          double* vi = box.getAtomVelocity(iFirstAtom+rigidAtoms[i]);
          for (int k=0; k<3; k++) {
            if (i<(int)bl.size()) vi[k] += rhs[i] * rijNew[i][k] / dt;
            if (i>0) vi[k] -= rhs[b]*rijNew[b][k]/dt;
            else if (ring) vi[k] -= rhs[bl.size()-1]*rijNew[bl.size()-1][k]/dt;
          }
        }
      }
      iFirstAtom += na;
    }
  }
}

// see LAPACK dgtsv
void ShakeListener::tridiagSolve(int n, double* dl, double* dd, double* du, double* rhs) {
  if (n==1) {
    rhs[0] /= dd[0];
    dl[0] = 0;
    return;
  }
  bool periodic = (du[0] == 0) && (dl[n-1] == 0);
  if (!periodic) {
    for (int i=0; i<n-2; i++) {
      if (fabs(dd[i]) >= fabs(dl[i])) {
        double f = dl[i] / dd[i];
        dd[i+1] -= f*du[i+1];
        rhs[i+1] -= f*rhs[i];
        dl[i] = 0;
      }
      else {
        // swap i and i+1
        double f = dd[i] / dl[i];
        dd[i] = dl[i];
        double t = dd[i+1];
        dd[i+1] = du[i+1] - f*t;
        dl[i] = du[i+2];  //dl no longer represents lower...
        du[i+2] *= -f*dl[i];
        du[i+1] = t;
        t = rhs[i];
        rhs[i] = rhs[i+1];
        rhs[i+1] = t - f*rhs[i+1];
      }
    }
    // last row
    int i = n-2;
    if (fabs(dd[i]) >= fabs(dl[i])) {
      double f = dl[i] / dd[i];
      dd[i+1] -= f*du[i+1];
      rhs[i+1] -= f*rhs[i];
    }
    else {
      // swap i and i+1
      double f = dd[i] / dl[i];
      dd[i] = dl[i];
      double t = dd[i+1];
      dd[i+1] = du[i+1] - f*t;
      du[i+1] = t;
      t = rhs[i];
      rhs[i] = rhs[i+1];
      rhs[i+1] = t - f*rhs[i+1];
    }
    // backsolve
    rhs[n-1] /= dd[n-1];
    rhs[n-2] = (rhs[n-2] - du[n-1]*rhs[n-1]) / dd[n-2];
    for (int i=n-3; i>=0; i--) {
      rhs[i] = (rhs[i] - du[i+1]*rhs[i+1] - dl[i]*rhs[i+2]) / dd[i];
    }
  }
  else {
  }
}
