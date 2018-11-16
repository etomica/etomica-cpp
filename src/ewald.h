/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include <vector>
#include <complex>

using namespace std;

class SpeciesList;
class Box;
class PotentialCallback;

class EwaldBase {
  protected:
    const SpeciesList& speciesList;
    Box& box;

    bool rigidMolecules;
    double* charges;

  public:
    EwaldBase(const SpeciesList &speciesList, Box& box);
    virtual ~EwaldBase();
    virtual void init() {}
    void setRigidMolecules(bool rigidMolecules);
    virtual void computeAllFourier(const bool doForces, const bool doPhi, const bool doDFDV, double &uTot, double &virialTot, double** force, vector<PotentialCallback*>* pairCallbacks) = 0;
    virtual void computeFourierIntramolecular(int iMolecule, const bool doForces, const bool doPhi, double &uTot, double &virialTot, double** force) = 0;
    virtual double oneMoleculeFourierEnergy(int iMolecule, bool oldEnergy) = 0;
    virtual void processAtomU(int coeff) = 0;
    virtual void resetAtomDU() = 0;
    virtual double uTotalFromAtoms() = 0;

    virtual void setCharge(int iType, double charge);
};

class EwaldFourier : public EwaldBase {
  protected:
    double kBasis[3];
    double kCut, alpha;
    complex<double> *sFacAtom;
    vector<complex<double> > sFac;
    vector<complex<double> > eik[3];
    vector<complex<double> > dsFacMolecule;
    vector<double> fExp;
    double phi[3][3];

  public:
    EwaldFourier(const SpeciesList &speciesList, Box& box);
    virtual ~EwaldFourier();
    void computeAllFourier(const bool doForces, const bool doPhi, const bool doDFDV, double &uTot, double &virialTot, double** force, vector<PotentialCallback*>* pairCallbacks);
    void computeFourierIntramolecular(int iMolecule, const bool doForces, const bool doPhi, double &uTot, double &virialTot, double** force);
    double oneMoleculeFourierEnergy(int iMolecule, bool oldEnergy);

    void setParameters(double kCut, double alpha);
    double uTotalFromAtoms();
    void processAtomU(int coeff);
    void resetAtomDU();

};
