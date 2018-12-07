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
    int* numAtomsByType;

    bool rigidMolecules;
    double* charges;
    double** B6;
    double** b6;

  public:
    EwaldBase(const SpeciesList &speciesList, Box& box);
    virtual ~EwaldBase();
    virtual void init() {}
    void setRigidMolecules(bool rigidMolecules);
    void setNumAtomsByType(int* numAtomsByType);
    virtual void computeAllFourier(const bool doForces, const bool doPhi, const bool doDFDV, const bool doVirialTensor, double &uTot, double &virialTot, double** force, double* virialTensor, vector<PotentialCallback*>* pairCallbacks) = 0;
    virtual void computeFourierIntramolecular(int iMolecule, const bool doForces, const bool doPhi, double &uTot, double &virialTot, double** force) = 0;
    virtual double oneMoleculeFourierEnergy(int iMolecule, bool oldEnergy) = 0;
    virtual double computeFourierAtom(int iAtom, bool oldEnergy) = 0;
    virtual void processAtomU(int coeff) = 0;
    virtual void resetAtomDU() = 0;
    virtual double uTotalFromAtoms() = 0;

    virtual void setCharge(int iType, double charge);
    virtual void setR6Coeff(int iType, double sigma, double epsilon);
};

class EwaldFourier : public EwaldBase {
  protected:
    double kBasis[3];
    double kCut, alpha, eta;
    vector<complex<double> > sFacAtom;
    vector<complex<double> > sFac;
    vector<complex<double> > sFacB[7];
    vector<complex<double> > eik[3];
    vector<complex<double> > dsFacMolecule;
    vector<complex<double> > dsFacBMolecule[7];
    vector<double> fExp, f6Exp;
    double phi[3][3];

  public:
    EwaldFourier(const SpeciesList &speciesList, Box& box);
    virtual ~EwaldFourier();
    void computeAllFourier(const bool doForces, const bool doPhi, const bool doDFDV, const bool doVirialTensor, double &uTot, double &virialTot, double** force, double* virialTensor, vector<PotentialCallback*>* pairCallbacks);
    void computeFourierIntramolecular(int iMolecule, const bool doForces, const bool doPhi, double &uTot, double &virialTot, double** force);
    double oneMoleculeFourierEnergy(int iMolecule, bool oldEnergy);
    double computeFourierAtom(int iAtom, bool oldEnergy);

    void setCutoff(double kCut);
    void setChargeAlpha(double alpha);
    void setR6eta(double eta);
    double uTotalFromAtoms();
    void processAtomU(int coeff);
    void resetAtomDU();
    void getOptimalAlpha(double s, double& alpha, double& rc, double& kc);

};
