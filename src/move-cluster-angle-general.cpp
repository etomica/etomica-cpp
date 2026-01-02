#include "move-virial.h"
#include "move.h"
#include "alloc2d.h"
/**
 * Monte Carlo move for Mayer sampling that explore bond angle degrees of
 * freedom.  One bond angle is varied in each molecule.  Atoms on each side of
 * the bond are rotated around the middle atom of the bond angle with the
 * overall molecule COM fixed.
 *
 * Originally developed by Arpit for alkanes.
 */

MCMoveClusterAngleGeneral::MCMoveClusterAngleGeneral(SpeciesList& sl, Box& b, PotentialMaster& p, bool oneSide, Random& r, double stepSize, Cluster& cluster) :
  MCMove(b, p, r, stepSize),
  speciesList(sl), atomInfo(atomInfo), oneSide(oneSide), cluster(cluster), iMolecule(0)
{
  maxStepSize = M_PI / 2;
}

// void MCMoveClusterAngleGeneral::setBox(Box& p)
// {
//
// }

bool MCMoveClusterAngleGeneral::doTrial() {
  if (tunable && numTrials >= adjustInterval) {
    adjustStepSize();
  }
  int nm = iSpecies<0 ? box.getTotalNumMolecules() : box.getNumMolecules(iSpecies);
  if (nm==0) {
    iMolecule = -1;
    return false;
  }
  iMolecule = random.nextInt(nm);
  int iMoleculeInSpecies = 0;
  if (iSpecies>=0) {
    iMoleculeInSpecies = iMolecule;
    iMolecule = box.getGlobalMoleculeIndex(iSpecies, iMolecule);
  }
  uOld = potentialMaster.oldMoleculeEnergy(iMolecule);
  wNew = wOld = 1;
  wOld = cluster.getValues()[0];
  if (false) {
    // see if we can recompute the molecule energy and get the same result
    // this messes up future computation, so abort afterwards
    printf("got old, now recompute\n");
    double uTmp = 0;
    potentialMaster.computeOneMolecule(iMolecule, uTmp);
    printf("%d uOlds %f %f\n", iMolecule, uTmp, uOld);
    abort();
  }
  int na = numOldPositions;
  mySpecies = iSpecies;
  if (iSpecies<0) {
    box.getMoleculeInfo(iMolecule, mySpecies, iMoleculeInSpecies, iAtomFirst, iAtomLast);
    na = iAtomLast-iAtomFirst+1;
    if (na==1) {
      iMolecule = -1;
      return false;
    }
    if (na>numOldPositions) {
      oldPositions = (double**)realloc2D((void**)oldPositions, na, 3, sizeof(double));
      numOldPositions = na;
    }
  }
  else {
    iAtomFirst = box.getFirstAtom(mySpecies, iMoleculeInSpecies);
    iAtomLast = iAtomFirst + na - 1;
  }
  modifiedIndex = 0;
  int d = random.nextInt(triplets.size());
  b = triplets[d][1];
  modified[modifiedIndex]=b;
  ++modifiedIndex;
  double axis[3];
  double temp[3];
  if (oneSide)
  {
    random.onSphere(axis);
  }
  else
  {
    Vector::Ev1Mv2(box.getAtomPosition(triplets[d][0]), box.getAtomPosition(b), axis);
    Vector::Ev1Mv2(box.getAtomPosition(triplets[d][2]), box.getAtomPosition(b), temp);
    Vector::cross(axis, temp, axis);
    Vector::normalize(axis);

  }
  double theta = stepSize*2*(random.nextDouble32()-0.5); //dt
  mat.setAxisAngle(axis, theta);
  Species* species = speciesList.get(mySpecies);
  double* center = species->getMoleculeCOM(box, iAtomFirst, iAtomLast);
  double shift[3];
  transform(mat, triplets[d][0], shift);
  transformBondedAtoms(mat, triplets[d][0], shift);
  if (!oneSide) {
    mat.setAxisAngle(axis, -dt);

    transform(mat, triplets[d][2], shift);

    transformBondedAtoms(mat, triplets[d][2], shift);
  }
  if (iMolecule==0 || (constraintMap != nullptr && constraintMap[iMolecule] == 0)) {
    double mt = 0;
    int na = iAtomLast-iAtomFirst+1;
    for (int i=0; i<na; i++)
    {
      int iAtom = iAtomFirst + i;
      int iType = box.getAtomType(iAtom);
      double m = atomInfo.getMass(iType);
      mt += m;
    }
    Vector::TE(-1.0/mt, shift);
    for (int i=0; i<na; i++)
    {
      int iAtom = iAtomFirst + i;
      Vector::PE(shift, box.getAtomPosition(iAtom));
    }

  }




  numTrials++;
  uNew = potentialMaster.oldMoleculeEnergy(iMolecule);
  wNew = cluster.getValues()[0];

  return true;
}

double MCMoveMoleculeRotate::getChi(double T) {
  uNew = 0;
  potentialMaster.computeOneMolecule(iMolecule, uNew);
  double chi = uNew<uOld ? 1 : exp(-(uNew-uOld)/T);
  //printf("%f => %f  ==>  %f\n", uOld, uNew, uNew-uOld);
  chiSum += chi;
  return chi;
}


double MCMoveClusterAngleGeneral::getChi(double temperature)
{
    return (wOld == 0 ? 1 : wNew / wOld) * exp(-(uNew - uOld) / temperature);
}

void MCMoveClusterAngleGeneral::acceptNotify()
{
    //printf("accepted\n");
    potentialMaster.processAtomU(1);
    int na = iAtomLast-iAtomFirst+1;
    for (int i=0; i<na; i++) {
        int iAtom = iAtomFirst + i;
        double *ri = box.getAtomPosition(iAtom);
        for (int j=0; j<3; j++) ri[j] = oldPositions[i][j];
        potentialMaster.updateAtom(iAtom);
    }
    double uTmp = 0;
    potentialMaster.computeOneMolecule(iMolecule, uTmp);
    // this call is designed to set up the next call.  uTmp won't necessarily be correct
    potentialMaster.processAtomU(-1);
    Species* species = speciesList.get(mySpecies);
    double* center = species->getMoleculeCOM(box, iAtomFirst, iAtomLast);
    for (int i=0; i<na; i++) {
        int iAtom = iAtomFirst + i;
        double *ri = box.getAtomPosition(iAtom);
        mat.transformAbout(ri, center, box);
        potentialMaster.updateAtom(iAtom);
    }
    numAccepted++;

}

void MCMoveClusterAngleGeneral::rejectNotify()
{
    //printf("rejected\n");
    if (iMolecule < 0) return;
    int na = iAtomLast-iAtomFirst+1;
    for (int i=0; i<na; i++) {
        int iAtom = iAtomFirst + i;
        double *ri = box.getAtomPosition(iAtom);
        for (int j=0; j<3; j++) ri[j] = oldPositions[i][j];
        potentialMaster.updateAtom(iAtom);
    }
    uNew = uOld;
    potentialMaster.resetAtomDU();

}

void MCMoveClusterAngleGeneral::transform(RotationMatrix rm, const int index, double* shift)
{
    double r[3];
    double* p = box.getAtomPosition(index);
    int iType = box.getAtomType(index);
    double m = atomInfo.getMass(iType);
    Vector::PEa1Tv1(-m, p, shift);
    Vector::Ev1Mv2(p, box.getAtomPosition(b), r);
    rm.transform(r);
    Vector::PE(box.getAtomPosition(b), r);
    Vector::E(r, p);
    Vector::PEa1Tv1(m, p, shift);
    modified[modifiedIndex] = index;
    modifiedIndex++;

}

void MCMoveClusterAngleGeneral::transformBondedAtoms(RotationMatrix rm, const int index, double* shift)
{
  for(int k : bonding[index]){
    bool rotated = false;
    for (int l = 0; l < modifiedIndex; l++) {
      if (k == modified[l]) {
        rotated = true;
        break;
      }
    }
    if (!rotated) {
      transform(rm, k, shift);
      transformBondedAtoms(rm, k, shift);
    }
  }

}
