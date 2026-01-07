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

MCMoveClusterAngleGeneral::MCMoveClusterAngleGeneral(Box& b, PotentialMaster& p, bool oneSide, Random& r, vector<vector <int>> bnd, vector<int*> t, SpeciesList& sl, double stepSize, Cluster& cluster) :
  MCMove(b, p, r, stepSize), speciesList(sl),
   oneSide(oneSide), cluster(cluster), triplets(t), bonding(bnd), iMolecule(0)
{
  maxStepSize = M_PI / 2;
}

MCMoveClusterAngleGeneral::~MCMoveClusterAngleGeneral() {
}


// void MCMoveClusterAngleGeneral::setBox(Box& p)
// {
//
// }

bool MCMoveClusterAngleGeneral::doTrial() {
  if (tunable && numTrials >= adjustInterval) {
    adjustStepSize();
  }
  int nm = mySpecies<0 ? box.getTotalNumMolecules() : box.getNumMolecules(mySpecies);
  if (nm==0) {
    iMolecule = -1;
    return false;
  }
  iMolecule = random.nextInt(nm);
  int iMoleculeInSpecies = 0;
  if (mySpecies>=0) {
    iMoleculeInSpecies = iMolecule;
    iMolecule = box.getGlobalMoleculeIndex(mySpecies, iMolecule);
  }
  uOld = potentialMaster.oldMoleculeEnergy(iMolecule);
  wNew = wOld = 1;
  wOld = fabs(cluster.getValues()[0]);
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
  if (mySpecies<0) {
    int iSpecies;
    box.getMoleculeInfo(iMolecule, iSpecies, iMoleculeInSpecies, iAtomFirst, iAtomLast);
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
    fprintf(stderr, "I don't know how to get atoms from a specific species \n");
    abort();
    iAtomFirst = box.getFirstAtom(mySpecies, iMoleculeInSpecies);
    iAtomLast = iAtomFirst + na - 1;
  }
  for (int i=0; i<na; i++) {
    int iAtom = iAtomFirst + i;
    double *ri = box.getAtomPosition(iAtom);

    for (int j=0; j<3; j++) oldPositions[i][j] = ri[j];
  }

  modified.clear();
  int d = random.nextInt(triplets.size());
  b = triplets[d][1];
  modified.push_back(iAtomFirst + b);
  double temp1[3];
  double axis[3];
  double temp[3];
  if (oneSide)
  {
    random.onSphere(axis);
  }
  else
  {
    Vector::Ev1Mv2(box.getAtomPosition(triplets[d][0]), box.getAtomPosition(b), temp1);
    Vector::Ev1Mv2(box.getAtomPosition(triplets[d][2]), box.getAtomPosition(b), temp);
    Vector::cross(temp1, temp, axis);
    Vector::normalize(axis);

  }
  double theta = stepSize*2*(random.nextDouble32()-0.5); //dt
  mat.setAxisAngle(axis, theta);
  double shift[3] = {0, 0, 0};
  transform(iAtomFirst + triplets[d][0], shift);
  transformBondedAtoms(triplets[d][0], shift);
  if (!oneSide) {
    mat.setAxisAngle(axis, -theta);

    transform(iAtomFirst + triplets[d][2], shift);

    transformBondedAtoms(triplets[d][2], shift);
  }

  if (iMolecule==0 || (constraintMap != nullptr && constraintMap[iMolecule] == 0)) {
    double mt = 0;
    int na = iAtomLast-iAtomFirst+1;
    for (int i=0; i<na; i++)
    {
      int iAtom = iAtomFirst + i;
      int iType = box.getAtomType(iAtom);
      double m = speciesList.getAtomInfo().getMass(iType);
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
  potentialMaster.computeOneMolecule(iMolecule, uNew);
  wNew = fabs(cluster.getValues()[0]);
  return true;
}



double MCMoveClusterAngleGeneral::getChi(double temperature)
{
    double chi =  (wOld == 0 ? 1 : wNew / wOld) * exp(-(uNew - uOld) / temperature);
    chiSum += chi;
    return chi;
}

void MCMoveClusterAngleGeneral::acceptNotify()
{
    //printf("accepted\n");
    numAccepted++;

}

void MCMoveClusterAngleGeneral::rejectNotify()
{
    // printf("rejected\n");
    if (iMolecule < 0) return;
    int na = iAtomLast-iAtomFirst+1;
    for (int i=0; i<na; i++) {
        int iAtom = iAtomFirst + i;
        double *ri = box.getAtomPosition(iAtom);

        for (int j=0; j<3; j++) ri[j] = oldPositions[i][j];
    }
}

void MCMoveClusterAngleGeneral::transform(const int index, double* shift)
{
    double* p = box.getAtomPosition(index);
    int iType = box.getAtomType(index);
    double m = speciesList.getAtomInfo().getMass(iType);
    Vector::PEa1Tv1(-m, p, shift);
    mat.transformAbout(p, box.getAtomPosition(iAtomFirst + b), box);
    Vector::PEa1Tv1(m, p, shift);
    modified.push_back(index);
}

void MCMoveClusterAngleGeneral::transformBondedAtoms(const int index, double* shift)
{
  for(int k : bonding[index]){
    bool rotated = false;
    for (int l : modified) {
      if (iAtomFirst + k == l) {
        rotated = true;
        break;
      }
    }
    if (!rotated) {
      transform(iAtomFirst + k, shift);
      transformBondedAtoms(k, shift);
    }
  }

}
double MCMoveClusterAngleGeneral::energyChange() {
  //bogus
  return 0;
}
