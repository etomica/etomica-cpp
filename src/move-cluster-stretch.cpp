#include "move-virial.h"
#include "move.h"
#include "alloc2d.h"
/**
 * Monte Carlo move for Mayer sampling that explore bond length degrees of
 * freedom.  One bond length is varied in each molecule.  Atoms on each side of
 * the bond are translated with the overall molecule COM fixed.
 *
 */
MCMoveClusterStretch::MCMoveClusterStretch(Box& b, PotentialMasterIntra& p, Random& r, vector<vector <int>> bnd, SpeciesList& sl, double stepSize, Cluster& cluster) :
  MCMove(b, r, stepSize), speciesList(sl), cluster(cluster), bonding(bnd), iMolecule(0), potentialMasterIntra(p)
{
  maxStepSize = 0.2;
}

MCMoveClusterStretch::~MCMoveClusterStretch() {
}


// void MCMoveClusterStretch::setBox(Box& p)
// {
//
// }

bool MCMoveClusterStretch::doTrial() {
  if (true)
  {
    return false;
  }

  if (tunable && numTrials >= adjustInterval) {
    adjustStepSize();
  }

  int nm = mySpecies<0 ? box.getTotalNumMolecules() : box.getNumMolecules(mySpecies);
  if (nm==0) {
    iMolecule = -1;
    return false;
  }
  iMolecule = random.nextInt(nm);
  // iMolecule = 1;
  int iMoleculeInSpecies = 0;
  if (mySpecies>=0) {
    iMoleculeInSpecies = iMolecule;
    iMolecule = box.getGlobalMoleculeIndex(mySpecies, iMolecule);
  }
  uOld = potentialMasterIntra.oldMoleculeEnergyIntra(iMolecule);
  sum += uOld;
  counter++;

  // if (counter % 100000 == 0) {
  //   printf("%d counter = %d %f \n", idx, counter, sum/counter);
  // }

  wNew = wOld = 1;
  wOld = fabs(cluster.getValues()[0]);
  if (false) {
    // see if we can recompute the molecule energy and get the same result
    // this messes up future computation, so abort afterwards
    printf("got old, now recompute\n");
    double uTmp = potentialMasterIntra.computeOneMoleculeIntra(iMolecule);
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
  step = stepSize * (random.nextDouble() - 0.5);
  do{
    b = random.nextInt(bonding.size());
  }while (bonding[b].size() < 1);
  // b = 0;
  modified.push_back(iAtomFirst + b);
  int a = random.nextInt(bonding[b].size());
  a = bonding[b][a];
  // a = 1;
  modified.push_back(iAtomFirst + a);
  double dr[3];
  Vector::Ev1Mv2(box.getAtomPosition(iAtomFirst + b), box.getAtomPosition(iAtomFirst + a), dr);
  Vector::TE(step / sqrt(Vector::squared(dr)), dr);
  double shift[3] = {0, 0, 0};
  double m = 0;
  m += transform(dr, iAtomFirst + b, shift);
  m += transformBondedAtoms(dr, b, shift);
  Vector::TE(-1, dr);
  m += transform(dr, iAtomFirst + a, shift);
  m += transformBondedAtoms(dr, a, shift);
  if (fixedCOM) {
    Vector::TE(-1.0/m, shift);

  }
  else {
    Vector::Ev1Mv2(oldPositions[0], box.getAtomPosition(iAtomFirst), shift);
  }
  na = iAtomLast-iAtomFirst+1;
  for (int i=0; i<na; i++)
  {
    int iAtom = iAtomFirst + i;
    Vector::PE(shift, box.getAtomPosition(iAtom));
  }

  numTrials++;
  uNew = potentialMasterIntra.computeOneMoleculeIntra(iMolecule);
  wNew = fabs(cluster.getValues()[0]);
  if (idx == 1)
  {
    double arr[3];
    Vector::Ev1Mv2(box.getAtomPosition(iAtomFirst + 0), box.getAtomPosition(iAtomFirst +1), arr);
    // printf("%f %f \n", sqrt(Vector::squared(arr)), wNew );
    // if (sqrt(Vector::dot(box.getAtomPosition(iAtomFirst + 0), box.getAtomPosition(iAtomFirst +1))) > 2) exit(0);

  }



  return true;
}



double MCMoveClusterStretch::getChi(double temperature)
{
    double chi =  (wOld == 0 ? 1 : wNew / wOld) * exp(-(uNew - uOld) / temperature);
    chiSum += chi;
    // if (idx == 0)
    // {
    //   printf("%d %f %f %f %f %f \n", idx, temperature, uOld, uNew, wOld, wNew);
    //
    // }
    return chi;
}

void MCMoveClusterStretch::acceptNotify()
{
    // if (idx == 0) printf("accepted\n");
    potentialMasterIntra.setMoleculeEnergy(iMolecule, uNew);

    numAccepted++;

}

void MCMoveClusterStretch::rejectNotify()
{
    // if (idx == 0) printf("rejected\n");
    if (iMolecule < 0) return;
    int na = iAtomLast-iAtomFirst+1;
    for (int i=0; i<na; i++) {
        int iAtom = iAtomFirst + i;
        double *ri = box.getAtomPosition(iAtom);

        for (int j=0; j<3; j++) ri[j] = oldPositions[i][j];
    }


}

double MCMoveClusterStretch::transform(double* dr, int index, double* shift)
{
    Vector::PE(dr, box.getAtomPosition(index));
    int iType = box.getAtomType(index);
    double m = speciesList.getAtomInfo().getMass(iType);
    Vector::PEa1Tv1(m, dr, shift);
    return m;
}

double MCMoveClusterStretch::transformBondedAtoms(double* dr, const int index, double* shift)
{
  double m = 0;
  for(int k : bonding[index]){
    bool alreadyMoved = false;
    for (int l : modified) {
      if (iAtomFirst + k == l) {
        alreadyMoved = true;
        break;
      }
    }
    if (!alreadyMoved) {
      m += transform(dr, iAtomFirst + k, shift);
      m += transformBondedAtoms(dr, k, shift);
    }
  }
  return m;
}
double MCMoveClusterStretch::energyChange() {
  //bogus
  return 0;
}
