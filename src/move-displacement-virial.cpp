#include "move-virial.h"

MCMoveDisplacementVirial::MCMoveDisplacementVirial(Box& b, PotentialMaster& p, Random& r, double ss, Cluster& c) : MCMove(b,p,r,ss), cluster(c) {
}

MCMoveDisplacementVirial::~MCMoveDisplacementVirial() {}

bool MCMoveDisplacementVirial::doTrial() {
  if (tunable && numTrials >= adjustInterval) {
    adjustStepSize();
  }
  int na = box.getNumAtoms();
  if (na<=1) {
    iAtom = -1;
    return false;
  }
  iAtom = random.nextInt(na-1) + 1;
  wOld = fabs(cluster.getValues()[0]);
  double* r = box.getAtomPosition(iAtom);
  std::copy(r, r+3, rOld);
  r[0] += 2*stepSize*(random.nextDouble32()-0.5);
  r[1] += 2*stepSize*(random.nextDouble32()-0.5);
  r[2] += 2*stepSize*(random.nextDouble32()-0.5);
  cluster.trialNotify();
  numTrials++;
  return true;
}

double MCMoveDisplacementVirial::getChi(double T) {
  wNew = fabs(cluster.getValues()[0]);
  double chi = wNew>wOld ? 1 : wNew/wOld;
  chiSum += chi;
  return chi;
}

void MCMoveDisplacementVirial::acceptNotify() {
  //printf("accepted\n");
  numAccepted++;
}

void MCMoveDisplacementVirial::rejectNotify() {
  //printf("rejected\n");
  double* r = box.getAtomPosition(iAtom);
  std::copy(rOld, rOld+3, r);
}

double MCMoveDisplacementVirial::energyChange() {
  //bogus
  return 0;
}
