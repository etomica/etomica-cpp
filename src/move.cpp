#include "move.h"

MCMove::MCMove(Box& b, PotentialMaster& p, Random& r, double ss) : box(b), potentialMaster(p), random(r), stepSize(ss) {
  init();
}

MCMove::~MCMove() {}

void MCMove::init() {
  numTrials = numAccepted = chiSum = 0;
  lastAdjust = 0;
  adjustInterval = 100;
  adjustStep = 1.05;
  minAdjustStep = 1;
  verboseAdjust = false;
  tunable = true;
  maxStepSize = 0;
  double* bs = box.getBoxSize();
  for (int i=0; i<3; i++) if (maxStepSize<bs[i]) maxStepSize = bs[i];
  maxStepSize /= 2;
}

void MCMove::setStepSize(double ss) {
  stepSize = ss;
}

double MCMove::getStepSize() {
  return stepSize;
}

double MCMove::getAcceptance() {
  if (numTrials==0) return 0;
  return chiSum/numTrials;
}

void MCMove::adjustStepSize() {
  double avg = chiSum/numTrials;
  if (avg > 0.5) {
    if (stepSize < maxStepSize) {
      if (lastAdjust < 0) {
        // back and forth
        adjustInterval *= 2;
        adjustStep = sqrt(adjustStep);
      }
      else if (lastAdjust == 5) {
        // sixth consecutive increase; increase adjustment step
        adjustStep *= adjustStep;
        if (adjustStep > 2) {
          adjustStep = 2;
        }
        lastAdjust = 3;
      }
      stepSize *= adjustStep;
      stepSize = std::min(stepSize, maxStepSize);
      if (verboseAdjust) {
        printf("move increasing step size: %f (<chi> = %f)\n", stepSize, avg);
      }
      if (lastAdjust < 1) lastAdjust = 1;
      else lastAdjust++;
    }
    else if (verboseAdjust) {
      printf("move step size: %f (<chi> = %f\n)", stepSize, avg);
    }
  }
  else {
    if (lastAdjust > 0) {
      // back and forth
      adjustInterval *= 2;
      adjustStep = sqrt(adjustStep);
    }
    else if (lastAdjust == -5) {
      // sixth consecutive increase; increase adjustment step
      adjustStep *= adjustStep;
      if (adjustStep > 2) {
        adjustStep = 2;
      }
      lastAdjust = -3;
    }
    stepSize /= adjustStep;
    if (verboseAdjust) {
      printf("move decreasing step size: %f (<chi> = %f)\n", stepSize, avg);
    }
    if (lastAdjust > -1) lastAdjust = -1;
    else lastAdjust--;
  }
  numTrials = numAccepted = 0;
  chiSum = 0;
}


MCMoveDisplacement::MCMoveDisplacement(Box& b, PotentialMaster& p, Random& r, double ss) : MCMove(b,p,r,ss) {
}

MCMoveDisplacement::~MCMoveDisplacement() {}

bool MCMoveDisplacement::doTrial() {
  if (tunable && numTrials >= adjustInterval) {
    adjustStepSize();
  }
  iAtom = random.nextInt(box.getNumAtoms());
  uOld = potentialMaster.oldEnergy(iAtom);
  double* r = box.getAtomPosition(iAtom);
  std::copy(r, r+3, rOld);
  r[0] += 2*stepSize*(random.nextDouble32()-0.5);
  r[1] += 2*stepSize*(random.nextDouble32()-0.5);
  r[2] += 2*stepSize*(random.nextDouble32()-0.5);
  box.nearestImage(r);
  numTrials++;
  return true;
}

double MCMoveDisplacement::getChi(double T) {
  uNew = 0;
  potentialMaster.computeOne(iAtom, box.getAtomPosition(iAtom), uNew, true);
  double chi = uNew<uOld ? 1 : exp(-(uNew-uOld)/T);
  chiSum += chi;
  return chi;
}

void MCMoveDisplacement::acceptNotify() {
  //printf("accepted\n");
  potentialMaster.processAtomU(1);
  double uTmp = 0;
  potentialMaster.computeOne(iAtom, rOld, uTmp, false);
  potentialMaster.processAtomU(-1);
  potentialMaster.updateAtom(iAtom);
  numAccepted++;
}

void MCMoveDisplacement::rejectNotify() {
  //printf("rejected\n");
  double* r = box.getAtomPosition(iAtom);
  std::copy(rOld, rOld+3, r);
  uNew = uOld;
  potentialMaster.resetAtomDU();
}

double MCMoveDisplacement::energyChange() {
  return uNew-uOld;
}


MCMoveInsertDelete::MCMoveInsertDelete(Box& b, PotentialMaster& p, Random& r, double m) : MCMove(b,p,r,0), mu(m) {
  tunable = false;
}

MCMoveInsertDelete::~MCMoveInsertDelete() {}

bool MCMoveInsertDelete::doTrial() {
  doInsert = random.nextInt(2) == 0;
  int n = box.getNumAtoms();
  if (doInsert) {
    uOld = 0;
    double *bs = box.getBoxSize();
    for (int j=0; j<3; j++) {
      rNew[j] = bs[j]*(random.nextDouble32()-0.5);
    }
    potentialMaster.computeOne(n, rNew, uNew, true);
  }
  else {
    // delete
    if (n>0) {
      iAtom = random.nextInt(n);
      uOld = potentialMaster.oldEnergy(iAtom);
      uNew = 0;
    }
  }
  numTrials++;
  return true;
}

double MCMoveInsertDelete::getChi(double T) {
  int n = box.getNumAtoms();
  if (n==0 && !doInsert) return 0;

  double* bs = box.getBoxSize();
  double vol = bs[0]*bs[1]*bs[2];

  double x = -(uNew-uOld);
  double a;
  if (doInsert) {
    x += mu;
    a = vol/(n+1);
  }
  else {
    x -= mu;
    a = (n-1)/vol;
  }
  double chi = a*exp(x/T);
  //printf("%d %d %f %f %f %f %f %f\n", doInsert, n, a, uOld, uNew, mu, x, chi);
  if (chi>1) chi = 1;
  chiSum += chi;
  return chi;
}

void MCMoveInsertDelete::acceptNotify() {
  int n = box.getNumAtoms();
  if (doInsert) {
    //printf("accept insert %d\n", n+1);
    box.setNumAtoms(n+1);
    double *ri = box.getAtomPosition(n);
    std::copy(rNew, rNew+3, ri);
    potentialMaster.newAtom();
    potentialMaster.processAtomU(1);
  }
  else {
    //printf("accept delete %d\n", n-1);
    double *ri = box.getAtomPosition(iAtom);
    double uTmp = 0;
    potentialMaster.computeOne(iAtom, ri, uTmp, false);
    potentialMaster.processAtomU(-1);
    // this removes iAtom (updates cell lists) and then
    // moves the last atom into its position
    potentialMaster.removeAtom(iAtom);
    double *rn = box.getAtomPosition(n-1);
    std::copy(rn, rn+3, ri);
    box.setNumAtoms(n-1);
  }

  numAccepted++;
}

void MCMoveInsertDelete::rejectNotify() {
  //if (doInsert) printf("reject insert\n");
  //else printf("reject delete\n");
  uNew = uOld;
  potentialMaster.resetAtomDU();
}

double MCMoveInsertDelete::energyChange() {
  return uNew-uOld;
}
