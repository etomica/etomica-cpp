#include "integrator.h"

IntegratorMC::IntegratorMC(PotentialMaster& p, Random& r) : potentialMaster(p), sfmt(r.sfmt), temperature(1), stepCount(0), pMoveSum(0), lastMove(NULL) {
  callFinished = true;
  selfPotentialCallbackVec.push_back(this);
}

IntegratorMC::IntegratorMC(PotentialMaster& p, sfmt_t& s) : potentialMaster(p), sfmt(s), temperature(1), stepCount(0), pMoveSum(0), lastMove(NULL) {
  callFinished = true;
  selfPotentialCallbackVec.push_back(this);
}

IntegratorMC::~IntegratorMC(){}

void IntegratorMC::allComputeFinished(double uTotNew, double virialTotNew, double **f) {
  energy = uTotNew;
}

void IntegratorMC::addMove(MCMove* move, double prob) {
  moves.push_back(move);
  moveProbabilities.push_back(prob);
  pMoveSum += prob;
}

void IntegratorMC::setTuning(bool doTune) {
  for (vector<MCMove*>::iterator it = moves.begin(); it!=moves.end(); it++) {
    (*it)->tunable = doTune;
  }
}

void IntegratorMC::doStep() {
  stepCount++;
  MCMove* m = nullptr;
  int nm = moves.size();
  if (nm>1) {
    double r = sfmt_genrand_real1(&sfmt)*pMoveSum;
    double s = 0;
    for (int i=0; i<nm; i++) {
      s += moveProbabilities[i];
      if (s >= r) {
        m = moves[i];
        break;
      }
    }
  }
  else {
    m = moves[0];
  }
  lastMove = m;
  m->doTrial();
  double chi = m->getChi(temperature);
  if (chi<1 && chi<sfmt_genrand_res53(&sfmt)) {
    //printf("chi %e 0\n", chi);
    m->rejectNotify();
  }
  else {
    //printf("chi %e 1\n", chi);
    m->acceptNotify();
    energy += m->energyChange();
  }
  if (false && stepCount%100==0) {
    double oldEnergy = energy;
    reset();
    if (fabs(oldEnergy-energy) > 1e-7) printf("%d: %e %e %e\n", stepCount, oldEnergy, energy, oldEnergy-energy);
  }
  for (vector<IntegratorListener*>::iterator it = listeners.begin(); it!=listeners.end(); it++) {
    (*it)->stepFinished();
  }
}

void IntegratorMC::doSteps(int steps) {
  for (int i=0; i<steps; i++) {
    doStep();
  }
}

void IntegratorMC::setTemperature(double T) {
  temperature = T;
}

double IntegratorMC::getTemperature() {
  return temperature;
}

void IntegratorMC::reset() {
  potentialMaster.computeAll(selfPotentialCallbackVec);
}

double IntegratorMC::getPotentialEnergy() {
  return energy;
}

void IntegratorMC::addListener(IntegratorListener* listener) {
  listeners.push_back(listener);
}

MCMove* IntegratorMC::getLastMove() {
  return lastMove;
}

IntegratorListener::IntegratorListener() {}
