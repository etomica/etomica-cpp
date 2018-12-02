/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "integrator.h"
#include "move.h"

IntegratorMC::IntegratorMC(PotentialMaster& p, Random& r) : Integrator(p), random(r), pMoveSum(0), lastMove(NULL) {
}

IntegratorMC::~IntegratorMC(){}

void IntegratorMC::addMove(MCMove* move, double prob) {
  moves.push_back(move);
  moveProbabilities.push_back(prob);
  pMoveSum += prob;
}

void IntegratorMC::removeMove(MCMove* move) {
  pMoveSum = 0;
  for (int i=0; i<(int)moves.size(); i++) {
    if (moves[i] == move) {
      moves.erase(moves.begin()+i);
      moveProbabilities.erase(moveProbabilities.begin()+i);
      for (int j=i; j<(int)moves.size(); j++) pMoveSum += moveProbabilities[j];
      return;
    }
    pMoveSum += moveProbabilities[i];
  }
  fprintf(stderr, "Could not find listener");
  abort();
}

void IntegratorMC::setTuning(bool doTune) {
  for (vector<MCMove*>::iterator it = moves.begin(); it!=moves.end(); it++) {
    (*it)->tunable = doTune;
  }
}

void IntegratorMC::doStep() {
  stepCount++;
  for (vector<IntegratorListener*>::iterator it = listenersStepStarted.begin(); it!=listenersStepStarted.end(); it++) {
    (*it)->stepStarted();
  }
  MCMove* m = nullptr;
  int nm = moves.size();
  if (nm>1) {
    double r = random.nextDouble32()*pMoveSum;
    double s = 0;
    for (int i=0; i<nm; i++) {
      s += moveProbabilities[i];
      if (s >= r) {
        m = moves[i];
        break;
      }
    }
  }
  else if (nm==0) {
    fprintf(stderr, "I need MC moves\n");
    abort();
  }
  else {
    m = moves[0];
  }
  lastMove = m;
  bool success = m->doTrial();
  double chi = success ? m->getChi(temperature) : 0;
  if (chi==0 || (chi<1 && chi<random.nextDouble())) {
    //printf("%ld chi %e rej\n", stepCount, chi);
    m->rejectNotify();
    for (vector<IntegratorListener*>::iterator it = listenersMoveRejected.begin(); it!=listenersMoveRejected.end(); it++) {
      (*it)->moveRejected(*m, chi);
    }
  }
  else {
    m->acceptNotify();
    double du = m->energyChange();
    energy += du;
    //printf("%ld chi %e acc %f\n", stepCount, chi, du);
    for (vector<IntegratorListener*>::iterator it = listenersMoveAccepted.begin(); it!=listenersMoveAccepted.end(); it++) {
      (*it)->moveAccepted(*m, chi);
    }
  }
#ifdef DEBUG
  if (fabs(energy-potentialMaster.uTotalFromAtoms()) > 1e-4) {
    double x = potentialMaster.uTotalFromAtoms();
    double y = energy;
    reset();
    printf("uAtoms trouble! step %ld:\n", stepCount);
    printf("incrementally updated:  % e\n", y);
    printf("uTotalFromAtoms:        % e\n", x);
    printf("actual energy:          % e\n", energy);
    printf("incremental off by:     % e\n", y-energy);
    printf("uTotalFromAtoms off by: % e\n", x-energy);
    abort();
  }
  if (stepCount%100==0) {
    double oldEnergy = energy;
    reset();
    if (fabs(oldEnergy-energy) > 1e-4) printf("%ld: %e %e %e\n", stepCount, oldEnergy, energy, oldEnergy-energy);
    if (fabs(oldEnergy-energy) > 1e-4) abort();
  }
#endif
  for (vector<IntegratorListener*>::iterator it = listenersStepFinished.begin(); it!=listenersStepFinished.end(); it++) {
    (*it)->stepFinished();
  }
}

MCMove* IntegratorMC::getLastMove() {
  return lastMove;
}

void IntegratorMC::addListener(IntegratorListener* listener) {
  Integrator::addListener(listener);
  if (listener->callAccept) listenersMoveAccepted.push_back(listener);
  if (listener->callReject) listenersMoveRejected.push_back(listener);
}

void IntegratorMC::removeListener(IntegratorListener* listener) {
  Integrator::removeListener(listener);
  if (listener->callAccept) {
    std::vector<IntegratorListener*>::iterator it = find (listenersMoveAccepted.begin(), listenersMoveAccepted.end(), listener);
    if (it != listenersMoveAccepted.end()) listenersMoveAccepted.erase(it);
  }
  if (listener->callReject) {
    std::vector<IntegratorListener*>::iterator it = find (listenersMoveRejected.begin(), listenersMoveRejected.end(), listener);
    if (it != listenersMoveRejected.end()) listenersMoveRejected.erase(it);
  }
}
