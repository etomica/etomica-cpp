/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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
  const double* bs = box.getBoxSize();
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
