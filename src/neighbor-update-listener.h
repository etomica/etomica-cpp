/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

#include "integrator.h"

class PotentialMasterList;

class NeighborUpdateListener : public IntegratorListener {
  private:
    PotentialMasterList& potentialMaster;
    int interval, intervalCountdown;

  public:
    NeighborUpdateListener(PotentialMasterList& potentialMaster, int interval);
    ~NeighborUpdateListener() {}
    void stepStarted();
    void setInterval(int newInterval);
};
