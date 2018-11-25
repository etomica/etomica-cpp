/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "neighbor-update-listener.h"
#include "potential-master.h"

NeighborUpdateListener::NeighborUpdateListener(PotentialMasterList& pm, int i) : potentialMaster(pm), interval(i), intervalCountdown(i) {
}

void NeighborUpdateListener::stepStarted() {
  intervalCountdown--;
  if (intervalCountdown) return;
  potentialMaster.checkUpdateNbrs();
  intervalCountdown = interval;
}

void NeighborUpdateListener::setInterval(int newInterval) {
  interval = newInterval;
  intervalCountdown = newInterval;
}
