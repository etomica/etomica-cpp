/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#pragma once

class Box;
class Random;

class ProtonDisorder {

  /**
   * Generates proton coordinates given a box of oxygen atoms.  Oxygen atoms
   * are considered neighbors (needing a hydrogen between them) if their
   * distance is less than drNbrOO.  The hydrogens are placed such that they
   * are bondLengthOH distance away from their owning oxygen and have an angle
   * of bondAngleHOH.  If offsetM is not 0, then an M coordinate is also
   * determined along the water bisector.
   */
  public:
    static double** go(Box& box, Random& rand, const double drNbrOO, const double bondLengthOH, const double bondAngleHOH, const double offsetM);
};
