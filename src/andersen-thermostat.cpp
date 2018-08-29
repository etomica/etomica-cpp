#include "integrator.h"
#include "box.h"
#include "math.h"

AndersenThermostat::AndersenThermostat(IntegratorMD& i, Random& r) : IntegratorListener(), integratorMD(i), random(r), mode(0), interval(-1), intervalCountdown(-1), pRandomize(0), zeroMomentum(false) {}

AndersenThermostat::~AndersenThermostat() {}

void AndersenThermostat::setSingle(int i) {
  mode = ANDERSEN_SINGLE;
  intervalCountdown = interval = i;
}

void AndersenThermostat::setPartial(double p) {
  mode = ANDERSEN_PARTIAL;
  pRandomize = p;
}

void AndersenThermostat::setFull(int i, bool zm) {
  mode = ANDERSEN_FULL;
  intervalCountdown = interval = i;
  zeroMomentum = zm;
}

int AndersenThermostat::getMode() {
  return mode;
}

void AndersenThermostat::stepFinished() {
  Box& box = integratorMD.getBox();
  int numAtoms = box.getNumAtoms();
  switch (mode) {
    case ANDERSEN_SINGLE:
      intervalCountdown--;
      if (intervalCountdown==0) {
        integratorMD.randomizeVelocity(random.nextInt(numAtoms));
      }
      break;
    case ANDERSEN_PARTIAL:
      // approach here works efficiently if frequency is low
      for (int i=0; i<numAtoms; i++) {
        // probabilty we will randomize an atom
        double p = pow(1-pRandomize, numAtoms-i);
        double x = random.nextDouble32();
        if (x > p) break;
        // we may end up randomizing an atom more than once
        integratorMD.randomizeVelocity(random.nextInt(numAtoms));
      }
      break;
    case ANDERSEN_FULL:
      intervalCountdown--;
      if (intervalCountdown==0) integratorMD.randomizeVelocities(zeroMomentum);
      break;
    default:
      break;
      // disabled
  }
}
