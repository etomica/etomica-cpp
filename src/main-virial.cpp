#include <stdio.h>

#include "potential-master.h"
#include "integrator.h"
#include "potential.h"
#include "move-virial.h"
#include "box.h"
#include "meter-virial.h"
#include "data-sink.h"
#include "random.h"
#include "util.h"
#include "cluster.h"

int main(int argc, char** argv) {
  int order = 2;
  double temperature = 1.0;
  long steps = 10000000;

  Random rand;
  printf("random seed: %d\n", rand.getSeed());

  PotentialLJ plj(1,1,TRUNC_NONE, 1.0/0.0);
  PotentialHS pHS(1.5);
  SpeciesList speciesList;
  speciesList.add(new SpeciesSimple(1,1));
  Box box(speciesList);
  box.setBoxSize(1,1,1);
  box.setNumMolecules(0, order);
  PotentialMasterVirial potentialMaster(speciesList, box);
  potentialMaster.setPairPotential(0, 0, &plj);
  PotentialMasterVirial potentialMasterHS(speciesList, box);
  potentialMasterHS.setPairPotential(0, 0, &pHS);
  IntegratorMC integrator(potentialMaster, rand);
  Cluster clusterLJ(potentialMaster, temperature, 0, true);
  Cluster clusterHS(potentialMasterHS, temperature, 0, false);
  MCMoveDisplacementVirial move(box, potentialMaster, rand, 0.2, clusterLJ);
  integrator.addMove(&move, 1);
  integrator.setTemperature(temperature);
  integrator.reset();
  integrator.doSteps(steps/10);
  MeterVirialDirect meter(clusterLJ, clusterHS);
  DataPump pumpVirial(meter, 1);
  integrator.addListener(&pumpVirial);

  double t1 = getTime();
  integrator.doSteps(steps);
  double t2 = getTime();

  double** stats = ((Average*)pumpVirial.getDataSink(0))->getStatistics();
  printf("target avg: %f  err: %f  cor: %f\n", stats[0][AVG_AVG], stats[0][AVG_ERR], stats[0][AVG_ACOR]);
  printf("reference avg: %f  err: %f  cor: %f\n", stats[1][AVG_AVG], stats[1][AVG_ERR], stats[1][AVG_ACOR]);
  printf("time: %4.3f\n", t2-t1);
}

