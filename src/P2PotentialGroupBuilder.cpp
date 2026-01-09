
#include "P2PotentialGroupBuilder.h"

#include "potential-master.h"


PotentialMasterVirial P2PotentialGroupBuilder::builder(PotentialMasterVirial& potentialGroup, const SpeciesList &species_list, int numAtomTypes, vector<double> sigma, vector<double> epsilon, vector<double> charge, Box &box){
    double sigmaHC = 0.3;
    for(int i = 0; i < numAtomTypes; i++) {
        for (int j = 0; j <= i; j++) {
            double sigmaij = 0.5*(sigma[i]+sigma[j]);
            double epsilonij = sqrt(epsilon[i]*epsilon[j]);
            double qiqj = charge[i]*charge[j];
            PotentialLJ *p2LJ = new PotentialLJ(epsilonij, sigmaij, TRUNC_NONE, 1.0/0.0);
            if (qiqj == 0) {
                potentialGroup.setPairPotential(i, j, p2LJ); //Set Atomic Potential

            }
            else {
                PotentialChargeBare *p2ES = new PotentialChargeBare(qiqj, sigmaHC, 1.0/0.0);
                P2SoftSphericalSum *p2 = new P2SoftSphericalSum(*p2LJ, *p2ES);
                potentialGroup.setPairPotential(i, j, p2); //Set Atomic Potential

            }
            // printf("%d, %d, %f, %f, %f\n", i, j, sigmaij, epsilonij, qiqj);
        }


    }
    return potentialGroup;

}


