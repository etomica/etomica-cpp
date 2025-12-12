
#include "P2PotentialGroupBuilder.h"

#include "potential-master.h"


PotentialMasterVirial P2PotentialGroupBuilder::builder(const SpeciesList &species_list, int numAtomTypes, vector<double> sigma, vector<double> epsilon, vector<double> charge, Box &box){
    PotentialMasterVirial potentialGroup(species_list, box);
    double sigmaHC = 0.3;
    for(int i = 0; i < numAtomTypes; i++) {
        for (int j = 0; j < numAtomTypes; j++) {
            double sigmaij = sigma[i];
            double epsilonij = epsilon[i];
            double qiqj = charge[i]*charge[j];
            PotentialLJ p2LJ(epsilonij, sigmaij, TRUNC_NONE, 1.0/0.0);
            // potentialGroup.setMoleculePairPotential(0, 0, &p2LJ); //Set Atomic Potential
            PotentialChargeBare p2ES(qiqj, sigmaHC, 1.0/0.0);
            P2SoftSphericalSum p2(p2LJ, p2ES);
            potentialGroup.setPairPotential(i, j, &p2); //Set Atomic Potential
        }


    }
    return potentialGroup;

}


