//
// Created by skwab on 12/6/2025.
//

#include "TraPPEParams.h"
#include "TraPPEParams.h"

#include "P2PotentialGroupBuilder.h"
#include "unit.h"

TraPPEParams::TraPPEParams() : species(3, 2), chemForm(CO2) {
            //TraPPE Parameters
            numAtomTypes = 2;
            double bondLengthCO = 1.160; // Angstrom
            double sigmaC = 2.800; // Angstrom
            double epsilonC = Kelvin::toSim(27.0);
            double qC = Electron::toSim(0.700);
            double sigmaO = 3.050; // Angstrom
            double epsilonO = Kelvin::toSim(79.0);
            double qO = Electron::toSim(-0.350);
            sigma = {sigmaC, sigmaO};
            epsilon = {epsilonC, epsilonO};
            charge = {qC, qO};
            species.addAtomType(12, 1);
            species.addAtomType(16, 2);
            species.setAtomPosition(0, 0, 0, 0);
            species.setAtomPosition(1, -bondLengthCO, 0, 0);
            species.setAtomPosition(2, bondLengthCO, 0, 0);

            speciesList.add(&species);
            // params = {numAtomTypes, sigma, epsilon, charge};

    }


// void TraPPEParams::buildPotentials() {
//     switch (chemForm) {
//         case CO2:
//             params params = {numAtomTypes, sigma, epsilon, charge};
//             potentialGroup = P2PotentialGroupBuilder::builder(speciesList, params, box);
//
//     }
// }


