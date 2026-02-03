//
// Created by skwab on 12/6/2025.
//

#include "TraPPEParams.h"
#include "unit.h"
#include <cmath>

TraPPEParams::TraPPEParams(ChemForm chemForm) : chemForm(chemForm), species(makeSpecies(chemForm)){
    if (chemForm == CO2)
    {
        //TraPPE Parameters
        isFlex = true;
        isPolar = true;
        double bondLengthCO = 1.160; // Angstrom
        double sigmaC = 2.800; // Angstrom
        double epsilonC = Kelvin::toSim(27.0);
        double qC = Electron::toSim(0.700);
        double sigmaO = 3.050; // Angstrom
        double epsilonO = Kelvin::toSim(79.0);
        double qO = Electron::toSim(-0.350);
        double k_r = Kelvin::toSim(1000000);
        double kCOO = Kelvin::toSim(62500);
        double thetaCOO = Degree::toSim(180);

        sigma = {sigmaC, sigmaO};
        epsilon = {epsilonC, epsilonO};
        charge = {qC, qO};
        k_b = {k_r};
        r_eq = {bondLengthCO};
        bonding = {{1, 2}, {0}, {0}};
        vector<int*> temp;
        int* p1 = (int*)malloc(2*sizeof(int));
        int* p2 = (int*)malloc(2*sizeof(int));

        p1[0] = 0;
        p1[1] = 1;
        p2[0] = 0;
        p2[1] = 2;
        temp.push_back(p1);
        temp.push_back(p2);
        pairs.push_back(temp);
        k_theta = {kCOO};
        theta_eq = {thetaCOO};
        int* t = (int*)malloc(3*sizeof(int));
        t[0] = 1;
        t[1] = 0;
        t[2] = 2;
        triplets.push_back(t);

        species.addAtomType(12, 1);
        species.addAtomType(16, 2);
        species.setAtomPosition(0, 0, 0, 0);
        species.setAtomPosition(1, -bondLengthCO, 0, 0);
        species.setAtomPosition(2, bondLengthCO, 0, 0);

        speciesList.add(&species);

    }
    else if (chemForm == propane)
    {
        //TraPPE-UA
        isFlex = true;
        //TraPPE Parameters
        double bondLengthCHxCHy = 1.54; // Angstrom
        double thetaCCH = Degree::toSim(114);
        double sigmaCH3 = 3.75; // Angstrom
        double epsilonCH3 = Kelvin::toSim(98);
        double qCH3 = Electron::toSim(0.0);
        double sigmaCH2 = 3.95; // Angstrom
        double epsilonCH2 = Kelvin::toSim(46);
        double qCH2 = Electron::toSim(0.0);
        double kCCC = Kelvin::toSim(62500);
        double yy = bondLengthCHxCHy*sin(thetaCCH) ;
        //Construct Arrays
        sigma = {sigmaCH3, sigmaCH2};
        epsilon = {epsilonCH3, epsilonCH2};
        charge = {qCH3, qCH2};
        k_theta = {kCCC};
        theta_eq = {thetaCCH};
        int* t = (int*)malloc(3*sizeof(int));
        t[0] = 0;
        t[1] = 2;
        t[2] = 1;
        triplets.push_back(t);
        bonding = {{2}, {2}, {0, 1}};

        species.addAtomType(12, 2);
        species.addAtomType(12, 1);

        //Get Coordinates
        double posCavg = (2 * bondLengthCHxCHy  - bondLengthCHxCHy * cos(thetaCCH))/3;
        species.setAtomPosition(0, 0 - posCavg, -yy / 3, 0);
        species.setAtomPosition(2, bondLengthCHxCHy - posCavg, -yy / 3, 0);
        species.setAtomPosition(1, bondLengthCHxCHy - bondLengthCHxCHy * cos(thetaCCH) - posCavg, 2 * yy / 3, 0);

        speciesList.add(&species);

    }


}

SpeciesSimple TraPPEParams::makeSpecies(ChemForm chemForm)
{
    if (chemForm == CO2)
    {
        return SpeciesSimple(3, 2);
    }
    if (chemForm == propane)
    {
        return SpeciesSimple(3, 2);
    }
}



