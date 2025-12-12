//
// Created by skwab on 12/6/2025.
//

#ifndef MYPROJECT_P2POTENTIALGROUPBUILDER_H
#define MYPROJECT_P2POTENTIALGROUPBUILDER_H
#include "box.h"
#include "potential-master.h"

#include "species.h"
#include "TraPPEParams.h"
class PotentialMasterVirial;

class P2PotentialGroupBuilder {

public:

    static PotentialMasterVirial builder(const SpeciesList &species_list, int numAtomTypes, vector<double> sigma, vector<double> epsilon, vector<double> charge, Box
                                         &box);
};


#endif //MYPROJECT_P2POTENTIALGROUPBUILDER_H