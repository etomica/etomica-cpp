//
// Created by raghu on 4/13/2026.
//
#include "cluster.h"
#include "potential-master.h"
#include "potential-molecular.h"


ClusterMBX::ClusterMBX(PotentialMasterVirial &pm, double t, int nd, bool cached, PotentialMolecularMBnrg3body &p3body) :
    ClusterVirial(pm, t, nd, cached), p3body(p3body)
{
}

ClusterMBX::~ClusterMBX() {
}

double ClusterMBX::computeExtraFQ(int i)
{
    vector<int> molecules;
    int l = 0;
    for (int a=0; a<numMolecules; a++) {
        if ((i & (1<<a)) != 0) {
            molecules.push_back(a);
            l++;
        }
    }
    if (l == 3)
    {
        //3 body stuff
        Box *box = potentialMaster.getBoxP();
        double u3 = p3body.u(*box, molecules[0], molecules[1], molecules[2]);
        return exp(-beta*u3);

    }
    return 1;
}
