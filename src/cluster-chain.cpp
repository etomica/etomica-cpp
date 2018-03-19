#include "cluster.h"

ClusterChain::ClusterChain(PotentialMasterVirial &pm, double t, double cFac, double rFac, bool cached) : Cluster(pm.getBox().getTotalNumMolecules(),1,cached), potentialMaster(pm), beta(1/t), chainFac(cFac), ringFac(rFac) {
}

ClusterChain::~ClusterChain() {
}

#define NF1  (1 << (numMolecules-1))
const double* ClusterChain::getValues() {
  if (useCache && !cacheDirty) return values;
  cacheDirty = false;
  const int n = numMolecules;
  double fValues[n][n];
  for(int iMol1=0; iMol1<n; iMol1++){
    moleculePair[0] = iMol1;
    for(int iMol2=iMol1+1; iMol2<n; iMol2++){
      moleculePair[1] = iMol2;
      double u12;
      potentialMaster.computeMolecules(moleculePair, 2, u12);
      fValues[iMol1][iMol2] = fValues[iMol2][iMol1] = fabs(exp(-beta*u12)-1);
    }
  }

  double nC[n-1][NF1];
  //nC[m][i] is the number of chains beginning at last vertex and ending at m, traversing all points in i

  //Start with all pairwise paths from last vertex to each other vertex
  for(int m=0; m<n-1; m++) {
    nC[m][(1 << m)] = fValues[m][n - 1];
  }

  //All other paths
  for(int i=1; i<NF1; i++) {//i excludes the last bit, which is implicit in all partitions
    //the following two loops generate all pairs formed by each bit in i with each bit not in i

    //loop over bits not in i; start with full complement of i (i^(NF1-1)), and in each iteration
    //get lowest bit (im=(iC&-iC)) and strip it from complement (iC^=im) until complement is empty (iC=0)
    for(int iC=i^(NF1-1), im=(iC&-iC); iC>0; iC^=im,im=(iC&-iC)) {
      int m = log2(im);
      int iim = i|im;
      nC[m][iim] = 0;
      //loop over bits in i, in same manner as loop over complement
      for (int it = i, ik = (it & -it); ik > 0; it ^= ik, ik = (it & -it)) {
        int k = log2(ik);
        nC[m][iim] += fValues[m][k] * nC[k][i];
      }
    }//end for(iC)
  }//end for(i)

  double ringValue = 0;
  double chainValue = 0;

  if (ringFac != 0.0) {
    for (int m = 0; m < n - 1; m++) {
      ringValue += nC[m][NF1-1] * fValues[m][n-1];
    }
  }

  if (chainFac != 0.0) {

    //Sum chains in which last (n-1) vertex is not a leaf.
    //Consider all partitions, counting paths beginning in one partition and ending in its complement 
    //Use same looping structure as employed above
    for (int iS = 1; iS < NF1; iS += 2) {//keep 1 in iS-partition to prevent double counting
      //loop over bits not in iS
      int iSComp = iS^(NF1-1);
      for (int iC = iSComp, im = (iC & -iC); iC > 0; iC ^= im, im = (iC & -iC)) {
        int m = log2(im);
        //loop over bits in iS
        for (int it = iS, ik = (it & -it); ik > 0; it ^= ik, ik = (it & -it)) {
          int k = log2(ik);
          chainValue += nC[m][iSComp] * nC[k][iS];
        }
      }
    }

    //Sum chains where last (n-1) vertex is a leaf
    for (int m = 0; m < n - 1; m++) {
      chainValue += nC[m][NF1-1];
    }

  }//end if(chainFrac)

  oldValues[0] = values[0];
  values[0] = chainFac*chainValue + ringFac*ringValue;

  return values;
}
