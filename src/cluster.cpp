#include "cluster.h"

Cluster::Cluster(PotentialMasterVirial &pm, double t, int nd, bool cached) : potentialMaster(pm), numMolecules(pm.getBox().getTotalNumMolecules()), beta(1/t), nDer(nd), useCache(cached), cacheDirty(true) {
  values = new double[nDer+1];
  oldValues = new double[nDer+1];
  binomial = new int*[nDer+1];
  int factorial[nDer+1];
  factorial[0] = factorial[1] = 1;
  for (int m=2; m<=nDer; m++) factorial[m] = m*factorial[m-1];
  for (int m=0; m<=nDer; m++) {
    binomial[m] = new int[m+1];
    for (int l=0; l<=m; l++) {
      binomial[m][l] = factorial[m]/(factorial[l]*factorial[m-l]);
    }
  }
  prefac = -(numMolecules-1.0)/factorial[numMolecules];
}

Cluster::~Cluster() {
  delete[] values;
  delete[] oldValues;
  for (int i=0; i<nDer; i++) {
    delete[] binomial[i];
  }
  delete[] binomial;
}

void Cluster::setCachingEnabled(bool enabled) {
  useCache = enabled;
}

int Cluster::numValues() {
  return nDer+1;
}

void Cluster::trialNotify() {
  for (int m=0; m<=nDer; m++) oldValues[m] = values[m];
  cacheDirty = true;
}

void Cluster::trialRejected() {
  for (int m=0; m<=nDer; m++) values[m] = oldValues[m];
  cacheDirty = false;
}

#define NF  (1 << numMolecules)
const double* Cluster::getValues() {
  if (useCache && !cacheDirty) return values;
  cacheDirty = false;

  double fQ[NF][nDer+1], fC[NF][nDer+1];
  double fA[NF][nDer+1], fB[NF][nDer+1];

  for(int iMol1=0; iMol1<numMolecules; iMol1++){
    int i = 1<<iMol1;
    fQ[i][0] = 1.0;
    for (int m=1; m<=nDer; m++) fQ[i][m] = 0;
    moleculePair[0] = iMol1;
    for(int iMol2=iMol1+1; iMol2<numMolecules; iMol2++){
      moleculePair[1] = iMol2;
      double u12;
      potentialMaster.computeMolecules(moleculePair, 2, u12);
      fQ[i|(1<<iMol2)][0] = exp(-beta*u12);
    }
  }

  //generate all partitions and compute
  for (int i=3; i<NF; i++){
    int j = i & -i; //lowest bit in i
    if (i == j) continue; //1-point set
    int k = i & ~j;
    if (k == (k & -k)) {
      // 2-point set
      if (fQ[i][0] == 0) {
        for (int m=1; m<=nDer; m++) fQ[i][m] = 0;
        continue;
      }
      double c = log(fQ[i][0])/beta;
      for (int m=1; m<=nDer; m++) fQ[i][m] = fQ[i][m-1]*c;
      continue;
    }
    fQ[i][0] = fQ[k][0];
    if (fQ[i][0] == 0) {
      for (int m=1; m<=nDer; m++) fQ[i][m] = 0;
      continue;
    }
    for (int l = (j<<1); l<i; l=(l<<1)){
      if ( (l&i) == 0 ) continue;
      fQ[i][0] *= fQ[l|j][0];
    }
    if (fQ[i][0] == 0) {
      for (int m=1; m<=nDer; m++) fQ[i][m] = 0;
      continue;
    }
    double c = log(fQ[i][0])/beta;
    for (int m=1; m<=nDer; m++) fQ[i][m] = fQ[i][m-1]*c;
  }

  //Compute the fC's
  for (int i=1; i<NF; i++){
    for (int m=0; m<=nDer; m++) fC[i][m] = fQ[i][m];
    int iLowBit = i & -i;
    int inc = iLowBit<<1;
    for (int j=iLowBit; j<i; j+=inc){
      int jComp = i & ~j;
      while ((j|jComp) != i && j<i){
        int jHighBits = j^iLowBit;
        int jlow = jHighBits & -jHighBits;
        j += jlow;
        jComp = (i & ~j);
      }
      if (j==i) break;
      for (int m=0; m<=nDer; m++) {
        for (int l=0; l<=m; l++) {
          fC[i][m] -= binomial[m][l] * fC[j][l] * fQ[jComp][m-l];
        }
      }
    }
  }

  //find fA1
  for (int i=2; i<NF; i+=2){
    //all even sets don't contain 1
    for (int m=0; m<=nDer; m++) fB[i][m] = fC[i][m];
  }

  for (int m=0; m<=nDer; m++) {
    fA[1][m] = 0;
    fB[1][m] = fC[1][m];
  }
  for (int i=3; i<NF; i+=2){
    //every set will contain 1
    for (int m=0; m<=nDer; m++) {
      fA[i][m] = 0;
      fB[i][m] = fC[i][m];
    }
    int ii = i - 1;//all bits in i but lowest
    int iLow2Bit = (ii & -ii);//next lowest bit
    int jBits = 1 | iLow2Bit;
    if (jBits == i) continue;

    int iii = ii ^ iLow2Bit; //i with 2 lowest bits off
    int jInc = (iii & -iii);//3rd lowest bit, alsso increment for j
    for (int j=jBits; j<i; j+=jInc){//sum over partitions of i containing j Bits
      int jComp = (i & ~j);//subset of i complementing j
      while ((j|jComp) != i && j<i){
        int jHighBits = j ^ jBits;
        int jlow = jHighBits & -jHighBits;
        j += jlow;
        jComp = (i & ~j);
      }
      if (j==i) break;
      for (int m=0; m<=nDer; m++) {
        for (int l=0; l<=m; l++) {
          fA[i][m] += binomial[m][l]*fB[j][l] * fC[jComp|1][m-l];
        }
      }
    }
    //remove from B graphs that contain articulation point 0.
    for (int m=0; m<=nDer; m++) fB[i][m] -= fA[i][m];
  }

  for (int v=1; v<numMolecules; v++){
    int vs1 = 1<<v;
    for (int i=vs1+1; i<NF; i++){
      for (int m=0; m<=nDer; m++) fA[i][m] = 0;

      if ( (i & vs1) == 0 ) continue;
      int iLowBit = (i & -i);
      if ( iLowBit == i ) continue;

      int jBits;
      int ii = i ^ iLowBit;
      int iLow2Bit = (ii & -ii);
      if ( iLowBit!=vs1 && iLow2Bit!=vs1 ){
        //v is not in the lowest 2 bits
        jBits = iLowBit | vs1;
        //we can only increment by the 2nd lowest
        int jInc = iLow2Bit;

        for (int j=jBits; j<i; j+=jInc){
          if ( (j&jBits) != jBits ){
            j |= vs1;
            if (j==i) break;
          }
          int jComp = i & ~j;
          while ((j|jComp) != i && j<i){
            int jHighBits = j^jBits;
            int jlow = jHighBits & -jHighBits;
            j += jlow;
            j |= vs1;
            jComp = (i & ~j);
          }
          if (j==i) break;
          for (int m=0; m<=nDer; m++) {
            for (int l=0; l<=m; l++) {
              fA[i][m] += binomial[m][l]*fB[j][l] * (fB[jComp|vs1][m-l] + fA[jComp|vs1][m-l]);
            }
          }
        }
      }else{
        //lowest 2 bits contain v
        jBits = iLowBit | iLow2Bit;
        if (jBits == i) continue; // no bits left jComp

        int iii = ii ^ iLow2Bit;
        int jInc = ( iii & -iii);

        //at this point jBits has (lowest bit + v)
        for (int j=jBits; j<i; j+=jInc){//sum over partitions of i
          int jComp = i & ~j;
          while ((j|jComp) != i && j<i){
            int jHighBits = j^jBits;
            int jlow = jHighBits & -jHighBits;
            j += jlow;
            jComp = (i & ~j);
          }
          if (j==i) break;
          for (int m=0; m<=nDer; m++) {
            for (int l=0; l<=m; l++) {
              fA[i][m] += binomial[m][l]*fB[j][l] * (fB[jComp|vs1][m-l] + fA[jComp|vs1][m-l]);
            }
          }
        }
      }
      for (int m=0; m<=nDer; m++) fB[i][m] -= fA[i][m];
    }
  }
  for (int m=0; m<=nDer; m++) {
    oldValues[m] = values[m];
  }

  values[0] = prefac*fB[NF-1][0];
  if (values[0] != 0 && fabs(values[0]) < 1.E-12 ){
    for (int m=1; m<=nDer; m++) {
      values[m] = 0;
    }
    return values;
  }
  for (int m=1; m<=nDer; m++) {
    values[m] = prefac*fB[NF-1][m];
  }
  return values;
}
