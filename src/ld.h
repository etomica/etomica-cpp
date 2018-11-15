class Potential;

class LatticeDynamics {
#ifdef LATTICE_DYNAMICS
  private:
    int numCells[3];
    double density;
    double rCut;
    double lrcFac;
    Potential* potential;
    bool doLRC;
    int status;
    int nBasis;
    double **basis;
    double cellSize[3];
    double strain[3];
    double ***selfSum;
    int wCount;
    double** waveVectors;
    std::complex<double> ***matrix;
    int* wvCount;
    int latticeShells[3];
    long long doneXYZ;

    double rNext;
    int xcellNext, ycellNext, zcellNext;
    int wvNext;

    double logSum;
    bool unstable;
    double uLat;
  public:
    LatticeDynamics(double d, Potential *p, bool doLRC, int nBasis);
    ~LatticeDynamics();
    void setBasis(int i, double x, double y, double z);
    void setUnitCell(double x, double y, double z);
    void setStrain(double x, double y, double z);
    void setNumCells(int x, int y, int z);
    int getStatus() {return status;}
    void setup();
    long long countLS();
    long long goLS(int nMax);
    int doSelfSum();
    int goEVD(int nwv);
    int getWaveVectorCount() {return wCount;}
    double getUnstable() {return unstable;}
    double getLogSum() {return logSum;}
    double getU();
#endif
};
