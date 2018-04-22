#include "move.h"
#include "meter.h"
#include "species.h"

class MCMoveVolume : public MCMove {
  private:
    double lnScale, scale;
    double pressure;
    double vOld;
    double uOld, uNew;
    SpeciesList& speciesList;
    Meter& oldMeterPE;
    vector<PotentialCallback*> callbacks;
    PotentialCallbackEnergy pce;
    void scaleVolume(double s);

  public:
    MCMoveVolume(Box& box, PotentialMaster& potentialMaster, Random& random, double pressure, double stepSize, SpeciesList& speciesList, Meter& oldMeterPE);
    ~MCMoveVolume() {}

    virtual bool doTrial();
    virtual double getChi(double temperature);
    virtual void acceptNotify();
    virtual void rejectNotify();
    virtual double energyChange();
};
