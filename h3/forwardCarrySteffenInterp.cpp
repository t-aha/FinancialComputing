#include "home3/home3.hpp"

using namespace cfl;
using namespace std;

cfl::Function
prb::forwardCarrySteffenInterp(double dSpot,
                        const std::vector<double> &rDeliveryTimes,
                        const std::vector<double> &rForwardPrices,
                        double dInitialTime)
    {
    std::vector<double> tims (rDeliveryTimes.size () + 1);
    tims.front () = dInitialTime;
    std::copy (rDeliveryTimes.begin (), rDeliveryTimes.end (), tims.begin () + 1);
    vector<double> qs;
    int n = tims.size();
    for (int i = 0; i < n; i++){
        qs.push_back(log(rForwardPrices[i]/dSpot)/(tims[i+1]-dInitialTime));
        if (i == 0) {qs.push_back(qs[0]);};
    }
    cfl::Interp rInterp = prb::steffen();
    rInterp.assign (tims.begin (), tims.end (), qs.begin ());
    Function qFunc = rInterp.interp ();

    Function out ([dSpot, dInitialTime, qFunc] (double dT) {
        return dSpot*exp(qFunc(dT)*(dT-dInitialTime));
    });

    return Function(out, dInitialTime);
    
    }