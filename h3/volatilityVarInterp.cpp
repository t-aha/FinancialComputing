#include "home3/home3.hpp"

using namespace cfl;
using namespace std;

cfl::Function
prb::volatilityVarInterp(const std::vector<double> &rTimes,
                        const std::vector<double> &rVols,
                        double dInitialTime, cfl::Interp &rInterp)
    {
    std::vector<double> tims (rTimes.size () + 1);
    tims.front () = dInitialTime;
    std::copy (rTimes.begin (), rTimes.end (), tims.begin () + 1);
    vector<double> vs;

    int n = tims.size();

    for (int i = 0; i < n; i++){
        vs.push_back(rVols[i]*rVols[i]*(tims[i+1]-dInitialTime));
        if (i == 0) {vs.push_back(vs[0]);};
    }

    rInterp.assign (tims.begin (), tims.end (), vs.begin ());
    Function vFunc = rInterp.interp ();

    Function out ([dInitialTime, vFunc, rVols] (double dT) {
        if (dT==dInitialTime) {return rVols[0];};
        return sqrt(vFunc(dT)/(dT-dInitialTime));
    });
    return Function(out, dInitialTime);
    }