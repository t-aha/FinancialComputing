#include "home4/home4.hpp"

using namespace cfl;
using namespace std;

cfl::Function
prb::discountYieldFit (const std::vector<double> &rTimes,
                    const std::vector<double> &rDF,
                    double dInitialTime, cfl::Fit &rFit,
                    cfl::Function &rErr)
    {

    vector<double> gs;

    for (size_t i = 0; i < rTimes.size(); i++) {
        gs.push_back(-log(rDF[i])/(rTimes[i]-dInitialTime));
    }
    rFit.assign(begin(rTimes), end(rTimes), begin(std::valarray<double> (gs.data(), gs.size())));
    
    Function ftr = rFit.fit();
    
    rErr = rFit.err();
    Function out ([rDF, dInitialTime, ftr] (double dT) {
        return exp(-ftr(dT)*(dT-dInitialTime));
    });

    Function out_err ([rDF, dInitialTime, ftr, rErr] (double dT) {
        return exp(-ftr(dT)*(dT-dInitialTime))*(dT-dInitialTime)*rErr(dT);
    });
    rErr = out_err;

    return Function(out, dInitialTime);
    }