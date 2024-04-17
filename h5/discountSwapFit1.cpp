#include "home5/home5.hpp"

using namespace cfl;
using namespace std;

cfl::Function 
prb::discountSwapFit (const std::vector<double> &rSwapRates,
                               double dPeriod, double dInitialTime,
                               cfl::Fit &rFit, cfl::Function &rErr)
{
    std::vector<double> tim (rSwapRates.size ());
    std::vector<double> ins (rSwapRates.size ());
    double tmp = 0.0;

    for (size_t i = 0; i < rSwapRates.size (); ++i)
    {
        tim[i] = dInitialTime + (i+1) * dPeriod;

        tmp *= rSwapRates [i];
        ins[i] = (1 -tmp) / (1 + rSwapRates[i]*dPeriod);
        tmp = tmp / rSwapRates[i] + dPeriod*ins [i];
    }
    
    std::transform(ins.begin(), ins.end(), tim.begin(), ins.begin(),
                   [dInitialTime](double rate, double t) 
                   { return log(rate) / (t - dInitialTime); });
    
    rFit.assign (tim.begin (), tim.end (), ins.begin());
    
    Function ftr = rFit.fit();
    rErr = rFit.err();

    Function out ([dInitialTime, ftr] (double dT) {
        return exp(-ftr(dT)*(dT-dInitialTime));
    });

    Function out_err ([dInitialTime, ftr, rErr] (double dT) {
        return exp(-ftr(dT)*(dT-dInitialTime))*(dT-dInitialTime)*rErr(dT);
    });
    
    rErr = out_err;
    return Function(out, dInitialTime);
}