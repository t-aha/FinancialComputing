#include "home5/home5.hpp"

using namespace cfl;
using namespace std;

cfl::Function 
prb::discountRateFit (const std::vector<double> &rPeriods,
                            const std::vector<double> &rRates,
                            double dInitialTime, cfl::Fit &rFit,
                            cfl::Function &rErr)
{
    std::vector<double> ins(rRates.size());
    std::vector<double> tim(rRates.size());
    std::transform(rPeriods.begin(), rPeriods.end(), tim.begin(),
                   [](double value) { return value + 1.0; });

    std::transform(rRates.begin(), rRates.end(), rPeriods.begin(), ins.begin(),
                   [dInitialTime](double rate, double t) 
                   { return -log(1/(1.+rate*t)) / (t); });

    rFit.assign(begin(tim), end(tim), 
        begin(std::valarray<double> (ins.data(), ins.size())));

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