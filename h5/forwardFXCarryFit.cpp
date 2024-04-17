#include "home5/home5.hpp"

using namespace cfl;
using namespace std;

cfl::Function 
prb::forwardFXCarryFit (double dSpotFX, const std::vector<double> &rTimes,
                   const std::vector<double> &rDomesticDiscountFactors,
                   const std::vector<double> &rForeignDiscountFactors,
                   double dInitialTime, cfl::Fit &rFit, cfl::Function &rErr)
{
    vector<double> ins(rDomesticDiscountFactors.size());
    std::transform(rForeignDiscountFactors.begin(), 
                   rForeignDiscountFactors.end(), 
                   rDomesticDiscountFactors.begin(), 
                   ins.begin(), std::divides<double>());
    
    std::transform(ins.begin(), ins.end(), rTimes.begin(), ins.begin(),
                   [dInitialTime](double rate, double t) 
                   { return log(rate) / (t - dInitialTime); });

    rFit.assign(begin(rTimes), end(rTimes), 
        begin(std::valarray<double> (ins.data(), ins.size())));
    

    Function ftr = rFit.fit();
    rErr = rFit.err();

    Function out ([dSpotFX, dInitialTime, ftr] (double dT) {
        return dSpotFX*exp(ftr(dT)*(dT-dInitialTime));
    });

    Function out_err ([dSpotFX, dInitialTime, ftr, rErr] (double dT) {
        return dSpotFX*exp(ftr(dT)*(dT-dInitialTime))*(dT-dInitialTime)*rErr(dT);
    });
    
    rErr = out_err;
    return Function(out, dInitialTime);
}