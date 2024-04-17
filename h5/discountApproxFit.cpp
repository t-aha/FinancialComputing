#include "home5/home5.hpp"

using namespace cfl;
using namespace std;

cfl::Function 
prb::discountApproxFit (const std::vector<double> &rDiscountTimes,
                        const std::vector<double> &rDiscountFactors,
                        double dInitialTime, cfl::Fit &rFit,
                        cfl::Function &rErr)
{
    std::vector<double> ins(rDiscountFactors.size());
    std::vector<double> w(rDiscountFactors.size());

    std::transform(rDiscountFactors.begin(), rDiscountFactors.end(), 
                rDiscountTimes.begin(), ins.begin(),
                   [dInitialTime](double rate, double t) 
                   { return -log(rate) / (t-dInitialTime); });

    std::transform(rDiscountFactors.begin(), rDiscountFactors.end(), 
                rDiscountTimes.begin(), w.begin(), 
                [dInitialTime](double rate, double t) 
                   { return rate * rate * (t-dInitialTime)*(t-dInitialTime);});

    rFit.assign(begin(rDiscountTimes), end(rDiscountTimes), 
            begin(std::valarray<double> (ins.data(), ins.size())),
            begin(std::valarray<double> (w.data(), w.size())));

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