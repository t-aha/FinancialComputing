#include "home4/home4.hpp"

using namespace cfl;
using namespace std;

cfl::Function
prb::discountHullWhiteFit(const std::vector<double> &rDiscountTimes,
                                   const std::vector<double> &rDiscountFactors,
                                   double dLambda, double dInitialTime,
                                   cfl::Function &rErr, cfl::FitParam &rParam) {
    std::vector<double> tims (rDiscountTimes.size () + 1);
    tims.front () = dInitialTime;
    std::copy (rDiscountTimes.begin (), rDiscountTimes.end (), tims.begin () + 1);

    vector<double> gs;

    for (size_t i = 0; i < rDiscountTimes.size(); i++) {
        gs.push_back(-log(rDiscountFactors[i])/(rDiscountTimes[i]-dInitialTime));
    }
    
    Function basisF ([dLambda, dInitialTime] (double dT) {
        if (dT==dInitialTime) {return 1.0;};
        return (1.-exp(-dLambda*(dT-dInitialTime)))/(dLambda*(dT-dInitialTime));
    });
    cfl::Fit rFit = cfl::NFit::linear_regression(basisF);
    rFit.assign(begin(rDiscountTimes), end(rDiscountTimes), begin(std::valarray<double> (gs.data(), gs.size())));
    rErr = rFit.err();
    Function ftr = rFit.fit();
    
    Function out ([rDiscountFactors, dInitialTime, ftr] (double dT) {
        return exp(-ftr(dT)*(dT-dInitialTime));
    });

    Function out_err ([rDiscountFactors, dInitialTime, rErr, ftr] (double dT) {
        return exp(-ftr(dT)*(dT-dInitialTime))*(dT-dInitialTime)*rErr(dT);
    });
    rErr = out_err;
    rParam = rFit.param();

    return Function(out, dInitialTime);
}