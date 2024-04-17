#include "home4/home4.hpp"

using namespace cfl;
using namespace std;

cfl::Function
prb::discountNelsonSiegelFit (const std::vector<double> &rTimes,
                            const std::vector<double> &rDF,
                            double dLambda, double dInitialTime,
                            cfl::Function &rErr,
                            cfl::FitParam &rParam)

    {
    
    Function basis1 ([dLambda, dInitialTime] (double dT) {
        if (dT==dInitialTime) {return 1.0;};
        return (1.-exp(-dLambda*(dT-dInitialTime)))/(dLambda*(dT-dInitialTime));
    });

    Function basis2 ([dLambda, dInitialTime] (double dT) {
        if (dT==dInitialTime) {return 1.0;};
        return (1.-exp(-dLambda*(dT-dInitialTime)))/(dLambda*(dT-dInitialTime))-exp(-dLambda*(dT-dInitialTime));
    });

    Fit rFit = prb::linear(vector<Function> {Function(1.0), basis1, basis2});
    Function ftr = prb::discountYieldFit(rTimes, rDF, dInitialTime, rFit, rErr);
    
    rParam = rFit.param();

    return ftr;
    }