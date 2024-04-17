#include "home4/home4.hpp"

using namespace cfl;
using namespace std;

cfl::Function
prb::discountSvenssonFit (const std::vector<double> &rTimes,
                        const std::vector<double> &rDF,
                        double dLambda1, double dLambda2,
                        double dInitialTime, cfl::Function &rErr,
                        cfl::FitParam &rParam)
    {
   
    Function basis1 ([dLambda1, dInitialTime] (double dT) {
        if (dT==dInitialTime) {return 1.0;};
        return (1.-exp(-dLambda1*(dT-dInitialTime)))/(dLambda1*(dT-dInitialTime));
    });

    Function basis2 ([dLambda1, dInitialTime] (double dT) {
        if (dT==dInitialTime) {return 1.0;};
        return (1.-exp(-dLambda1*(dT-dInitialTime)))/(dLambda1*(dT-dInitialTime))-exp(-dLambda1*(dT-dInitialTime));
    });

    Function basis3 ([dLambda2, dInitialTime] (double dT) {
        if (dT==dInitialTime) {return 1.0;};
        return (1.-exp(-dLambda2*(dT-dInitialTime)))/(dLambda2*(dT-dInitialTime))-exp(-dLambda2*(dT-dInitialTime));
    });

    Fit rFit = prb::linear(vector<Function> {Function(1.0), basis1, basis2, basis3});
    Function ftr = prb::discountYieldFit(rTimes, rDF, dInitialTime, rFit, rErr);
    
    rParam = rFit.param();

    return ftr;
    }