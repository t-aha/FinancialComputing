#include "home2/home2.hpp"

using namespace cfl;
using namespace std;

double
prb::yieldTMCashFlow(const std::vector<double> &rPayments,
                         const std::vector<double> &rPaymentTimes,
                         double dValue, double dInitialTime, double dY0,
                         double dY1, const cfl::Root &rRoot) 
{
    Function uF ([rPayments, rPaymentTimes, dInitialTime, dValue] (double dYTM) {
        int n = rPayments.size();
        double outs = 0.0;
        for (int i = 0; i < n; i++) {
            outs += rPayments[i]*exp(-dYTM*(rPaymentTimes[i]-dInitialTime));
        }
        return outs-dValue;
    });
    return rRoot.find(uF, dY0, dY1);
}