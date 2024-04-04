#include "home2/home2.hpp"

using namespace cfl;
using namespace std;

cfl::Function
prb::forwardYTMCashFlow(const std::vector<double> &rPayments,
                            const std::vector<double> &rPaymentTimes,
                            const cfl::Function &rForward,
                            double dInitialTime, double dY0, double dY1,
                            const cfl::Root &rRoot) 
{
    std::function<double(double)> uF = [rPayments, rPaymentTimes, rForward, dInitialTime, dY0, dY1, rRoot](double dT) {
        double val = rForward(dT);
        vector<double> fPayments;
        vector<double> fPaymentTimes;
        for (int i = 0; i < int(rPayments.size()); i++) {
            if (rPaymentTimes[i] > dT) {
                fPayments.push_back(rPayments[i]);
                fPaymentTimes.push_back(rPaymentTimes[i]);
            }
        }
        double g0 = yieldTMCashFlow(fPayments, fPaymentTimes, val, dT, dY0, dY1, rRoot);
        return g0;
    };
    return Function(uF, dInitialTime);
}