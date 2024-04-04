#include "home1/home1.hpp"

using namespace cfl;
using namespace std;

cfl::Function
prb::forwardStockDividends(double dSpot,
                       const std::vector<double> &rDividendsTimes,
                       const std::vector<double> &rDividends,
                       const cfl::Function &rDiscount, double dInitialTime)
{
    std::function<double (double)> uF
        = [dSpot, rDividendsTimes, rDividends, rDiscount, dInitialTime] (double dT) {
        double dR = -log(rDiscount(dInitialTime+1));
        double out = dSpot*exp(dR*(dT-dInitialTime));
        int i = 0;
        while (dT > rDividendsTimes[i]) {
            out -= rDividends[i]*exp(dR*(dT - rDividendsTimes[i]));
            i++;
        }
        return out;
    };
    return Function (uF, dInitialTime);
}
