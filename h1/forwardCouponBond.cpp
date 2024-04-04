#include "home1/home1.hpp"

using namespace cfl;
using namespace std;

cfl::Function
prb::forwardCouponBond(double dRate, double dPeriod, double dMaturity,
                       const cfl::Function &rDiscount, double dInitialTime, bool bClean)
{
    std::function<double(double)> forwardPrice = [dRate, dPeriod, dMaturity, rDiscount, bClean](double dT) {
        double dPayTime = dMaturity;
        double dSum = 0.0;
        
        while (dPayTime > dT) {
            dSum += rDiscount(dPayTime);
            dPayTime -= dPeriod;
        }
        
        double dPayment = dRate * dPeriod;
        dSum *= dPayment;
        dSum += rDiscount(dMaturity);
        double dF = dSum / rDiscount(dT);
        
        if (bClean) {
            dF -= dRate * (dT - dPayTime);
        }
        
        return dF;
    };
    return Function(forwardPrice);
}