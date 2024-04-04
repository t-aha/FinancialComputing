#include "home2/home2.hpp"

using namespace cfl;
using namespace std;

double
prb::durationCouponBond(double dRate, double dPeriod, double dMaturity,
                            double dYTM, double dInitialTime) 
{
    double del = 1e-6;
    double tmp1 = prb::couponBond(dRate, dPeriod, dMaturity,dYTM+del,dInitialTime,true);
    double tmp2 = prb::couponBond(dRate, dPeriod, dMaturity,dYTM,dInitialTime,true);
    return -(tmp1-tmp2)/del;
}