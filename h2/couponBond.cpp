#include "home2/home2.hpp"

using namespace cfl;
using namespace std;

double
prb::couponBond(double dRate, double dPeriod, double dMaturity, double dYTM,
                   double dInitialTime, bool bClean)
{
    double dSum = 0.0;
    double totT = (dMaturity - dInitialTime)/dPeriod;
    int numPayments = floor(totT);
    double aI = totT - numPayments;
    int i = 1;
    for (i = 1; i <= numPayments+1; i++) {
        dSum += exp(-dYTM*(i*dPeriod+1-dInitialTime));
    }
    dSum *= dRate*dPeriod;
    dSum += exp(-dYTM*(dMaturity-dInitialTime));
    if (!bClean) 
    {
        return dSum;
    }

    return dSum - dRate*(dPeriod*aI);
}
