#include "home2/home2.hpp"

using namespace cfl;
using namespace std;
namespace prb {
    double 
    yieldTMCouponBond(double dRate, double dPeriod, double dMaturity,
                            double dInitialTime, double dPrice, double dYield0,
                            double dErr) 
    {
        Function uF ([dRate, dPeriod, dMaturity, dInitialTime, dPrice] (double dYTM) {
            return couponBond(dRate, dPeriod, dMaturity, dYTM, dInitialTime, false)-dPrice;
        });
        Function uDF ([dRate, dPeriod, dMaturity, dInitialTime] (double dYTM) {
            return -durationCouponBond (dRate, dPeriod, dMaturity, dYTM, dInitialTime);
        });
        RootD uRoot = NRootD::newton (dErr);
        return yieldTMCouponBond(dRate, dPeriod, dMaturity, dInitialTime, dPrice, dYield0, uRoot);
    }

    double
    yieldTMCouponBond(double dRate, double dPeriod, double dMaturity,
                            double dInitialTime, double dPrice, double dYield0,
                            const cfl::RootD &rRootD) {
        
        Function uF ([dRate, dPeriod, dMaturity, dInitialTime, dPrice] (double dYTM) {
            return couponBond(dRate, dPeriod, dMaturity, dYTM, dInitialTime, false)-dPrice;
        });
        Function uDF ([dRate, dPeriod, dMaturity, dInitialTime] (double dYTM) {
            return -durationCouponBond (dRate, dPeriod, dMaturity, dYTM, dInitialTime);
        });

        return rRootD.find (uF, uDF, dYield0);
    }
}