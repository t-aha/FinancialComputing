#include "home1/home1.hpp"

using namespace cfl;
using namespace std;

cfl::Function
prb::forwardRate(double dPeriod, const cfl::Function &rDiscount)
{
    std::function<double (double)> uF
        = [dPeriod, rDiscount] (double dT) {
        return (rDiscount(dT)/rDiscount(dT+dPeriod)-1)/dPeriod;
    };
    return Function (uF);
}
