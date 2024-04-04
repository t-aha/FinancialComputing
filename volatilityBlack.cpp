#include "home1/home1.hpp"

using namespace cfl;
using namespace std;

cfl::Function
prb::volatilityBlack(double dSigma, double dLambda,
                               double dInitialTime)
{
  // PRECONDITION ((dLambda >= 0) && (dSigma >= 0));

  std::function<double (double)> uV = [dSigma, dLambda, dInitialTime] (double dT) {
      double diffT = dT - dInitialTime;
      if (diffT == 0)
      {
        return dSigma;
      }
      return dSigma * sqrt((1.-exp(-2.*dLambda*diffT))/(2.*dLambda*diffT));
    };
  return Function (uV, dInitialTime);
}

