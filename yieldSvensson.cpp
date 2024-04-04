#include "home1/home1.hpp"

using namespace cfl;
using namespace std;

cfl::Function
prb::yieldSvensson(double dC0, double dC1, double dC2, double dC3,
                             double dLambda1, double dLambda2,
                             double dInitialTime)
{
  // PRECONDITION ((dLambda >= 0) && (dSigma >= 0));

  std::function<double (double)> uY = [dC0, dC1, dC2, dC3, dLambda1, dLambda2, dInitialTime] (double dT) {
      double diffT = dT - dInitialTime;
      if (diffT == 0)
      {
        return dC0 + dC1;
      }
      return dC0 + dC1*(1.-exp(-dLambda1*diffT))/(dLambda1*diffT) + dC2*((1.-exp(-dLambda1*diffT))/(dLambda1*diffT)-exp(-dLambda1*diffT)) + dC3*((1.-exp(-dLambda2*diffT))/(dLambda2*diffT)-exp(-dLambda2*diffT));
    };
  return Function (uY, dInitialTime);
}

