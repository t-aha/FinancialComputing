#include "home1/home1.hpp"

using namespace cfl;
using namespace std;

cfl::Function
prb::carryBlack(double dTheta, double dLambda, double dSigma,
                              double dInitialTime)
{
  PRECONDITION ((dLambda >= 0) && (dSigma >= 0));

  std::function<double (double)> uC = [dTheta, dLambda, dSigma, dInitialTime] (double dT) {
      double diffT = dT - dInitialTime;
      if (diffT == 0) {
        return dTheta+dSigma*dSigma*0.5;
      }
      return dTheta*(1-exp(-dLambda*diffT))/(dLambda*diffT) + 0.5*(dSigma*dSigma)*(1.-exp(-2.*dLambda*diffT))/(2.*dLambda*diffT);
    };
  return Function (uC, dInitialTime);
}

