#include "home3/Output.hpp"
#include "home3/home3.hpp"
#include "test/Data.hpp"
#include "test/Main.hpp"
#include "test/Print.hpp"

using namespace test;
using namespace cfl;
using namespace std;
using namespace test::Data;

void
testInterp (cfl::Interp &rInterp, const std::string &sInterp)
{
  print (sInterp.c_str ());
  std::string sInterpNodes ("Comparison at the nodes of interpolation");
  std::string sRandNodes ("Comparison at the random arguments");
  double dL = 0;
  double dR = 1;
  unsigned iN = 8;
  std::valarray<double> uX = getArg (dL, dR, iN);
  std::valarray<double> uY = getRandArg (dL, dR, iN);
  Function uArg ([] (double dX) { return dX; });
  double dC = 1;
  double dT = 2;
  double dP = 3;
  std::vector<std::string> uNames
      = { "interp", "deriv", "deriv2", "err", "err_d", "err_d2" };
  std::vector<Function> uF (uNames.size ());

  std::vector<unsigned> uColumns (uNames.size (), 10);
  unsigned iColumn = 12;
  unsigned iArg = 10;
  unsigned iSpace = 8;

  print (dC, "Interpolation of const function: f(x)");
  std::valarray<double> uV (dC, uX.size ());
  rInterp.assign (begin (uX), end (uX), begin (uV));
  uF[0] = rInterp.interp ();
  uF[1] = rInterp.deriv ();
  uF[2] = rInterp.deriv2 ();
  uF[3] = cfl::abs (uF[0] - dC);
  uF[4] = cfl::abs (uF[1]);
  uF[5] = cfl::abs (uF[2]);
  printTable (uF, uNames, uX, sInterpNodes, iColumn, iArg, iSpace);
  printTable (uF, uNames, uY, sRandNodes, iColumn, iArg, iSpace);

  cout << "Interpolation  of linear function: f(x) = " << dC << " + " << dT
       << "x" << endl
       << endl;
  uV = dC + dT * uX;
  rInterp.assign (begin (uX), end (uX), begin (uV));
  uF[0] = rInterp.interp ();
  uF[1] = rInterp.deriv ();
  uF[2] = rInterp.deriv2 ();
  uF[3] = cfl::abs (uF[0] - (dC + dT * uArg));
  uF[4] = cfl::abs (uF[1] - dT);
  uF[5] = cfl::abs (uF[2]);
  printTable (uF, uNames, uX, sInterpNodes, iColumn, iArg, iSpace);
  printTable (uF, uNames, uY, sRandNodes, iColumn, iArg, iSpace);

  cout << "Interpolation  of quadratic function: f(x) = " << dC << " + " << dT
       << "x + " << dP << "x^2" << endl
       << endl;
  uV = dC + dT * uX + dP * uX * uX;
  rInterp.assign (begin (uX), end (uX), begin (uV));
  uF[0] = rInterp.interp ();
  uF[1] = rInterp.deriv ();
  uF[2] = rInterp.deriv2 ();
  uF[3] = cfl::abs (uF[0] - (dC + dT * uArg + dP * uArg * uArg));
  uF[4] = cfl::abs (uF[1] - (dT + 2. * dP * uArg));
  uF[5] = cfl::abs (uF[2] - 2. * dP);
  printTable (uF, uNames, uX, sInterpNodes, iColumn, iArg, iSpace);
  printTable (uF, uNames, uY, sRandNodes, iColumn, iArg, iSpace);

  cout << "Interpolation  of cubic function: f(x) = " << dC << " + 0.5x^2 + "
       << dT << "x^3" << endl
       << endl;
  uV = dC + 0.5 * uX * uX + dT * uX * uX * uX;
  rInterp.assign (begin (uX), end (uX), begin (uV));
  uF[0] = rInterp.interp ();
  uF[1] = rInterp.deriv ();
  uF[2] = rInterp.deriv2 ();
  uF[3]
      = cfl::abs (uF[0] - (dC + 0.5 * uArg * uArg + dT * uArg * uArg * uArg));
  uF[4] = cfl::abs (uF[1] - (uArg + 3. * dT * uArg * uArg));
  uF[5] = cfl::abs (uF[2] - (1. + 6. * dT * uArg));
  printTable (uF, uNames, uX, sInterpNodes, iColumn, iArg, iSpace);
  printTable (uF, uNames, uY, sRandNodes, iColumn, iArg, iSpace);

  print ("Interpolation of exponential function: f(x) = exp(x)");
  uV = exp (uX);
  rInterp.assign (begin (uX), end (uX), begin (uV));
  uF[0] = rInterp.interp ();
  uF[1] = rInterp.deriv ();
  uF[2] = rInterp.deriv2 ();
  uF[3] = cfl::abs (uF[0] - exp (uArg));
  uF[4] = cfl::abs (uF[1] - exp (uArg));
  uF[5] = cfl::abs (uF[2] - exp (uArg));
  printTable (uF, uNames, uX, sInterpNodes, iColumn, iArg, iSpace);
  printTable (uF, uNames, uY, sRandNodes, iColumn, iArg, iSpace);

  print ("Interpolation of itself. Should get zero errors everywhere.");
  uV = getValues (uF.front (), uX);
  auto uG (uF);
  rInterp.assign (begin (uX), end (uX), begin (uV));
  uF[0] = rInterp.interp ();
  uF[1] = rInterp.deriv ();
  uF[2] = rInterp.deriv2 ();
  uF[3] = cfl::abs (uF[0] - uG[0]);
  uF[4] = cfl::abs (uF[1] - uG[1]);
  uF[5] = cfl::abs (uF[2] - uG[2]);
  printTable (uF, uNames, uX, sInterpNodes, iColumn, iArg, iSpace);
  printTable (uF, uNames, uY, sRandNodes, iColumn, iArg, iSpace);
};

void
cspline ()
{
  Interp uSpline = prb::cspline ();
  testInterp (uSpline, "CUBIC SPLINE INTERPOLATION WITH GSL");
}

void
akima ()
{
  Interp uSpline = prb::akima ();
  testInterp (uSpline, "AKIMA INTERPOLATION WITH GSL");
}

void
steffen ()
{
  Interp uSpline = prb::steffen ();
  testInterp (uSpline, "STEFFEN INTERPOLATION WITH GSL");
}

void
forwardCarryInterp ()
{
  test::print ("FORWARD PRICES BY INTERPOLATION OF COST-OF-CARRY RATES");

  double dSpot = 100;
  double dInitialTime = 1.;

  auto uF = test::Data::getForward (dSpot, dInitialTime);

  double dInitialCarryRate = std::log (uF.second.front () / dSpot)
                             / (uF.first.front () - dInitialTime);
  print (dInitialCarryRate, "initial carry rate", true);
  Interp uInterp = cfl::NInterp::linear ();
  print ("We use linear interpolation");

  Function uResult = prb::forwardCarryInterp (
      dSpot, uF.first, uF.second, dInitialCarryRate, dInitialTime, uInterp);

  double dInterval = uF.first.back () - dInitialTime;
  test::Data::print (uResult, dInitialTime, dInterval);
}

void
forwardCarrySteffenInterp ()
{
  test::print (
      "FORWARD PRICES  BY STEFFEN INTERPOLATION OF COST-OF-CARRY RATES");

  double dSpot = 100;
  double dInitialTime = 1.;
  auto uF = test::Data::getForward (dSpot, dInitialTime);

  Function uResult = prb::forwardCarrySteffenInterp (dSpot, uF.first,
                                                     uF.second, dInitialTime);

  double dInterval = uF.first.back () - dInitialTime;
  test::Data::print (uResult, dInitialTime, dInterval);
}

void
volatilityVarInterp ()
{
  test::print ("VOLATILITY CURVE BY INTERPOLATION OF VARIANCE CURVE");

  double dInitialTime = 1.;
  auto uV = test::Data::getVol (dInitialTime);
  Interp uInterp = NInterp::linear ();
  print ("We use linear interpolation");

  Function uResult
      = prb::volatilityVarInterp (uV.first, uV.second, dInitialTime, uInterp);

  double dInterval = uV.first.back () - dInitialTime;
  test::Data::print (uResult, dInitialTime, dInterval);
}

std::function<void ()>
test_home3 ()
{
  return [] () {
    print ("ONE-DIMENSIONAL INTERPOLATION");

    cspline ();
    steffen ();
    akima ();
    forwardCarryInterp ();
    forwardCarrySteffenInterp ();
    volatilityVarInterp ();
  };
}

int
main ()
{
  project (test_home3 (), PROJECT_NAME, PROJECT_NAME, "Homework 3");
}
