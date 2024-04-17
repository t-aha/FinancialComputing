#include "home5/Output.hpp"
#include "home5/home5.hpp"
#include "test/Data.hpp"
#include "test/Main.hpp"
#include "test/Print.hpp"

using namespace test;
using namespace cfl;
using namespace std;
using namespace test::Data;

void
testFit (Fit &rFit, const std::valarray<double> &rFitNodes,
         const std::valarray<double> &rRandNodes)
{
  std::vector<string> uNames = { "exact", "numeric", "est err", "act err" };
  Function uFunc ([] (double dX) { return 2 * std::exp (dX); });

  print ("We fit the function f(x) = 2*exp(x)");
  std::valarray<double> uFitValues = getValues (uFunc, rFitNodes);
  rFit.assign (begin (rFitNodes), end (rFitNodes), begin (uFitValues));

  test::Data::printFit (rFit.param ());

  std::vector<Function> uF
      = { uFunc, rFit.fit (), rFit.err (), abs (uFunc - rFit.fit ()) };

  unsigned iColumn = 12;
  unsigned iArg = 10;
  unsigned iSpace = 8;

  test::printTable (uF, uNames, rFitNodes, "Fit at the nodes of fitting",
                    iColumn, iArg, iSpace);
  printChi2 (rFit.err (), uF[uF.size () - 1], rFitNodes);

  test::printTable (uF, uNames, rRandNodes, "Fit at random nodes", iColumn,
                    iArg, iSpace);
  printChi2 (rFit.err (), uF[uF.size () - 1], rRandNodes);
}

void
bspline ()
{
  print ("LEAST-SQUARES FITTING WITH BASIS SPLINES");
  double dL = -1.;
  double dR = 2.;
  unsigned iBreakpoints = 4;
  unsigned iOrder = 2;
  unsigned iPoints = iBreakpoints + iOrder + 4;
  std::valarray<double> uFitNodes = getArg (dL, dR, iPoints);
  std::valarray<double> uRandNodes = getRandArg (dL, dR, iPoints);

  print ("BASIS SPLINES WITH GENERAL BREAKPOINTS");

  std::valarray<double> uBrP = getRandArg (dL, dR, iBreakpoints - 2);
  std::vector<double> uBreakpoints (iBreakpoints);
  uBreakpoints.front () = dL;
  uBreakpoints.back () = dR;
  std::copy (begin (uBrP), end (uBrP), uBreakpoints.begin () + 1);
  print (iOrder, "order");
  print (uBreakpoints.begin (), uBreakpoints.end (), "breakpoints");
  Fit uFit = prb::bspline (iOrder, uBreakpoints);
  testFit (uFit, uFitNodes, uRandNodes);

  print ("BASIS SPLINES WITH UNIFORM BREAKPOINTS");

  print (iOrder, "order");
  cout << "interval = [" << dL << "," << dR << "]" << endl;
  print (iBreakpoints, "number of breakpoints", true);
  uFit = prb::bspline (iOrder, dL, dR, iBreakpoints);
  testFit (uFit, uFitNodes, uRandNodes);
}

const std::string
    c_sForwardFX ("Fitted forward exchange rates and their errors:");

void
forwardFXCarryFit ()
{
  test::print ("FORWARD PRICES FOR EXCHANGE RATES BY LEAST-SQUARES FIT");

  double dSpotFX = 100;
  double dInitialTime = 1.;

  auto uDF = test::Data::getDiscountFX (dSpotFX, dInitialTime);

  double dL = dInitialTime;
  double dR = get<0> (uDF).back ();
  unsigned iOrder = 2;
  unsigned iBreakpoints = 3;
  print ("Basis spline with uniform breakpoints");
  print (iOrder, "order");
  cout << "interval = [" << dL << "," << dR << "]" << endl;
  print (iBreakpoints, "number of breakpoints", true);
  Fit uFit = prb::bspline (iOrder, dL, dR, iBreakpoints);

  Function uErr;

  Function uForwardFX
      = prb::forwardFXCarryFit (dSpotFX, get<0> (uDF), get<1> (uDF),
                                get<2> (uDF), dInitialTime, uFit, uErr);

  double dInterval = get<0> (uDF).back () - dInitialTime;
  test::Data::printFit (uFit.param ());
  test::Data::print (c_sForwardFX, uForwardFX, uErr, dInitialTime, dInterval);
}

const std::string c_sDF ("Fitted discount factors and their errors:");

void
discountRateFit ()
{
  test::print ("LEAST-SQUARES FIT OF DISCOUNT CURVE TO FLOAT RATES");

  double dInitialTime = 1.;

  auto uL = test::Data::getRate (dInitialTime);

  double dL = dInitialTime;
  double dR = dInitialTime + uL.first.back ();
  unsigned iOrder = 2;
  unsigned iBreakpoints = 3;
  print ("Basis spline with uniform breakpoints");
  print (iOrder, "order");
  cout << "interval = [" << dL << "," << dR << "]" << endl;
  print (iBreakpoints, "number of breakpoints", true);
  Fit uFit = prb::bspline (iOrder, dL, dR, iBreakpoints);

  Function uErr;
  Function uDiscount
      = prb::discountRateFit (uL.first, uL.second, dInitialTime, uFit, uErr);

  double dInterval = uL.first.back () / 1.1;
  test::Data::printFit (uFit.param ());
  test::Data::print (c_sDF, uDiscount, uErr, dInitialTime, dInterval);
}

void
discountApproxFit ()
{
  test::print ("APPROXIMATE LEAST-SQUARES FIT OF DISCOUNT CURVE");

  double dInitialTime = 1.;
  auto uDF = test::Data::getDiscount (dInitialTime);

  double dL = dInitialTime;
  double dR = uDF.first.back ();
  unsigned iOrder = 2;
  unsigned iBreakpoints = 3;
  print ("Basis spline with uniform breakpoints");
  print (iOrder, "order");
  cout << "interval = [" << dL << "," << dR << "]" << endl;
  print (iBreakpoints, "number of breakpoints", true);
  Fit uFit = prb::bspline (iOrder, dL, dR, iBreakpoints);

  Function uErr;
  Function uDiscount = prb::discountApproxFit (uDF.first, uDF.second,
                                               dInitialTime, uFit, uErr);
  double dInterval = uDF.first.back () - dInitialTime;
  test::Data::printFit (uFit.param ());
  test::Data::print (c_sDF, uDiscount, uErr, dInitialTime, dInterval);
}

void
discountSwapFit1 ()
{
  test::print ("DISCOUNT CURVE OBTAINED FROM SWAP RATES BY LEAST-SQUARES FIT");

  double dPeriod = 1. / 12;
  double dInitialTime = 1.;

  std::vector<double> uSwapRates
      = test::Data::getSwapRates (dPeriod, dInitialTime);

  double dL = dInitialTime;
  double dR = dInitialTime + uSwapRates.size () * dPeriod;
  unsigned iOrder = 2;
  unsigned iBreakpoints = 3;
  print ("Basis spline with uniform breakpoints");
  print (iOrder, "order");
  cout << "interval = [" << dL << "," << dR << "]" << endl;
  print (iBreakpoints, "number of breakpoints", true);
  Fit uFit = prb::bspline (iOrder, dL, dR, iBreakpoints);

  Function uErr;
  Function uDiscount
      = prb::discountSwapFit (uSwapRates, dPeriod, dInitialTime, uFit, uErr);

  double dInterval = uSwapRates.size () * dPeriod / 1.1;
  test::Data::printFit (uFit.param ());
  test::Data::print (c_sDF, uDiscount, uErr, dInitialTime + 0.001, dInterval);
}

std::function<void ()>
test_home5 ()
{
  return [] () {
    print ("LEAST-SQUARES FITTING");

    bspline ();
    forwardFXCarryFit ();
    discountRateFit ();
    discountSwapFit1 ();
    discountApproxFit ();
  };
}

int
main ()
{
  project (test_home5 (), PROJECT_NAME, PROJECT_NAME, "Homework 5");
}
