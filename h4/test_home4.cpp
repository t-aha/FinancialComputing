#include "home4/Output.hpp"
#include "home4/home4.hpp"
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
multiLinear ()
{
  print ("MULTI-DIM LINEAR LEAST-SQUARES FITTING");

  double dL = -1.;
  double dR = 2.;
  unsigned iPoints = 10;
  std::valarray<double> uFitNodes = getArg (dL, dR, iPoints);
  std::valarray<double> uRandNodes = getRandArg (dL, dR, iPoints);
  unsigned iOrder = 3; // iPoints / 3;
  Function uLinFunc ([] (double dX) { return dX; });

  std::vector<Function> uBase (iOrder);
  uBase.front () = Function (1.);
  for (unsigned i = 0; i + 1 < iOrder; i++)
    {
      uBase[i + 1] = uBase[i] * uLinFunc;
    }
  print (iOrder - 1, "Fitting with polynomials of order");
  print ("free function = exp(x)");
  Function uFree ([] (double dX) { return exp (dX); });

  Fit uFit = prb::linear (uBase, uFree);
  testFit (uFit, uFitNodes, uRandNodes);
}

const std::string c_sDF ("Fitted discount factors and their errors:");
const std::string c_sConstYield ("We fit with constant yield.");

void
discountYieldFit ()
{
  test::print ("DISCOUNT CURVE OBTAINED BY LEAST-SQUARES FIT OF YIELD CURVE");

  double dInitialTime = 1.;

  auto uDF = test::Data::getDiscount (dInitialTime);

  print (c_sConstYield);
  cfl::Fit uFit = cfl::NFit::linear (cfl::Function (1.));

  Function uErr;
  Function uDiscount = prb::discountYieldFit (uDF.first, uDF.second,
                                              dInitialTime, uFit, uErr);

  double dInterval = uDF.first.back () - dInitialTime;
  test::Data::printFit (uFit.param ());
  test::Data::print (c_sDF, uDiscount, uErr, dInitialTime, dInterval);
}
void
discountHullWhiteFit ()
{
  test::print ("LEAST-SQUARES FIT OF DISCOUNT CURVE FOR HULL AND WHITE MODEL");

  double dLambda = 0.05;
  double dInitialTime = 1.;

  print (dLambda, "lambda");
  auto uDF = test::Data::getDiscount (dInitialTime);

  Function uErr;
  FitParam uParam;
  Function uDiscount = prb::discountHullWhiteFit (
      uDF.first, uDF.second, dLambda, dInitialTime, uErr, uParam);
  double dInterval = uDF.first.back () - dInitialTime;
  test::Data::printFit (uParam);
  test::Data::print (c_sDF, uDiscount, uErr, dInitialTime, dInterval);
}

void
discountNelsonSiegelFit ()
{
  test::print ("LEAST-SQUARES FIT OF DISCOUNT CURVE IN NELSON-SIEGEL MODEL");

  double dLambda = 0.05;
  double dInitialTime = 1.;

  print (dLambda, "lambda");
  auto uDF = test::Data::getDiscount (dInitialTime);

  Function uErr;
  FitParam uParam;
  Function uDiscount = prb::discountNelsonSiegelFit (
      uDF.first, uDF.second, dLambda, dInitialTime, uErr, uParam);
  double dInterval = uDF.first.back () - dInitialTime;
  test::Data::printFit (uParam);
  test::Data::print (c_sDF, uDiscount, uErr, dInitialTime, dInterval);
}

void
discountSvenssonFit ()
{
  test::print ("LEAST-SQUARES FIT OF DISCOUNT CURVE IN SVENSSON MODEL");

  double dLambda1 = 0.02;
  double dLambda2 = 0.1;
  double dInitialTime = 1.;

  print (dLambda1, "lambda1");
  print (dLambda2, "lambda2");
  auto uDF = test::Data::getDiscount (dInitialTime);

  Function uErr;
  FitParam uParam;
  Function uDiscount = prb::discountSvenssonFit (
      uDF.first, uDF.second, dLambda1, dLambda2, dInitialTime, uErr, uParam);
  double dInterval = uDF.first.back () - dInitialTime;
  test::Data::printFit (uParam);
  test::Data::print (c_sDF, uDiscount, uErr, dInitialTime, dInterval);
}

std::function<void ()>
test_home4 ()
{
  return [] () {
    print ("LEAST-SQUARES FITTING");

    multiLinear ();
    discountYieldFit ();
    discountHullWhiteFit ();
    discountNelsonSiegelFit ();
    discountSvenssonFit ();
  };
}

int
main ()
{
  project (test_home4 (), PROJECT_NAME, PROJECT_NAME, "Homework 4");
}
