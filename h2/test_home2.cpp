#include "home2/Output.hpp"
#include "home2/home2.hpp"
#include "test/Data.hpp"
#include "test/Main.hpp"
#include "test/Parameters.hpp"
#include "test/Print.hpp"
#include <numeric>

using namespace test;
using namespace cfl;
using namespace std;
using namespace test::Data;

void
couponBond ()
{
  test::print ("PRICE FOR A COUPON BOND");

  cfl::Data::CashFlow uBond = test::swapParameters ();
  uBond.notional = 1.;
  double dYTM = uBond.rate * 1.1;
  double dIssueTime = 1.;
  double dInitialTime = dIssueTime + 0.5 * uBond.period;
  double dMaturity = dIssueTime + uBond.period * uBond.numberOfPayments;
  Function uDiscount = cfl::Data::discount (dYTM, dInitialTime);
  test::print (dYTM, "yield to maturity", true);

  test::printCashFlow (uBond, "bond parameters");

  test::print (dIssueTime, "issue time");
  test::print (dInitialTime, "initial time");
  test::print (dMaturity, "maturity", true);

  bool bClean = false;
  double dPD = prb::couponBond (uBond.rate, uBond.period, dMaturity, dYTM,
                                dInitialTime, bClean);
  test::print (dPD, "dirty price");

  bClean = true;
  double dPC = prb::couponBond (uBond.rate, uBond.period, dMaturity, dYTM,
                                dInitialTime, bClean);
  test::print (dPC, "clean price");

  print (dPD - dPC, "accrued interest", true);
}

void
durationCouponBond ()
{
  test::print ("DURATION FOR A COUPON BOND");

  cfl::Data::CashFlow uBond = test::swapParameters ();
  uBond.notional = 1.;
  double dYTM = uBond.rate * 1.1;
  double dIssueTime = 1.;
  double dInitialTime = dIssueTime + 0.5 * uBond.period;
  double dMaturity = dIssueTime + uBond.period * uBond.numberOfPayments;
  Function uDiscount = cfl::Data::discount (dYTM, dInitialTime);
  test::print (dYTM, "yield to maturity", true);

  test::printCashFlow (uBond, "bond parameters");

  test::print (dIssueTime, "issue time");
  test::print (dInitialTime, "initial time");
  test::print (dMaturity, "maturity", true);

  double dD = prb::durationCouponBond (uBond.rate, uBond.period, dMaturity,
                                       dYTM, dInitialTime);
  test::print (dD, "dollar duration", true);
}

void
yieldTMCouponBond_simple ()
{
  test::print ("YIELD TO MATURITY FOR COUPON BOND");

  cfl::Data::CashFlow uBond = test::swapParameters ();
  uBond.notional = 1.;
  double dIssueTime = 1.;
  double dInitialTime = dIssueTime + 0.5 * uBond.period;
  double dMaturity = dIssueTime + uBond.period * uBond.numberOfPayments;
  double dPrice = 1. + (dInitialTime - dIssueTime) * uBond.rate;
  double dY0 = uBond.rate / 2;
  double dErr = cfl::EPS;

  test::printCashFlow (uBond, "bond parameters");

  print (dMaturity, "maturity");
  test::print (dPrice, "price");
  test::print (dInitialTime, "initial time");
  test::print (dY0, "initial guess for YTM");
  print (dErr, "expected yield error", true);

  double dYTM = prb::yieldTMCouponBond (uBond.rate, uBond.period, dMaturity,
                                        dInitialTime, dPrice, dY0, dErr);
  test::print (dYTM, "YTM", true);
}

void
yieldTMCouponBond ()
{
  test::print ("YIELD TO MATURITY FOR COUPON BOND WITH cfl::RootD");

  cfl::Data::CashFlow uBond = test::swapParameters ();
  uBond.notional = 1.;
  double dIssueTime = 1.;
  double dInitialTime = dIssueTime + 0.5 * uBond.period;
  double dMaturity = dIssueTime + uBond.period * uBond.numberOfPayments;
  double dPrice = 1. + (dInitialTime - dIssueTime) * uBond.rate;
  double dY0 = uBond.rate / 2;
  double dErrR = 10E-6;
  double dErrA = 10E-8;

  test::printCashFlow (uBond, "bond parameters");

  print (dMaturity, "maturity");
  test::print (dPrice, "price");
  test::print (dInitialTime, "initial time");
  test::print (dY0, "initial guess for YTM");
  print (dErrR, "relative yield error");
  print (dErrA, "absolute yield error");
  print ("method = secant");

  RootD uRoot = cfl::NRootD::secant (dErrR, dErrA);

  double dYTM = prb::yieldTMCouponBond (uBond.rate, uBond.period, dMaturity,
                                        dInitialTime, dPrice, dY0, uRoot);
  test::print (dYTM, "YTM", true);
}

void
yieldTMCashFlow ()
{
  test::print ("YIELD TO MATURITY FOR A CASH FLOW");

  double dInitialTime = 1.;
  std::vector<double> uPayments = { 1, -2, 3, 5, 6, -2 };
  std::vector<double> uTimes
      = getTimes (dInitialTime + 0.5, dInitialTime + 3., uPayments.size ());

  double dY = 0.05;
  double dPrice = 0.;
  dPrice = std::inner_product (
      uPayments.begin (), uPayments.end (), uTimes.begin (), dPrice,
      std::plus<double> (),
      [dInitialTime, dY] (double dPayment, double dTime) {
        return dPayment * std::exp (-dY * (dTime - dInitialTime));
      });

  test::print (dPrice, "initial price");
  test::print (dInitialTime, "initial time", true);
  test::Data::print ("cash flow:", uTimes, uPayments);
  double dY0 = 0.0;
  double dY1 = 0.2;
  test::print (dY0, "lower bound");
  test::print (dY1, "upper bound");
  double dEpsR = 10E-6;
  double dEpsA = 10E-8;
  test::print (dEpsR, "relative error for YTM");
  test::print (dEpsA, "absolute error for YTM");
  test::print ("method = brent");

  cfl::Root uRoot = prb::brent (dEpsR, dEpsA);

  double dYTM = prb::yieldTMCashFlow (uPayments, uTimes, dPrice, dInitialTime,
                                      dY0, dY1, uRoot);
  test::print (dYTM, "YTM");
  test::print (std::abs (dYTM - dY), "error", true);
}

void
forwardYTMCashFlow ()
{
  test::print ("FORWARD YIELD TO MATURITY FOR A CASH FLOW");

  double dInitialTime = 1.;
  std::vector<double> uPayments = { 1, -2, -3, 5, 6, 2 };
  double dYield = 0.05;

  cfl::Function uYield (
      [dYield] (double dT) { return dYield * (1. + 0.1 * dT); }, dInitialTime);
  cfl::Function uDiscount = cfl::Data::discount (uYield, dInitialTime);

  std::vector<double> uTimes
      = getTimes (dInitialTime + 0.5, dInitialTime + 3., uPayments.size ());

  std::function<double (double)> uF = [&uDiscount, &uPayments,
                                       &uTimes] (double dT) {
    unsigned iTime = std::lower_bound (uTimes.begin (), uTimes.end (), dT)
                     - uTimes.begin ();
    double dPrice = 0.;
    dPrice = std::inner_product (uPayments.begin () + iTime, uPayments.end (),
                                 uTimes.begin () + iTime, dPrice,
                                 std::plus<double> (),
                                 [&uDiscount] (double dPayment, double dTime) {
                                   return dPayment * uDiscount (dTime);
                                 });
    dPrice /= uDiscount (dT);

    return dPrice;
  };

  cfl::Function uForward (uF, dInitialTime);

  test::print (dInitialTime, "initial time", true);
  test::Data::print ("cash flow:", uTimes, uPayments);
  double dY0 = 0;
  double dY1 = 0.2;
  test::print (dY0, "lower bound");
  test::print (dY1, "upper bound");
  double dEpsR = 10E-6;
  double dEpsA = 10E-8;
  test::print (dEpsR, "relative error for YTM");
  test::print (dEpsA, "absolute error for YTM");
  test::print ("method = brent");

  cfl::Root uRoot = prb::brent (dEpsR, dEpsA);

  Function uForwardYTMCashFlow = prb::forwardYTMCashFlow (
      uPayments, uTimes, uForward, dInitialTime, dY0, dY1, uRoot);
  double dInterval = (uTimes.back () - dInitialTime) / 1.01;
  test::Data::print (uForwardYTMCashFlow, dInitialTime, dInterval);
}

void
testRoot (const cfl::Root &rRoot, const cfl::Function &rF, double dL,
          double dR, double dRoot)
{
  print (dL, "initial left bound");
  print (dR, "initial right bound");
  double dX = rRoot.find (rF, dL, dR);
  print (dX, "root");
  print (std::abs (dX - dRoot), "absolute error for root");
  print (std::abs (rF (dX) - rF (dRoot)), "absolute error for value at root",
         true);
}

void
testRoot (const cfl::Root &rRoot, const std::string &sRoot)
{
  print (sRoot.c_str ());

  double dL = 0.25;
  double dR = 0.75;

  print ("Linear function: f(x) = 2x -1");
  Function uF ([] (double dX) { return 2 * dX - 1; });
  double dRoot = 0.5;
  testRoot (rRoot, uF, dL, dR, dRoot);

  print ("Quadtratic function: f(x) = x^2 - 0.25");
  uF = Function ([] (double dX) { return dX * dX - 0.25; });
  testRoot (rRoot, uF, dL, dR, dRoot);

  print ("Exponential function: f(x) = e^x - e^{0.5}");
  uF = Function ([] (double dX) { return std::exp (dX) - std::exp (0.5); });
  testRoot (rRoot, uF, dL, dR, dRoot);
};

std::function<void ()>
test_home2 ()
{
  return [] () {
    print ("ONE-DIMENSIONAL ROOT FINDING");

    couponBond ();
    durationCouponBond ();
    yieldTMCouponBond_simple ();
    yieldTMCouponBond ();

    double dAbsErr = cfl::EPS;
    double dRelErr = cfl::EPS;
    double dAbsFErr = cfl::EPS;
    unsigned iSteps = 3; // cfl::IMAX;

    for (unsigned i = 0; i < 2; i++)
      {
        print ("TEST OF ROOT BRACKETING ALGORITHMS");
        print (dAbsErr, "absolute error");
        print (dRelErr, "relative error");
        print (dAbsFErr, "residual error");
        print (iSteps, "maximal number of steps", true);

        testRoot (prb::bisection (dAbsErr, dRelErr, iSteps),
                  "BISECTION WITH INTERVAL ERRORS");
        testRoot (prb::bisection (dAbsFErr, iSteps),
                  "BISECTION WITH RESIDUAL ERRORS");
        testRoot (prb::falsepos (dAbsErr, dRelErr, iSteps),
                  "FALSE POSITION WITH INTERVAL ERRORS");
        testRoot (prb::falsepos (dAbsFErr, iSteps),
                  "FALSE POSITION WITH RESIDUAL ERRORS");
        testRoot (prb::brent (dAbsErr, dRelErr, iSteps),
                  "BRENT WITH INTERVAL ERRORS");
        testRoot (prb::brent (dAbsFErr, iSteps), "BRENT WITH RESIDUAL ERRORS");
        iSteps = 200;
      }

    yieldTMCashFlow ();
    forwardYTMCashFlow ();
  };
}

int
main ()
{
  project (test_home2 (), PROJECT_NAME, PROJECT_NAME, "Homework 2");
}
