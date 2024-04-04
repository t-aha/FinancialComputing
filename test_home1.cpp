#include "home1/Output.hpp"
#include "home1/home1.hpp"
#include "test/Data.hpp"
#include "test/Main.hpp"
#include "test/Parameters.hpp"
#include "test/Print.hpp"

using namespace test;
using namespace cfl;
using namespace std;
using namespace test::Data;

void
carryBlack ()
{
  test::print ("COST-OF-CARRY RATE IN BLACK MODEL");

  double dTheta = 0.03;
  double dLambda = 0.05;
  double dSigma = 0.2;
  double dInitialTime = 0.75;

  print (dTheta, "theta");
  print (dLambda, "lambda");
  print (dSigma, "sigma");
  print (dInitialTime, "initial time", true);

  Function uCostOfCarry
      = prb::carryBlack (dTheta, dLambda, dSigma, dInitialTime);
  double dInterval = 1;
  test::Data::print (uCostOfCarry, dInitialTime, dInterval);
}

void
yieldSvensson ()
{
  test::print ("SVENSSON YIELD CURVE");

  double dLambda1 = 0.05;
  double dLambda2 = 0.07;
  double dC0 = 0.02;
  double dC1 = 0.04;
  double dC2 = 0.06;
  double dC3 = 0.03;
  double dInitialTime = 1.5;

  print (dC0, "c0");
  print (dC1, "c1");
  print (dC2, "c2");
  print (dC3, "c3");
  print (dLambda1, "lambda 1");
  print (dLambda2, "lambda 2");
  print (dInitialTime, "initial time", true);

  Function uYield = prb::yieldSvensson (dC0, dC1, dC2, dC3, dLambda1, dLambda2,
                                        dInitialTime);
  double dInterval = 5;
  test::Data::print (uYield, dInitialTime, dInterval);
}

void
volatilityBlack ()
{
  test::print ("STATIONARY IMPLIED VOLATILITY IN BLACK MODEL");
  double dLambda = 0.05;
  double dSigma = 0.2;
  double dInitialTime = 0.75;

  print (dSigma, "sigma");
  print (dLambda, "lambda");
  print (dInitialTime, "initial time", true);

  Function uVol = prb::volatilityBlack (dSigma, dLambda, dInitialTime);
  double dInterval = 1;
  test::Data::print (uVol, dInitialTime, dInterval);
}

void
forwardStockDividends ()
{
  test::print ("FORWARD PRICES FOR A STOCK WITH DIVIDENDS");

  double dSpot = 100;
  double dRate = 0.12;
  double dInitialTime = 1.;
  Function uDiscount = cfl::Data::discount (dRate, dInitialTime);
  unsigned iTimes = 10;
  std::vector<double> uTimes (iTimes);
  double dPeriod = 0.5;
  uTimes.front () = dInitialTime + dPeriod;
  std::transform (uTimes.begin (), uTimes.end () - 1, uTimes.begin () + 1,
                  [dPeriod] (double dX) { return dX + dPeriod; });
  double dDividend = 5.;
  double dStep = 0.5;
  std::vector<double> uDividends (iTimes, dDividend);
  std::transform (uDividends.begin (), uDividends.end () - 1,
                  uDividends.begin () + 1,
                  [dStep] (double dD) { return dD + dStep; });

  print (dInitialTime, "initial time");
  print (dRate, "interest rate");
  print (dSpot, "spot", true);
  test::Data::print ("Stock dividends:", uTimes, uDividends);

  Function uForwardStockDividends = prb::forwardStockDividends (
      dSpot, uTimes, uDividends, uDiscount, dInitialTime);
  double dInterval = iTimes * dPeriod / 1.01;
  test::Data::print (uForwardStockDividends, dInitialTime, dInterval);
}

void
forwardCouponBond ()
{
  test::print ("FORWARD PRICES FOR A COUPON BOND");

  cfl::Data::CashFlow uBond = test::swapParameters ();
  uBond.notional = 1.;
  double dRate = uBond.rate;
  double dInitialTime = 1.;
  Function uDiscount = cfl::Data::discount (dRate, dInitialTime);
  test::print (dRate, "interest rate");
  test::print (dInitialTime, "initial time", true);

  test::printCashFlow (uBond, "bond parameters");

  for (int iI = 0; iI < 2; iI++)
    {
      bool bClean = (iI == 0) ? true : false;
      if (bClean)
        {
          print ("clean prices:");
        }
      else
        {
          print ("dirty prices:");
        }
      double dRate = uBond.rate;
      double dPeriod = uBond.period;
      double dMaturity = dInitialTime + dPeriod * uBond.numberOfPayments;
      Function uForwardCouponBond = prb::forwardCouponBond (
          dRate, dPeriod, dMaturity, uDiscount, dInitialTime, bClean);
      double dInterval = uBond.period * uBond.numberOfPayments / 1.1;
      test::Data::print (uForwardCouponBond, dInitialTime, dInterval);
    }
}

void
forwardRate ()
{
  test::print ("FORWARD RATES");

  double dRate = 0.03;
  double dInitialTime = 1.5;
  double dPeriod = 0.25;

  print (dPeriod, "rate period");
  print (dInitialTime, "initial time", true);

  Function uDiscount = cfl::Data::discount (dRate, dInitialTime);
  uDiscount += cfl::Data::discount (dRate * 4, dInitialTime);
  uDiscount *= 0.5;

  Function uForwardRate = prb::forwardRate (dPeriod, uDiscount);

  double dInterval = 5;
  test::Data::print (uForwardRate, dInitialTime, dInterval);
}

std::function<void ()>
test_home1 ()
{
  return [] () {
    print ("DATA CURVES FOR FINANCIAL MODELS");

    carryBlack ();
    yieldSvensson ();
    volatilityBlack ();
    forwardStockDividends ();
    forwardCouponBond ();
    forwardRate ();
  };
}

int
main ()
{
  project (test_home1 (), PROJECT_NAME, PROJECT_NAME, "Homework 1");
}
