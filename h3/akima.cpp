#include "cfl/Interp.hpp"
#include "cfl/Error.hpp"
#include <gsl/gsl_spline.h>
#include "home3/home3.hpp"

using namespace cfl;
using namespace std;

class AkimaInterp : public IInterp
{
public:
  AkimaInterp() {}

  AkimaInterp(const std::vector<double> &rArg, const std::vector<double> &rVal)
      : m_dL(rArg.front()), m_dR(rArg.back())
  {
    PRECONDITION((rArg.size() == rVal.size()) && (rArg.size() >= 2));
    PRECONDITION(
        is_sorted(rArg.begin(), rArg.end(), std::less_equal<double>()));

    m_uS.reset(gsl_spline_alloc(gsl_interp_akima, rArg.size()),
               &gsl_spline_free);
    // copies of rArg and rVal will be created
    gsl_spline_init(m_uS.get(), rArg.data(), rVal.data(), rArg.size());
  }

  IInterp *
  newObject(const std::vector<double> &rArg,
            const std::vector<double> &rVal) const
  {
    return new AkimaInterp(rArg, rVal);
  }

  Function
  interp() const
  {
    std::shared_ptr<gsl_interp_accel> uAcc(gsl_interp_accel_alloc(),
                                           &gsl_interp_accel_free);

    std::function<double(double)> uF = [uS = m_uS, uAcc](double dX) {
      return gsl_spline_eval(uS.get(), dX, uAcc.get());
    };

    return Function(uF, m_dL, m_dR);
  }

  Function
  deriv() const
  {
    std::shared_ptr<gsl_interp_accel> uAcc(gsl_interp_accel_alloc(),
                                           &gsl_interp_accel_free);

    std::function<double(double)> uF = [uS = m_uS, uAcc](double dX) {
      return gsl_spline_eval_deriv(uS.get(), dX, uAcc.get());
    };

    return Function(uF, m_dL, m_dR);
  }

  Function
  deriv2() const
  {
    std::shared_ptr<gsl_interp_accel> uAcc(gsl_interp_accel_alloc(),
                                           &gsl_interp_accel_free);

    std::function<double(double)> uF = [uS = m_uS, uAcc](double dX) {
      return gsl_spline_eval_deriv2(uS.get(), dX, uAcc.get());
    };

    return Function(uF, m_dL, m_dR);
  }

private:
  std::shared_ptr<gsl_spline> m_uS;
  double m_dL, m_dR;
};

// functions from NInterp

cfl::Interp
prb::akima()
{
  return Interp(new AkimaInterp());
}