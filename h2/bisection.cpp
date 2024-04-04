#include "home2/home2.hpp"
#include "home2/Output.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

struct tfunction
{
  cfl::Function rF;
};
double
value (double dX, void *param)
{
  return ((tfunction *)param)->rF (dX);
}

class AbsRoot : public cfl::IRoot
{
public:
  AbsRoot (double dFuncErr, unsigned iMaxSteps,
           const gsl_root_fsolver_type *type)
      : dFuncErr (dFuncErr), iMaxSteps (iMaxSteps),
        type (type)
  {
  }

    double
    find(const cfl::Function &rF, double dL, double dR) const
    {
        const gsl_root_fsolver_type *T = this->type;
        gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
        gsl_function gf;
        cfl::Function ug(rF);
        gf.function = &value;
        gf.params = &ug;
        gsl_root_fsolver_set(s, &gf, dL, dR);

        int iStatus = GSL_CONTINUE;
        double r;
        unsigned M = 0;
        while (iStatus == GSL_CONTINUE)
        {
            if (M++ >= this->iMaxSteps)
                break;
            gsl_root_fsolver_iterate(s);
            r = gsl_root_fsolver_root(s);
            double residual = std::abs(rF(r));
            iStatus = gsl_root_test_residual(residual, this->dFuncErr);
        }
        gsl_root_fsolver_free(s);
        return r;
    }

private:
  double dFuncErr, iMaxSteps;
  const gsl_root_fsolver_type *type;
};

cfl::Root
prb::bisection (double dFuncErr, unsigned iMaxSteps)
{
  return cfl::Root (

      new AbsRoot (dFuncErr, iMaxSteps, gsl_root_fsolver_bisection));
}

class RelRoot : public cfl::IRoot
{
public:
  RelRoot (double dAbsErr, double dRelErr, unsigned iMaxSteps,
           const gsl_root_fsolver_type *type)
      : dAbsErr (dAbsErr), dRelErr (dRelErr), iMaxSteps (iMaxSteps),
        type (type)
  {
  }

  double
  find (const cfl::Function &rF, double dL, double dR) const
  {
    const gsl_root_fsolver_type *T = this->type;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc (T);
    gsl_function gf;
    cfl::Function ug (rF);
    gf.function = &value;
    gf.params = &ug;
    gsl_root_fsolver_set (s, &gf, dL, dR);

    int iStatus = GSL_CONTINUE;
    double r, x_lo, x_hi;
    unsigned M = 0;
    while (iStatus == GSL_CONTINUE)
      {
        if (M++ >= this->iMaxSteps)
          break;
        gsl_root_fsolver_iterate (s);
        r = gsl_root_fsolver_root (s);
        x_lo = gsl_root_fsolver_x_lower (s);
        x_hi = gsl_root_fsolver_x_upper (s);
        iStatus = gsl_root_test_interval (x_lo, x_hi, this->dAbsErr,
                                          this->dRelErr);
      }
    gsl_root_fsolver_free (s);
    return r;
  }

private:
  double dAbsErr, dRelErr, iMaxSteps;
  const gsl_root_fsolver_type *type;
};

cfl::Root
prb::bisection (double dAbsErr, double dRelErr, unsigned iMaxSteps)
{
  return cfl::Root (

      new RelRoot (dAbsErr, dRelErr, iMaxSteps, gsl_root_fsolver_bisection));
}