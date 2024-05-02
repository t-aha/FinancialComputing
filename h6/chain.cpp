#include "home6/home6.hpp"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "cfl/GaussRollback.hpp"

using namespace cfl;
using namespace std;

class ChainFit : public IGaussRollback
{
public:
  ChainFit (unsigned iExplSteps, const cfl::GaussRollback &rFast,
               unsigned iImplSteps, double dExplP, double dImplP)
      : n_expl (iExplSteps), n_impl (iImplSteps), fast (rFast),
        explP (dExplP), implP (dImplP)
  {
  }

  ChainFit (const unsigned iSize, const double dh, const double var,
               unsigned iExplSteps, const cfl::GaussRollback &rFast,
               unsigned iImplSteps, double dExplP, double dImplP)
      : ChainFit (iExplSteps, rFast, iImplSteps, dExplP, dImplP)
  {

    m_var = var;
    h = dh;
    Te = 2 * dh * dh * dExplP * iExplSteps;
    Ti = 2 * dh * dh * dImplP * iImplSteps;
    flag = (var > Te + Ti);
    if (flag){fast.assign (iSize, dh, m_var - (Te + Ti));}
  }

  ChainFit *
  newObject (const unsigned iSize, const double dh, const double var) const
  {
    return new ChainFit (iSize, dh, var, n_expl, fast,
                            n_impl, explP, implP);
  }

  void
  rollback_expl (std::valarray<double> &rValues) const
  {
    double q = Te / (2.0 * h * h * n_expl);
    size_t N = rValues.size ();
    std::valarray<double> rDelta (N);
    for (size_t m = 0; m < n_expl; ++m)
      {
        for (size_t n = 1; n < N - 1; ++n)
          {
            rDelta[n] = rValues[n - 1] - 2 * rValues[n] + rValues[n + 1];
          }
        rDelta[0] = rDelta[1];
        rDelta[N - 1] = rDelta[N - 2];

        rValues += q * rDelta;
      }
  }

  void
  rollback_impl (std::valarray<double> &rValues) const
  {
    double q = Ti / (2.0 * h * h * n_impl);
    size_t N = rValues.size ();
    gsl_vector *diag = gsl_vector_alloc (N);
    gsl_vector *sub = gsl_vector_alloc (N - 1);

    for (size_t i = 0; i < N; ++i)
      {
        gsl_vector_set (diag, i, 1 + 2 * q);
        if (i < N - 1)
          {
            gsl_vector_set (sub, i, -q);
          }
      }

    gsl_vector *b = gsl_vector_alloc(N);
    gsl_vector *x = gsl_vector_alloc(N);
    gsl_vector_view rV_view = gsl_vector_view_array (&rValues[0], N);
    gsl_vector_memcpy(b, &rV_view.vector);

    for (size_t j = 0; j < n_impl; ++j)
      {
        gsl_linalg_solve_tridiag (diag, sub, sub, b, x);
        gsl_vector_memcpy (b, x);
      }

    gsl_vector_memcpy (&rV_view.vector, x);

    gsl_vector_free (b);
    gsl_vector_free (x);
    gsl_vector_free (diag);
    gsl_vector_free (sub);
  }

  void
  rollback (std::valarray<double> &rValues) const
  {
    if (m_var > Te + Ti)
      {
        rollback_expl (rValues);
        fast.rollback (rValues);
        rollback_impl (rValues);
      }
    else
      {
        rollback_expl (rValues);
      }
  }

private:
  bool flag;
  unsigned n_expl, n_impl;
  cfl::GaussRollback fast;
  double explP, implP, Te, Ti, m_var, h;
};

cfl::GaussRollback
prb::chain (unsigned iExplSteps, const cfl::GaussRollback &rFast,
            unsigned iImplSteps, double dExplP, double dImplP)
{
  return GaussRollback (
      new ChainFit (iExplSteps, rFast, iImplSteps, dExplP, dImplP));
}