#include "cfl/GaussRollback.hpp"
#include "home6/home6.hpp"
#include <cmath>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fft_real.h>
using namespace cfl;
using namespace std;

class FFT : public IGaussRollback
{
public:
  FFT () {}

  FFT (const unsigned iSize, const double dH, const double dVar)
      : FFT ()
  {
    dF = dVar / (dH * dH);
  }

  FFT *newObject (const unsigned iSize, const double dH, const double dVar) const
  {
    return new FFT(iSize, dH, dVar);
  }

  void rollback (std::valarray<double> &rValues) const {
    size_t n = rValues.size ();
    double *coeffs = new double[n];
    rValues *= n;   
    std::copy (std::begin (rValues), std::end (rValues), coeffs);

    gsl_fft_real_workspace *w = gsl_fft_real_workspace_alloc(n);
    gsl_fft_real_wavetable *wt = gsl_fft_real_wavetable_alloc(n);
    gsl_fft_real_transform(coeffs, 1, n, wt, w);

    for (size_t i = 1; i < (n+1)/2; ++i)
      {
        double tmp = exp((-2*M_PI*M_PI*i*i*dF) / (n*n));
        coeffs[2*i] *= tmp;
        coeffs[2*i-1] *= tmp;
      }
    if (n % 2 == 0)
      {
        int k = n/2;
        coeffs[n-1] *= exp (-2*M_PI*M_PI*k*k*dF) / (n*n);
      }

    gsl_fft_halfcomplex_wavetable *hc = gsl_fft_halfcomplex_wavetable_alloc (n);
    gsl_fft_halfcomplex_inverse(coeffs, 1, n, hc, w);

    for (size_t j = 0; j < n; ++j) 
      {
        rValues[j] = coeffs[j]/n;
      }

    gsl_fft_real_wavetable_free(wt);
    gsl_fft_halfcomplex_wavetable_free(hc);
    gsl_fft_real_workspace_free(w);
    delete[] coeffs;
  }

private:
  double dF;
};

cfl::GaussRollback
prb::fft ()
{
  return GaussRollback (new FFT ());
}