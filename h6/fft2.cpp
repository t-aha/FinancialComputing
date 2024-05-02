#include "home6/home6.hpp"
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include "cfl/GaussRollback.hpp"

using namespace cfl;
using namespace std;

class FFT2 : public IGaussRollback {
public:
    FFT2() {}

    FFT2(const unsigned iSize, const double dH, const double dVar): 
    FFT2() {

        int M = ceil(dVar / (2 * dH * dH));
        dF = M*2;
    }

    void rollback (std::valarray<double> &rValues) const {
        size_t n = rValues.size();
        rValues *= n;
        double *coeffs = new double[n];
        copy(begin(rValues), end(rValues), coeffs);

        gsl_fft_real_workspace *w = gsl_fft_real_workspace_alloc(n);
        gsl_fft_real_wavetable *wt = gsl_fft_real_wavetable_alloc(n);
        gsl_fft_real_transform(coeffs, 1, n, wt, w);

        for (size_t i = 1; i < (n+1)/2; ++i) {
            double tmp = exp((-2*M_PI*M_PI*i*i*dF)/(n*n));
            coeffs[2*i] *= tmp;
            coeffs[2*i-1] *= tmp;
        }
        if (n % 2 == 0)
        {
            int k = n/2;
            coeffs[n-1] *= exp((-2*M_PI*M_PI*k*k*dF)/(n*n));
        }

        gsl_fft_halfcomplex_wavetable *hc = gsl_fft_halfcomplex_wavetable_alloc(n);
        gsl_fft_halfcomplex_inverse(coeffs, 1, n, hc, w);

        for (size_t j = 0; j < n; ++j)
        {
            rValues[j] = coeffs[j]/n;
        }

        gsl_fft_real_wavetable_free (wt);
        gsl_fft_halfcomplex_wavetable_free (hc);
        gsl_fft_real_workspace_free(w);
        delete[] coeffs;
    }

    IGaussRollback *
    newObject(const unsigned iSize, const double dH, const double dVar)  const {
        return new FFT2(iSize, dH, dVar);
    }

private:
    double dF;
    int M;
    double q;
};

cfl::GaussRollback prb::fft2() {
    return GaussRollback(new FFT2());
}

