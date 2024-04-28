#include "home6/home6.hpp"
#include <gsl/gsl_matrix.h>
#include "cfl/GaussRollback.hpp"

using namespace cfl;
using namespace std;

class FFT2 : public IGaussRollback {
public:
    FFT2() {}

    FFT2(const unsigned iSize, const double dH, const double dVar): 
    FFT2() {

        int M = static_cast<int>(dVar / (2 * dH * dH));
        double q = dVar / (2 * dH * dH * M);
        
        this->M = M;
        this->q = q;
    }

    void rollback (std::valarray<double> &rValues) const {
        int N = rValues.size();
        std::valarray<double> delta(N);
        
        for (int n = 1; n < N - 1; ++n) {
            delta[n] = rValues[n - 1] - 2 * rValues[n] + rValues[n + 1];
        }

        delta[0] = delta[1];
        delta[N - 1] = delta[N - 2];

        for (int n = 0; n < N; ++n) {
            rValues[n] += q * delta[n];
        }
    }

    IGaussRollback *
    newObject(const unsigned iSize, const double dH, const double dVar)  const {
        return new FFT2(iSize, dH, dVar);
    }

private:
    double p;
    int M;
    double q;
};

cfl::GaussRollback prb::fft2() {
    return GaussRollback(new FFT2());
}

