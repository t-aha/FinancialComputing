#include "home6/home6.hpp"
#include <gsl/gsl_matrix.h>
#include "cfl/GaussRollback.hpp"

using namespace cfl;
using namespace std;

class Explicit : public IGaussRollback {
public:
    Explicit(const double dP) : p (dP) {}

    Explicit(const double dP, const unsigned iSize, const double dH, const double dVar): 
    Explicit(dP) {

        int M = ceil(dVar / (2 * dH * dH * dP));
        double q = dVar / (2 * dH * dH * M);
        this->m_dH = dH;
        this->m_dVar = dVar;
        this->N = iSize;
        this->M = M;
        this->q = q;
    }

    void rollback (std::valarray<double> &rValues) const {
        std::valarray<double> delta(N);

        for (int m = 0; m < M; ++m) {
            for (unsigned int n = 1; n < N-1; ++n) {
                delta[n] = rValues[n-1] - 2*rValues[n]+ rValues[n+1];
            }
            delta[0] = delta[1];
            delta[N-1] = delta[N-2];
            rValues = rValues + q*delta;
        }
        
    }

    IGaussRollback *
    newObject(const unsigned iSize, const double dH, const double dVar)  const {
        return new Explicit(p, iSize, dH, dVar);
    }

private: 
    double m_dH, m_dVar;
    unsigned int N;
    double p, q;
    int M;
};

cfl::GaussRollback prb::expl(double dP) {
    return GaussRollback(new Explicit(dP));
}

