#include "home6/home6.hpp"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "cfl/GaussRollback.hpp"

using namespace cfl;
using namespace std;

class Chain : public IGaussRollback {
public:
    Chain(const unsigned iExplSteps, 
            const cfl::GaussRollback &rFast,
            unsigned iImplSteps, double dExplP,
            double dImplP) 
        : n_expl (iExplSteps), fast (rFast), n_impl (iImplSteps),
          explP (dExplP), implP (dImplP) {}

    Chain(const unsigned iExplSteps, 
            const cfl::GaussRollback &rFast,
            unsigned iImplSteps, double dExplP,
            double dImplP, const unsigned iSize, 
            const double dH, const double dVar): 
    Chain(iExplSteps, rFast, iImplSteps, dExplP, dImplP) {
        Te = 2*dH*dH*dExplP*iExplSteps;
        Ti = 2*dH*dH*dImplP*iImplSteps;
        flag = dVar > (Te + Ti);
        size = iSize;
        var = dVar;
        h = dH;
    }

    void rollback (std::valarray<double> &rValues) const {
        if (flag) {
            // start w explicit scheme
            rollback_expl(rValues, explP, n_expl, h, Te);
            // continue with fast scheme
            GaussRollback tmp = fast;
            tmp.assign(size, h, var - (Te + Ti));
            tmp.rollback(rValues);
            // finish with implicit scheme
            rollback_impl(rValues, implP, n_impl, h, Ti);

        } else {
            rollback_expl(rValues, explP, size, h, var);
        }
    }

    void rollback_expl (std::valarray<double> &rValues, double dP, unsigned imsize, double dH, double dVar) const {
        std::valarray<double> delta(imsize);
        int M = ceil(dVar / (2 * dH * dH * dP));
        double q = dVar / (2 * dH * dH * M);
        for (int m = 0; m < M; ++m) {
            for (unsigned int n = 1; n < imsize-1; ++n) {
                delta[n] = rValues[n-1] - 2*rValues[n]+ rValues[n+1];
            }
            delta[0] = delta[1];
            delta[imsize-1] = delta[imsize-2];
            rValues = rValues + q*delta;
        }
    
    }

    void rollback_impl (std::valarray<double> &rValues, double dP, unsigned imsize, double dH, double dVar) const {
        int M = ceil(dVar / (2 * dH * dH * dP));
        size_t N = imsize;
        for (int m = 0; m < M; ++m) {
            // Define the tridiagonal matrix coefficients
            gsl_vector *diag = gsl_vector_alloc(N);
            gsl_vector *sup = gsl_vector_alloc(N - 1);
            gsl_vector *sub = gsl_vector_alloc(N - 1);

            // Define the right-hand side vector
            gsl_vector *b = gsl_vector_alloc(N);

            // Initialize the tridiagonal matrix coefficients and the right-hand side vector
            for (size_t i = 0; i < N; ++i) {
                gsl_vector_set(diag, i, 1 + 2 * dP);
                gsl_vector_set(b, i, rValues[i]);
            }
            for (size_t i = 0; i < N - 1; ++i) {
                gsl_vector_set(sup, i, -dP);
                gsl_vector_set(sub, i, -dP);
            }

            // Allocate memory for the solution vector
            gsl_vector *x = gsl_vector_alloc(N);

            // Solve the tridiagonal linear system
            gsl_linalg_solve_tridiag(diag, sup, sub, b, x);

            // Update rValues with the solution
            for (size_t i = 0; i < N; ++i) {
                rValues[i] = gsl_vector_get(x, i);
            }

            // Free allocated memory
            gsl_vector_free(diag);
            gsl_vector_free(sup);
            gsl_vector_free(sub);
            gsl_vector_free(b);
            gsl_vector_free(x);
        }
    }

    IGaussRollback *
    newObject(const unsigned iSize, const double dH, const double dVar)  const {
        return new Chain(n_expl, fast, n_impl, explP, implP, iSize, dH, dVar);
    }

private:
    bool flag;
    unsigned n_expl, n_impl, size;
    cfl::GaussRollback fast;
    double explP, implP, Te, Ti, var, h;
};

cfl::GaussRollback prb::chain(unsigned iExplSteps, const cfl::GaussRollback &rFast,
                          unsigned iImplSteps, double dExplP, double dImplP) {
    return GaussRollback(new Chain(iExplSteps, rFast, iImplSteps, dExplP, dImplP));
}

