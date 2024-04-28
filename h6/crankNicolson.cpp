#include "home6/home6.hpp"
#include <gsl/gsl_linalg.h>
#include "cfl/GaussRollback.hpp"

using namespace cfl;
using namespace std;

class CrankNicolson : public IGaussRollback {
public:
    CrankNicolson(const double dR) : r (dR) {}

    CrankNicolson(const double dR, const unsigned iSize, const double dH, const double dVar): 
    CrankNicolson(dR) {

        int M = ceil(dVar / (dH * dR));
        double q = dVar / (2 * dH * dH * M);
        
        this->M = M;
        this->q = q;
    }

    void rollback (std::valarray<double> &rValues) const {
        size_t N = rValues.size();
        for (int m = 0; m < M; ++m) {
            // Define the tridiagonal matrix coefficients
            gsl_vector *diag = gsl_vector_alloc(N);
            gsl_vector *sup = gsl_vector_alloc(N - 1);
            gsl_vector *sub = gsl_vector_alloc(N - 1);

            // Define the right-hand side vector
            gsl_vector *b = gsl_vector_alloc(N);

            // Initialize the tridiagonal matrix coefficients and the right-hand side vector
            for (size_t i = 0; i < N; ++i) {
                gsl_vector_set(diag, i, 1 + q);
            }
            for (size_t i = 0; i < N - 1; ++i) {
                gsl_vector_set(sup, i, -0.5*q);
                gsl_vector_set(sub, i, -0.5*q);
            }

            for (size_t i = 1; i < N - 1; ++i) {
                gsl_vector_set(b, i, 
                q/2*rValues[i+1] + (1-q)*rValues[i] + q/2 * rValues[i-1]);
            }
            gsl_vector_set(b, 0, 0);
            gsl_vector_set(b, N-1, 0);

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
        return new CrankNicolson(r, iSize, dH, dVar);
    }

private:
    double r;
    int M;
    double q;
};

cfl::GaussRollback prb::crankNicolson(double dR) {
    return GaussRollback(new CrankNicolson(dR));
}

