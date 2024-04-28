#include "home6/home6.hpp"
#include "cfl/GaussRollback.hpp"
#include <gsl/gsl_linalg.h>

using namespace cfl;
using namespace std;

class Implicit : public IGaussRollback {
public:
    Implicit(const double dP) : p (dP) {}

    Implicit(const double dP, const unsigned iSize, const double dH, const double dVar): 
    Implicit(dP) {

        int M = ceil(dVar / (2 * dH * dH * dP));
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
                gsl_vector_set(diag, i, 1 + 2 * p);
                gsl_vector_set(b, i, rValues[i]);
            }
            for (size_t i = 0; i < N - 1; ++i) {
                gsl_vector_set(sup, i, -p);
                gsl_vector_set(sub, i, -p);
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
        return new Implicit(p, iSize, dH, dVar);
    }

private:
    unsigned int N;
    double p, q;
    int M;
};

cfl::GaussRollback prb::impl(double dP) {
    return GaussRollback(new Implicit(dP));
}

