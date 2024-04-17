#include "home4/home4.hpp"

#include "cfl/Fit.hpp"
#include "cfl/Error.hpp"
#include <algorithm>
#include <cmath>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_multifit.h>
#include <numeric>

using namespace cfl;
using namespace std;


class MultiDimFit : public IFit {
public:
    MultiDimFit(const std::vector<cfl::Function> &rBasisF,
                const cfl::Function &rFreeF)
        : m_uBasisF(rBasisF), m_uFreeF(rFreeF) {}

    MultiDimFit(const std::vector<cfl::Function> &rBasisF,
                const Function rFreeF, const std::vector<double> &rArg,
                const std::vector<double> &rVal, const std::vector<double> &rWt,
                bool bChi2)
        : MultiDimFit(rBasisF, rFreeF) {
        size_t n = rArg.size();
        size_t p = rBasisF.size(); // number of basis functions

        gsl_matrix *X = gsl_matrix_alloc(n, p);
        gsl_vector *y = gsl_vector_alloc(n);
        gsl_vector *w = gsl_vector_alloc(n);

        // Populate X matrix with basis functions evaluated at data points
        for (size_t i = 0; i < n; ++i) {
            gsl_vector_set(y, i, rVal[i]-rFreeF(rArg[i]));
            gsl_vector_set(w, i, rWt[i]);
            for (size_t j = 0; j < p; ++j) {
                gsl_matrix_set(X, i, j, rBasisF[j](rArg[i]));
            }
        }

        gsl_vector *c = gsl_vector_alloc(p);
        gsl_matrix *cov = gsl_matrix_alloc(p, p);
        double chisq;

        gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(n, p);
        gsl_multifit_wlinear(X, w, y, c, cov, &chisq, work);

        gsl_multifit_linear_free(work);
        gsl_matrix_free(X);
        gsl_vector_free(y);
        gsl_vector_free(w);

        #define C(i) (gsl_vector_get(c, (i)))
        #define COV(i,j) (gsl_matrix_get(cov,(i),(j)))
        m_dChi2 = chisq;
        
        for (size_t b = 0; b < p; b++) {
          m_dC.push_back(C(b));
          for (size_t q = 0; q < p; q++) {
            if (bChi2) {
                double dVar = m_dChi2 / (n - p);
                m_dCov.push_back(COV(b,q)*dVar);
            } else {
              m_dCov.push_back(COV(b,q));
            }
          }
        }
        gsl_vector_free(c);
        gsl_matrix_free(cov);
    }

    IFit *newObject(const std::vector<double> &rArg, const std::vector<double> &rVal,
                    const std::vector<double> &rWt, bool bChi2) const {
        return new MultiDimFit(m_uBasisF, m_uFreeF, rArg, rVal, rWt, bChi2);
    }

    Function fit() const {
        //return m_uParam.fit[0] + m_uParam.fit[1] * m_uBaseF + m_uFreeF;
        Function fittedFunction ([*this](double x) {
            double result = 0.0;
            for (size_t i = 0; i < m_dC.size(); i++) {
                result += m_dC[i] * m_uBasisF[i](x); // Add contribution from each basis function
            }
            return result + m_uFreeF(x); // Add the contribution from the free function
        });
        return fittedFunction;
    }

    Function err() const {
        Function errorFunction ([*this](double x) {
            double result = 0.0;
            for (size_t i = 0; i < m_dC.size(); ++i) {
                for (size_t j = 0; j < m_dC.size(); ++j) {
                    result +=  m_uBasisF[i](x) * m_dCov[i*(m_dC.size())+j] * m_uBasisF[j](x); 
                }
            }
            return sqrt(result);
        });
        return errorFunction;
    }

    FitParam param() const {
        FitParam uParam;
        uParam.fit = std::valarray<double> (m_dC.data(), m_dC.size());
        uParam.cov = std::valarray<double> (m_dCov.data(), m_dCov.size());
        uParam.chi2 = m_dChi2;

        return uParam;
    }

private:
    std::vector<Function> m_uBasisF;
    Function m_uFreeF;
    std::vector<double> m_dC;
    std::vector<double> m_dCov; 
    double m_dChi2;
};

// linear one-dim

cfl::Fit prb::linear(const std::vector<cfl::Function> &rBasisF,
                     const cfl::Function &rFreeF) {
    return Fit(new MultiDimFit(rBasisF, rFreeF));
}
