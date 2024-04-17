#include "home5/home5.hpp"
#include "cfl/Fit.hpp"
#include "cfl/Error.hpp"
#include <algorithm>
#include <cmath>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <numeric>

using namespace cfl;
using namespace std;

class BSpline : public IFit {
    public:
        BSpline(unsigned iOrder, const std::vector<double> &rBreakpoints)
            : order(iOrder), bp(rBreakpoints) {}

        BSpline(unsigned iOrder, const std::vector<double> &rBreakpoints,
            const std::vector<double> &rArg, const std::vector<double> &rVal,
            const std::vector<double> &rWt, bool bChi2)
            : BSpline(iOrder, rBreakpoints) {
            size_t n = rArg.size();
            size_t nbreak = bp.size();
            size_t ncoeffs = nbreak + order - 2;

            gsl_bspline_workspace *bw = gsl_bspline_alloc(order, nbreak);
            gsl_vector *B = gsl_vector_alloc(ncoeffs);

            gsl_vector *x = gsl_vector_alloc(n);
            gsl_vector *y = gsl_vector_alloc(n);
            gsl_matrix *X = gsl_matrix_alloc(n, ncoeffs);
            gsl_vector *c = gsl_vector_alloc(ncoeffs);
            gsl_vector *w = gsl_vector_alloc(n);
            gsl_matrix *cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
            gsl_multifit_linear_workspace *mw = gsl_multifit_linear_alloc(n, ncoeffs);
            gsl_vector *brkpnt = gsl_vector_alloc(bp.size());

            for (size_t i = 0; i < n; ++i) {
                gsl_vector_set(x, i, rArg[i]);
                gsl_vector_set(y, i, rVal[i]);
                gsl_vector_set(w, i, rWt[i]);
                if (i >= bp.size()) {continue;};
                gsl_vector_set(brkpnt, i, bp[i]);
            }

            gsl_bspline_knots(brkpnt, bw);

            for (size_t i = 0; i < n; ++i) {
                double xi = gsl_vector_get(x, i);
                gsl_bspline_eval(xi, B, bw);
                for (size_t j = 0; j < ncoeffs; ++j) {
                    double Bj = gsl_vector_get(B, j);
                    gsl_matrix_set(X, i, j, Bj);
                }
            }

            gsl_multifit_wlinear(X, w, y, c, cov, &m_dChi2, mw);

            #define C(i) (gsl_vector_get(c, (i)))
            #define COV(i,j) (gsl_matrix_get(cov,(i),(j)))
            
            for (size_t b = 0; b < ncoeffs; b++) {
                m_dC.push_back(C(b));
                for (size_t q = 0; q < ncoeffs; q++) {
                    if (bChi2) {
                        double dVar = m_dChi2 / (n - ncoeffs);
                        m_dCov.push_back(COV(b,q)*dVar);
                    } else {
                    m_dCov.push_back(COV(b,q));
                    }
                }
            }
    

            gsl_bspline_free(bw);
            gsl_vector_free(B);
            gsl_vector_free(x);
            gsl_vector_free(y);
            gsl_matrix_free(X);
            gsl_vector_free(c);
            gsl_vector_free(w);
            gsl_matrix_free(cov);
            gsl_multifit_linear_free(mw);
        
    }

    IFit *newObject(const std::vector<double> &rArg, const std::vector<double> &rVal,
                    const std::vector<double> &rWt, bool bChi2) const {
        return new BSpline(order, bp, rArg, rVal, rWt, bChi2);
    }

    Function fit() const {
        //return m_uParam.fit[0] + m_uParam.fit[1] * m_uBaseF + m_uFreeF;
        
        Function fittedFunction ([*this](double x) {
            double result = 0.0;
            size_t ncoeffs = m_dC.size();
            gsl_vector *B = gsl_vector_alloc(ncoeffs);
            gsl_bspline_workspace *bw = gsl_bspline_alloc(order, bp.size());
            gsl_vector *brkpnt = gsl_vector_alloc(bp.size());
            
            for (size_t i = 0; i < bp.size(); ++i) {
                gsl_vector_set(brkpnt, i, bp[i]);
            }
            
            gsl_bspline_knots(brkpnt, bw);
            gsl_bspline_eval(x, B, bw);

            for (size_t i = 0; i < ncoeffs; ++i) {
                double Bi = gsl_vector_get(B, i);
                result += m_dC[i] * Bi;
            }

            gsl_vector_free(B);
            return result;
        });
        return fittedFunction;
    }

    Function err() const {
        Function errorFunction ([*this](double x) {
            size_t ncoeffs = m_dC.size();
            gsl_vector *c = gsl_vector_alloc(ncoeffs);
            gsl_matrix *cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
            gsl_vector *B = gsl_vector_alloc(ncoeffs);
            gsl_bspline_workspace *bw = gsl_bspline_alloc(order, bp.size());
            gsl_vector *brkpnt = gsl_vector_alloc(bp.size());
            
            for (size_t i = 0; i < bp.size(); ++i) {
                gsl_vector_set(c, i, m_dC[i]);
                gsl_vector_set(brkpnt, i, bp[i]);
                for (size_t j = 0; j < bp.size(); ++j){
                    gsl_matrix_set(cov, i, j, m_dCov[i*ncoeffs+j]);
                }
            }
            double yi, yerr;
            gsl_bspline_knots(brkpnt, bw);
            gsl_bspline_eval(x, B, bw);
            gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
            return yerr;
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
    unsigned order;
    std::vector<double> bp;

    std::vector<Function> bsfun;
    std::vector<double> m_dC;
    std::vector<double> m_dCov; 
    double m_dChi2;
};

// linear one-dim

cfl::Fit prb::bspline(unsigned iOrder, const std::vector<double> &rBreakpoints) {
    return Fit(new BSpline(iOrder, rBreakpoints));
}

cfl::Fit prb::bspline (unsigned iOrder, double dL, double dR,
                  unsigned iBreakpoints)
{
    std::vector<double> breakpoints;
    breakpoints.reserve(iBreakpoints);
    for (unsigned i = 0; i < iBreakpoints; ++i) {
        double value = dL + i*(dR - dL) / (iBreakpoints-1);
        breakpoints.push_back(value);
    }
    return Fit(new BSpline(iOrder, breakpoints));
}