/*
 * lsqr.cc - LSQR implementation
 * 
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#include <cmath>
#include <algorithm>
#include <numeric>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "axb_solver.hpp"
#include "in_house_solver.hpp"
#include "error.hpp"

namespace {
    bool
    is_nan(double const& d) {
        return d != d;
    }
    
    template<typename T, typename A, typename Op>
    void
    omp_transform(std::vector<T,A>& v, Op op) {
        std::size_t const n = v.size();
        T* d = v.data();
#pragma omp parallel for 
        for(int i = 0; i < n; ++i) {
            d[i] = op(d[i]);
        }
    }
    
    template<class A>
    double
    normalize(std::vector<double,A>& x) {
        double norm=0;
        int const n = x.size();
        double* d = x.data();
#pragma omp parallel for reduction(+:norm)
        for (int i = 0; i < n; ++i) {
            assert(!is_nan(x[i]));
            norm += std::pow(d[i], 2);
        }
        norm = std::sqrt(norm);
        double ratio = 1/norm;
#pragma omp parallel for 
        for (int i = 0; i < n; ++i) {
            d[i] = ratio*d[i];
        }
        return norm;
    }

    template<class T>
    T const& cst(T& v) {
        return const_cast<T const&>(v);
    }
    
    class in_house_omp_solver 
        : public axb_solver {
    public:
        in_house_omp_solver() = default;
        virtual ~in_house_omp_solver() {}
        
        std::string name() const override { return "<in house OMP solver>"; }

        int
        solve_impl(SparseRectangular const& mat, std::vector<double> const& rshv,
                   std::vector<double>& xvar, double tolerance, int& max_iter) const override {
            // check inputs
            int nnode = mat.nCol();
            int ndata = mat.nRow();
            if (nnode != xvar.size() || ndata != rshv.size()){
                std::cerr << "iterativeSolver_LSQR::Size mismatch "
                          << nnode << " " << xvar.size() << ", "
                          << ndata << " " << rshv.size() << '\n';
                std::exit(1);
            }
            
            aligned_vector<double> v;
            aligned_vector<double> w;
            aligned_vector<double> u(rshv.begin(), rshv.end());
            
            // initialization
            double beta = normalize(u);
            v = mat.omp_mul_transposed(u);
            double alpha = normalize(v); 
            w = v;
            std::fill(xvar.begin(), xvar.end(),0.0);
            double phibar = beta;
            double rhobar = alpha;
            
            double bnorm = beta; // |b|
            double Bnorm = 0.0; // |B|
            
            int istop = 0;
            int iter = 0;
            
            while (istop == 0){
                ++iter;
                double nalpha = -alpha;
                omp_transform(u, [=](double const& d) { return d*nalpha; });
                mat.omp_mul(v, u, false);
                beta = normalize(u); 
                double nbeta = -beta;
                omp_transform(v, [=](double const& d) { return d*nbeta; });
                mat.omp_mul_transposed(cst(u),v,false);
                alpha = normalize(v);
                
                Bnorm += std::pow(alpha,2)+std::pow(beta,2); // |B|
                
                double inv_rho = 1/std::sqrt(std::pow(rhobar,2)+ std::pow(beta,2));
                double c = rhobar*inv_rho;
                double s = beta*inv_rho;
                double theta = s*alpha;
                rhobar = -c*alpha;
                double phi = c*phibar;
                phibar *= s;
                
                double phirho = phi*inv_rho; // update x and w
                double thetarho = theta*inv_rho;
#pragma omp parallel for 
                for (int i=0; i<nnode; ++i){
                    double tmp = w[i];
                    xvar[i] += phirho*tmp;
                    w[i] = v[i]-thetarho*tmp;
                }
                
                // check stoping criteria
                double Anorm = sqrt(Bnorm); // |a|
                // c id a ratio of two sqrt, how could it be < 0 ???
                double conv = (alpha*abs(c))/Anorm; // |Atr|/(|A|*|r|)
                if (is_nan(conv)) {
                    std::cout << std::scientific
                              << "alpha: " << alpha
                              << "\nc: " << c
                              << "\nAnorm: " << Anorm
                              << "\nrhobar: " << rhobar
                              << "\n1/rho: " << inv_rho
                              << "\nbeta: " << beta
                              << '\n';
                    std::abort();
                }
                if (iter%100 == 0) {
                    std::cerr << "LSQR iter= " << iter
                              << " nnode= " << nnode << ", delta = " << conv << "\r";
                }
                assert(conv > 0);
                
                if (iter > max_iter) {
                    istop = 4;
                }
                if (conv < tolerance) {
                    istop = 30;
                }
            }
            
            std::cerr << "LSQR iter= " << iter
                      << " nnode= " << nnode << " ndata= " << ndata << "\n";
            max_iter = iter;
            return istop;
        }
    };
}

bool
select_in_house_omp_solver() {
    auto solver = new in_house_omp_solver();
    std::cout << "Installing " << solver->name() << '\n';

    return true;
}
