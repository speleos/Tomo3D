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
    
    double
    normalize(std::vector<double>& x) {
        double norm=0;
        int n = x.size();
        for (int i = 0; i < n; ++i){
            assert(!is_nan(x[i]));
            norm += std::pow(x[i], 2);
        }
        norm = std::sqrt(norm);
        double ratio = 1/norm;
        std::transform(x.begin(), x.end(), x.begin(), std::bind2nd(std::multiplies<double>(), ratio));
        return norm;
    }

    template<class T>
    T const& cst(T& v) {
        return const_cast<T const&>(v);
    }
    
    class in_house_solver 
        : public axb_solver {
    public:
        in_house_solver() = default;
        virtual ~in_house_solver() {}
        
        std::string name() const override { return "<default in house solver>"; }

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
            
            std::vector<double> v;
            std::vector<double> w;
            std::vector<double> u(rshv);
            
            // initialization
            double beta = normalize(u);
            mat.Atx(cst(u),v);
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
                std::transform(u.begin(), u.end(), u.begin(), 
                              [=](double const& d) { return d*nalpha; });
                mat.Ax(cst(v),u,false);
                beta = normalize(u); 
                double nbeta = -beta;
                std::transform(v.begin(), v.end(), v.begin(),
                               [=](double const& d) { return d*nbeta; });
                mat.Atx(cst(u),v,false);
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
select_in_house_solver() {
    auto solver = new in_house_solver();
    std::cout << "Installing " << solver->name() << '\n';
    return true;
}
