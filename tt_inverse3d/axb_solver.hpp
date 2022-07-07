#if !defined(TOMO_AXB_SOLVER_HPP)
#define TOMO_AXB_SOLVER_HPP

#include <vector>
#include <string>
#include "sparse_rect.hpp"

class axb_solver {
public:
    axb_solver();
    axb_solver(axb_solver const&) = delete;
    axb_solver& operator=(axb_solver const&) = delete;
    virtual ~axb_solver() = 0;
    
    virtual std::string name() const = 0;
    
    virtual int operator()(SparseRectangular const& a, std::vector<double> const& b,
                           std::vector<double>& x, double tolerance, int& itermax) const;
    
    static axb_solver& instance();
    
    void computation_statistics(bool v) { my_computation_statistics = v; }
    void trace_computation_above(int nb_seconds) {
        my_computation_tracing_limit = nb_seconds; 
    }
    void computation_trace_dir(std::string dirname) {
        my_computation_tracing_storage = dirname;
    }
    void do_not_trace_computations(int nb_seconds) { my_computation_tracing_limit = -1; }

private:
    virtual int solve_impl(SparseRectangular const& a, std::vector<double> const& b,
                           std::vector<double>& x, double tolerance, int& max_iter) const = 0;
    

private:
    bool        my_computation_statistics = false;
    int         my_computation_tracing_limit = -1;
    std::string my_computation_tracing_storage = ".";
    int         my_index = 0;
};

#endif

