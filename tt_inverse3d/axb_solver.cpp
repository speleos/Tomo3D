#include <memory>
#include <mutex>
#include <chrono>
#include <sstream>
#include <fstream>

#include "boost/archive/binary_oarchive.hpp"

#include "axb_solver.hpp"

#include "sparse_rect.hpp"

namespace {
    std::unique_ptr<axb_solver> our_instance;
}

axb_solver::axb_solver() {
    static std::mutex m;
    std::lock_guard<std::mutex> g(m);
    our_instance.reset(this);
}

axb_solver::~axb_solver() {}

axb_solver&
axb_solver::instance()  {
    return *our_instance;
}

int
axb_solver::operator()(SparseRectangular const& a, std::vector<double> const& b,
                       std::vector<double>& x, double tolerance, int& max_iter) const {
    int index = const_cast<axb_solver*>(this)->my_index++;
    auto start = std::chrono::system_clock::now();
    if (my_computation_statistics) {
        std::cout << "Entering solver with     x.sz() = " << x.size() 
                  << ", b.sz() = "      << b.size() << " max iter: " << max_iter 
                  << ", tolerance: " << tolerance << '\n';
    }
    int res = solve_impl(a, b, x, tolerance, max_iter);
    std::chrono::duration<double> whole = std::chrono::system_clock::now() - start;
    if (my_computation_statistics) {
        std::cout << "Solved in "  << whole.count() << "s and " << max_iter << " iter.\n";
    }
    if (my_computation_tracing_limit >0 
        && my_computation_tracing_limit < whole.count()) {
        std::ostringstream dump_fname;
        dump_fname << my_computation_tracing_storage 
                   << "/" << index << "_" << int(whole.count()) << ".dump";
        std::string fname(dump_fname.str());
        std::cout << "Dumping computation into " << fname << "..." << std::flush;
        std::ofstream ofs(fname, std::ios::binary | std::ios::trunc |std::ios::out);
        if (!ofs) {
            std::cerr << "Could not create " << fname << ", bye\n";
            std::abort();
        }
        boost::archive::binary_oarchive ar(ofs);
        ar << a;
        ar << b;
        ar << tolerance;
        ar << max_iter;
        std::cout << "Done.\n";
    }
    return 0;
}
