#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <numeric>

#include <omp.h>

#include "boost/mpi.hpp"
#include "boost/optional.hpp"
#include "boost/program_options.hpp"
#include "boost/archive/binary_iarchive.hpp"

#include "axb_solver.hpp"
#include "in_house_solver.hpp"
#include "in_house_omp_solver.hpp"


boost::optional<double>
compare(std::vector<double> const& x, std::string const& solution_ref_fname) {
    
    std::vector<double> xref(x.size());
    std::ifstream ifs(solution_ref_fname, std::ios::in | std::ios::binary);
    if (!ifs) { 
        std::cerr << "Could not open " << solution_ref_fname << '\n';
        return boost::optional<double>();
    } else {
        ifs.read(reinterpret_cast<char*>(xref.data()), xref.size()*sizeof(double));
        if (ifs.gcount()/sizeof(double) != xref.size()) {
            std::cerr << "Only " << (ifs.gcount()/sizeof(double)) << " fp number could be loaded from "
                      << solution_ref_fname << " while " << x.size() << " were expected.\n";
            return boost::optional<double>();
        } else {
            std::vector<double> diff(x.size());
            std::transform(xref.begin(), xref.end(), x.begin(), diff.begin(),
                           [](double const& d1, double const& d2) { return std::abs(d1-d2); });
            double accdiff = std::accumulate(diff.begin(), diff.end(), 0.0);
            double accx    = std::accumulate(x.begin(), x.end(), 0.0,
                                             [](double d1, double d2) { return std::abs(d1+d2); });
            double accxref = std::accumulate(xref.begin(), xref.end(), 0.0,
                                             [](double d1, double d2) { return std::abs(d1+d2); });
            std::cout << "Accumulated diff is: " << accdiff 
                      << " (avg: "<< (accdiff/diff.size()) << ")\n";
            auto max_iter = std::max_element(diff.begin(), diff.end());
            int max_pos = max_iter - diff.begin();
            std::cout << "Max diff is: " << *max_iter 
                      << " at position " << max_pos
                      << " x: " << x[max_pos] 
                      << ", xref " << xref[max_pos]
                      << '\n';
            return boost::optional<double>(*max_iter);
        }
    }
}

int
main(int argc, char* argv[]) {
    namespace po  = boost::program_options;
    namespace mpi = boost::mpi;
    namespace ar  = boost::archive;
    
    mpi::environment env;
    mpi::communicator world;

    std::string solver_type = "omp";
    std::string archve_fname, solution_fname, solution_ref_fname;
    int nb_threads = 1;
    bool dump_placement = false;
    bool debug_wait_asked = false;

    std::string solver_name      = "default";
    std::string solver_trace_dir = ".";
    bool        solver_stats     = false;
    int         solver_trace_limit = -1;
    std::string archive_fname;
    double user_tolerance = -1;
    int    user_max_iter  = -1;

    po::options_description opts("Among supported options:");
    opts.add_options() 
        ("help,h", "Print this help message.")
        ("placement", "Print placement of the processes at startup.")
        ("wait-db", "wait for the debugger.")
        ("tolerance",  po::value<double>(&user_tolerance), "Override the convergence tolerance.")
        ("max-iterations",  po::value<int>(&user_max_iter), "Override the number of authorized iterations.")
        ("nb-threads", po::value<int>(&nb_threads), "The number of thread to use.")
        ("solver", po::value<std::string>(&solver_name)->default_value("omp"), "select the solver name. Supported values are: \"omp\", \"default\".")
        ("archive,a", po::value<std::string>(&archive_fname), "The archive containing the system to solve.")
        ("solution-output", po::value<std::string>(&solution_fname), "The output file for solution. Not written by default.")
        ("solution-ref", po::value<std::string>(&solution_ref_fname), "The archive containing a reference solution. If provided, comparison will be displayed wrt to the computed solution,")
        ("solver-stats", "print solver time along with problem sizes.");
    
    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv)
              .options(opts).run(),
              vm);
    po::notify(vm);
    dump_placement   = vm.count("placement")   > 0;
    debug_wait_asked = vm.count("wait-db")     > 0;
    solver_stats     = vm.count("solver-stats") > 0;
    if (vm.count("help") > 0) {
        std::cerr << opts;
        return 1;
    }
    
    if (vm.count("archive") == 0) {
        std::cerr << "Error: archive not specified.\n"
                  << opts;
        return 1;
    }
    if (solver_name == "omp")  {
        select_in_house_omp_solver();
    } else if (solver_name == "nothread" || solver_name == "default" )  {
        select_in_house_solver();
    } else {
        std::cerr << "unknown solver " << solver_name << ". Bye\n";
        return -1;
    }
    axb_solver& solver = axb_solver::instance();
    solver.computation_statistics(true);
    std::ifstream ifs(archive_fname, std::ios::binary | std::ios::in);
    if (!ifs) {
        std::cerr << "Could not read " << archive_fname << ", bye\n";
        return -1;
    }
    SparseRectangular   matrix;
    std::vector<double> rhs, x;
    float omega= -1;
    int max_iter = -1;
    float tolerance   = -1;
    boost::archive::binary_iarchive iar(ifs);
    iar >> matrix;
    iar >> rhs;
    iar >> tolerance;
    if (user_tolerance >= 0) {
        tolerance = user_tolerance;
    }
    iar >> max_iter;
    if (user_max_iter >= 0) {
        max_iter = user_max_iter;
    }
    x.resize(matrix.nCol(), 0);
    omp_set_num_threads(nb_threads);
    int status = solver(matrix,rhs,x,tolerance,max_iter);
    if (solution_fname.length() != 0) {
        std::cout << "Writing solution to '" << solution_fname << "'.." << std::flush;
        std::ofstream ofs(solution_fname, std::ios::trunc | std::ios::out | std::ios::binary);
        ofs.write(reinterpret_cast<char*>(x.data()), x.size()*sizeof(double));
        std::cout << "'.done\n";
    }
    if (solution_ref_fname.length() != 0) {
        boost::optional<double> diff = compare(x, solution_ref_fname);
        if (!diff) {
            return -1;
        }
    }
    
    std::cout << "Solver status: " << status << '\n';
    return 0;    
} 
