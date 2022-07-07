/*
 * inverse3d.h based on inverse.h
 *
 * Adria Melendez and Jun Korenaga
 * Fall 2011
 */

#ifndef _TOMO_INVERSE3D_H_
#define _TOMO_INVERSE3D_H_

#include <iostream>
#include <fstream>
#include <map>
#include <memory>
#include <cassert>

#include "boost/mpi.hpp"
#include "boost/shared_ptr.hpp"
#include "boost/scoped_ptr.hpp"
#include "boost/optional.hpp"

#include "array.hpp"
#include "smesh3d.hpp"
#include "ani3d.hpp"
#include "graph3d.hpp"
#include "betaspline3d.hpp"
#include "bend3d.hpp"
#include "sparse_rect.hpp"
#include "interface3d.hpp"
#include "corrlen3d.hpp"

class TomographicInversion3d {
public:
    enum ray_type {
        RAY_REFRACTION = 0,
        RAY_REFLECTION = 1,
        RAY_MULTIPLE_REFRACTION = 2,
        RAY_MULTIPLE_REFLECTION = 3        
    };

private:
    static int  identity(int i)   { return i; }
    static bool unfiltered(int i) { return true;}
public:

    TomographicInversion3d( boost::mpi::communicator& com, SlownessMesh3d& m, const char* datafn,
                            int xorder, int yorder, int zorder, double crit_len,
                            int nintp, double cg_tol, double br_tol);
    ~TomographicInversion3d() {}
    void solve(int niter);
    void fixing(bool, bool, bool, bool);
    void doRobust(double);
    void removeOutliers();
    void setLSQR_TOL(double);
    void add_reflections(boost::shared_ptr<Interface3d> const& intfp);
    void doFullRefl();
    void setReflWeight(double);

    void tiltAngle(double);

    void add_anisotropy(boost::shared_ptr<AnisotropyMesh3d> const& anidp,
			boost::shared_ptr<AnisotropyMesh3d> const& aniep);

    void usingVh();

    void onlyForward(bool);
    
    void SmoothVelocity(const char*, double, double, double,
			bool logscale=false);
    void SmoothDepth(double, double, double,
		     bool logscale=false);
    void SmoothDepth(const char*, double, double, double,
		     bool logscale=false);
    void SmoothAniD(const char*, double, double, double,
		    bool logscale=false);
    void SmoothAniE(const char*, double, double, double,
		    bool logscale=false);
    void applyFilter(const char*);
    void DampVelocity(double);
    void DampDepth(double);
    void DampAniD(double);
    void DampAniE(double);
    void FixDamping(double, double, double, double);
    void Squeezing(const char*);
    void SqueezingD(const char*);
    void SqueezingE(const char*);
    void targetChisq(double);
    void doJumping();

    Array1d<double> const&  velocity_model() const;
    void filter_velocity_perturbations();

    void set_output_prefix(std::string r) { my_output_prefix = r; }
    void set_verbosity_level(int level) { my_verbosity_level = level;}
    /// \brief Check verbosity level
    /// \param any if true, is verbose whatever the mpi process. Only one otherwise
    bool is_verbose(int level, bool any = false) const;
    void outStepwise(int level);
    void outFinal(int level);
    void set_log(std::string const& fname);
    void outMask();
    void outMaskRefl();
    void outMaskAnid();
    void outMaskAnie();
    void setVerbose(int);

    void printSynTime(ostream&, double) const;
    
    void set_nb_threads(int i) {nb_threads = i; }
    bool use_multithreading() const { return nb_threads > 1; }
    boost::mpi::communicator& com() {return my_com;}
    boost::mpi::communicator const& com() const {return my_com;}
    SlownessMesh3d const& slowness_mesh() const { return my_graph_solver.slowness_mesh(); }
    SlownessMesh3d& slowness_mesh() { /* we *do* own it*/ return const_cast<SlownessMesh3d&>(my_graph_solver.slowness_mesh()); }
    bool need_reflections() const;
    bool xflat() const { return slowness_mesh().xflat(); }
    bool yflat() const { return slowness_mesh().yflat(); }

private:
    void trace_path_processing(char mark, int bend_iter, int src, int rcv, int step) const;
    int read_file(std::string const& ifn);
    void limit_range(GraphSolver3d& g, int src) const;
    void reset_kernel();
    void add_kernel(map<int,double>& A_i, const Array1d<Point3d>&) const;
    void add_kernel_refl(map<int,double>& A_i, const Array1d<Point3d>&, int, int) const;
    void calc_averaging_matrix();
    void calc_damping_matrix();
    void calc_refl_averaging_matrix();
    void calc_refl_damping_matrix();
    void calc_anid_averaging_matrix();
    void calc_anid_damping_matrix();
    void calc_anie_averaging_matrix();
    void calc_anie_damping_matrix();
    void add_global(const Array2d<double>&, const std::vector<int>&,
		    sparse_matrix&);
    int _solve(bool,double,bool,double,bool,double,bool,double,
	       bool,double,bool,double,bool,double,bool,double);
    void auto_damping(int&, int&, double&, double&, double&, double&);
    void fixed_damping(int&, int&, double&, double&, double&, double&);
    double auto_damping_depth(double, double, double, int&, int&);
    double auto_damping_anie(double, double, double, int&, int&);
    double auto_damping_anid(double, double, double, int&, int&);
    double auto_damping_vel(double, double, double, int&, int&);
    double calc_ave_dmv();
    double calc_ave_dmd();
    double calc_ave_dmad();
    double calc_ave_dmae();

    int nb_data() const  { assert(my_nb_data >= 0);  return my_nb_data; }
    int nb_noded() const { assert(my_nb_noded >= 0); return my_nb_noded; }

    typedef boost::optional<double> opt_double;
    void calc_Lm(opt_double& lmx, opt_double& lmy, double& lmz, opt_double& lmdx, opt_double& lmdy,
		 opt_double& lmadx, opt_double& lmady, double& lmadz,
		 opt_double& lmaex, opt_double& lmaey, double& lmaez);
    double calc_chi();
    ray_type get_ray_type(int src, int rcv) const { return my_ray_types[src-1][rcv-1]; }
    ray_type set_ray_type(int src, int rcv, ray_type t) { return my_ray_types[src-1][rcv-1] = t; }
    static bool is_refraction(ray_type rt) { return (rt == RAY_REFRACTION || rt == RAY_MULTIPLE_REFRACTION); }
    static bool is_reflection(ray_type rt) { return (rt == RAY_REFLECTION || rt == RAY_MULTIPLE_REFLECTION); }
    static bool is_multiple(ray_type rt)   { return (rt == RAY_MULTIPLE_REFRACTION || rt == RAY_MULTIPLE_REFLECTION); }
    bool is_refraction(int src, int rcv) const { return is_refraction(get_ray_type(src,rcv)); }
    bool is_reflection(int src, int rcv) const { return is_reflection(get_ray_type(src,rcv)); }
    bool is_multiple(int src, int rcv) const { return is_multiple(get_ray_type(src,rcv)); }
    
    sparse_matrix& rvx() { return *my_rvx; };
    sparse_matrix& rdx() { return *my_rdx; };
    sparse_matrix& radx() { return *my_radx; };
    sparse_matrix& raex() { return *my_raex; };
    sparse_line& rvx(int i) { return rvx()(i); };
    sparse_line& rdx(int i) { return rdx()(i); };
    sparse_line& radx(int i) { return radx()(i); };
    sparse_line& raex(int i) { return raex()(i); };

    sparse_matrix& rvy() { return *my_rvy; };
    sparse_matrix& rdy() { return *my_rdy; };
    sparse_matrix& rady() { return *my_rady; };
    sparse_matrix& raey() { return *my_raey; };
    sparse_line& rvy(int i) { return rvy()(i); };
    sparse_line& rdy(int i) { return rdy()(i); };
    sparse_line& rady(int i) { return rady()(i); };
    sparse_line& raey(int i) { return raey()(i); };
    
    sparse_matrix& rvz() { return my_rvz; };
    sparse_matrix& radz() { return my_radz; };
    sparse_matrix& raez() { return my_raez; };
    sparse_line& rvz(int i) { return rvz()(i); };
    sparse_line& radz(int i) { return radz()(i); };
    sparse_line& raez(int i) { return raez()(i); };

    boost::shared_ptr<AnisotropyMesh3d> anid, anie;

    boost::mpi::communicator& my_com;
    GraphSolver3d   my_graph_solver;
    BetaSpline3d betasp;
    double cg_tolerance;
    double br_tolerance;
    
    Array1d<Point3d> src;
    Array1d< Array1d<Point3d> > rcv;
    std::vector<std::vector<ray_type> > my_ray_types;
    Array1d< Array1d<double> > obs_ttime;
    Array1d< Array1d<double> > obs_dt;
    Array1d< Array1d<double> > res_ttime;
    Array1d< Array1d<double> > syn_ttime;
    Array1d<double> r_dt_vec, path_wt;

    const Interface3d *bathyp;
    boost::shared_ptr<Interface3d> reflp;
    double refl_weight;
    bool do_full_refl;
    
    int nnodev, my_nb_noded, my_nb_data, ndata_valid, nx, ny, nz;
    int nnoded = 0;
    int nnodee = 0;
    double rms_tres[4], init_chi[4], rms_tres_total, init_chi_total;
    int ndata_in[4];
    bool robust;
    double crit_chi;
    sparse_matrix A, Tv, Td, Tad, Tae;
    sparse_matrix                  my_rvz, my_rdz, my_radz, my_raez;
    std::unique_ptr<sparse_matrix> my_rvx, my_rvy, my_rdx, my_rdy,
	my_radx, my_raex, my_rady, my_raey;
    
    Array2d<double> Td_common;
    Array1d<double> data_vec, total_data_vec;
    Array1d<double> modelad, modelae;
    Array1d<double> modeld, dmodel_total;
    Array1d<double> my_velocity_model;
    Array1d<double> my_full_reflection_velocity_model;
    int itermax_LSQR;
    double LSQR_ATOL;

    bool jumping;
    Array1d<double> dmodel_total_sum, mvscale, mdscale, madscale, maescale;

    bool smooth_velocity, logscale_vel;
    double wsv_min, wsv_max, dwsv;
    double weight_s_v;
    CorrelationLength3d *corr_vel_p;
    bool my_apply_filter;
    Interface3d *uboundp;

    bool ani = false;
    double phi;

    bool fv, fd, fad, fae;

    bool useVh = false;

    bool only_forward = false;
    bool first_it = false;

    bool smooth_anid, logscale_anid;
    double wsad_min, wsad_max, dwsad;
    double weight_s_ad;
    CorrelationLength3d *corr_anid_p;

    bool smooth_anie, logscale_anie;
    double wsae_min, wsae_max, dwsae;
    double weight_s_ae;
    CorrelationLength3d *corr_anie_p;

    bool smooth_depth, logscale_dep;
    double wsd_min, wsd_max, dwsd;
    double weight_s_d;
    CorrelationLength2d *corr_dep_p;

    bool damp_velocity, damp_depth;
    bool damp_anid, damp_anie;
    double target_dv, target_dd;
    double target_dad, target_dae;
    bool damping_is_fixed; 
    double weight_d_v, weight_d_d;
    double weight_d_ad, weight_d_ae;
    bool do_squeezing, do_squeezingD, do_squeezingE;
    DampingWeight3d *damping_wt_p, *dampingD_wt_p, *dampingE_wt_p;

    int nnode_total, nnode_subtotal;
    std::vector<int> row_index;
    std::vector<int> tmp_nodedc, tmp_nodedr;
    std::vector<int> tmp_nodeadc, tmp_nodeadr, tmp_nodeaec, tmp_nodeaer;
    bool out_mask, out_mask_refl, out_mask_anid, out_mask_anie;
    Array1d<double> dws, dwsr, dwsda, dwsea;

    ofstream* vmesh_os_p;

    double target_chisq;
    
    bool printTransient;
    bool printFinal;
    std::string my_output_prefix; 
    typedef boost::shared_ptr<std::ostream> ostream_ptr;
    std::string indexed_filename(std::string base, int idx1, int idx2) const;
    ostream_ptr indexed_trace(std::string base, int idx1, int idx2, bool cond = true) const;
    void compute_full_reflection_model();
    
    /// \brief Pseudo randomly elect one MPI process
    bool pick_one_process() const;

    ostream_ptr log_stream() const;

    int my_verbosity_level;
    int my_output_level;

    sparse_matrix B;
    std::string my_log_filename;
    ofstream* dump_os_p;

    class solver;
    class source_processor;
    class source_solver;
    class path_processor;
    class path_solver;
    friend class solver;
    friend class source_solver;
    friend class source_processor;
    friend class path_solver;
    friend class path_processor;

    void compute_paths(solver& g);

    int nb_threads;
    void process_source(source_solver& l) const;
    void process_receiver(source_solver& l, int recv) const;
};

#endif /* _TOMO_INVERSE3D_H_ */
