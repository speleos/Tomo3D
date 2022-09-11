/*
 * inverse3d.cc based on inverse.cc
 * 
 * Adria Melendez and Jun Korenaga
 * Fall 2011
 */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>

#include "boost/scoped_ptr.hpp"
#include "boost/thread.hpp"
#include "boost/optional/optional_io.hpp"

#include "geom3d.hpp"
#include "inverse3d.hpp"
#include "sparse_rect.hpp"
#include "axb_solver.hpp"
#include "d_nr.hpp"

typedef boost::lock_guard<boost::mutex> lock;

boost::mutex&
unconst(boost::mutex const& m) {
    return const_cast<boost::mutex&>(m);
}

class TomographicInversion3d::solver {
public:
    friend class source_solver;
    friend class source_processor;

    solver(TomographicInversion3d& inv, int nb)
        : inversion(inv),
          my_nb_step(nb),
          my_step(1),
          my_graph_time(0),
          my_bend_time(0),
          my_matrix_bits(),
          my_time_mutex(),
          my_next_source(-1),
          my_producer_mutex(),
          my_used_threads(0),
          my_used_threads_mutex()
	{}
    
    int  step()       const { return my_step; }
    int  nb_step()    const { return my_nb_step; }
    bool solving()    const { return my_step >= 1 && my_step <= my_nb_step; }
    bool last_step()  const { assert(solving()); return my_step == my_nb_step; }
    void force_last_step() { assert(solving()); my_nb_step = my_step; }
    bool first_step() const { assert(solving()); return my_step == 1; }
    void next() { assert(solving());  ++my_step; }
    void init_step() {
        my_graph_time = 0;
        my_bend_time  = 0;
        my_next_source = inversion.com().rank();
        // should not be necessary:
        std::vector<sparse_matrix> fresh(inversion.src.size());
        my_matrix_bits.swap(fresh);
    }
    // thread safe:
    void add_time(double g, double b) {
        lock guard(my_time_mutex);
        my_graph_time += g;
        my_bend_time  += b;
    }
    double graph_time() const { return my_graph_time; }
    double bend_time() const { return my_bend_time; }
    std::vector<sparse_matrix>& matrix_bits() { return  my_matrix_bits; }
    
    bool all_sources_processed() const {
        lock guard(unconst(my_producer_mutex));
        return my_next_source >= inversion.src.size();
    }

    int  next_source() {
        lock guard(my_producer_mutex);
        int src;
        if (my_next_source < inversion.src.size()) {
            // +1 for fortran compatibility
            src = my_next_source + 1;
            my_next_source += inversion.com().size();
        } else {
            src = -1;
        }
        return src;
    }
    
    int available_threads() const {
        lock guard(const_cast<boost::mutex&>(my_used_threads_mutex));
        return inversion.nb_threads - my_used_threads;
    }

    int allocate_thread() {
        lock guard(my_used_threads_mutex);
        ++my_used_threads;
        if (my_used_threads > my_max_used_threads) {
            my_max_used_threads = my_used_threads;
        }
        return inversion.nb_threads - my_used_threads;;
    }

    int release_thread() {
        lock guard(my_used_threads_mutex);
        --my_used_threads;
        return inversion.nb_threads - my_used_threads;
    }
    
    struct thread_user {
        solver& owner;
        bool    do_lock;
        thread_user(TomographicInversion3d::solver& s, bool doit = true) 
            : owner(s), do_lock(doit)  { if (do_lock) {owner.allocate_thread();} }
        ~thread_user() { if (do_lock) {owner.release_thread();} }
    };

    TomographicInversion3d&     inversion;
    int max_used_threads() const { return my_max_used_threads; }
    
    void reset_max_used_threads() {
        lock guard(my_used_threads_mutex);
        my_max_used_threads = my_used_threads;
    }

private:
    int                         my_step;
    int                         my_nb_step;
    double                      my_graph_time;
    double                      my_bend_time;
    std::vector<sparse_matrix>  my_matrix_bits;
    boost::mutex                my_time_mutex;
    int                         my_next_source;
    boost::mutex                my_producer_mutex;
    int                         my_used_threads;
    boost::mutex                my_used_threads_mutex;
    int                         my_max_used_threads;
};

struct TomographicInversion3d::source_solver {
    source_solver(solver& g, int s) 
        : global(g), src(s),
          graph_time(0), bend_time(0), res_ttime(), syn_ttime(), matrix_bit(),
          refraction_solver(), reflection_solver(),
          tres_os_mutex(), ray_os_mutex(),
          my_receiver_index(0)
	{}
    
    bool all_receivers_processed() const {
        int const nb_receiver = global.inversion.rcv(src).size();
        lock g(unconst(rcv_mutex)); // could be useless, not sure about the atomicity of this:
        return (my_receiver_index >= nb_receiver);
    }

    int next_receiver() {
        int const nb_receiver = global.inversion.rcv(src).size();
        lock g(rcv_mutex);
        if (my_receiver_index < nb_receiver) {
            ++my_receiver_index;
            // +1 for fortran compatibility...
            return my_receiver_index;
        } else {
            return -1;
        }        
    }
    
    solver&              global;
    int                  src;
    double               graph_time;
    double               bend_time;
    Array1d<double>      res_ttime;
    Array1d<double>      syn_ttime;
    sparse_matrix        matrix_bit;
    typedef boost::shared_ptr<GraphSolver3d const> solver_ptr;
    solver_ptr           refraction_solver;
    solver_ptr           reflection_solver;
    boost::mutex         tres_os_mutex;
    boost::mutex         ray_os_mutex;
    boost::mutex         rcv_mutex;
private:
    int my_receiver_index;
};

struct TomographicInversion3d::path_solver {
    path_solver(source_solver& s, int r) 
        : parent(s), receiver(r) {}
    
    void operator()() {
        parent.global.inversion.process_receiver(parent, receiver);
    }
    source_solver&       parent;
    int                  receiver;
};

class TomographicInversion3d::source_processor {
public:
    source_processor(solver& g) : global(g) {}
    void operator()();

    solver&  global;
};

class TomographicInversion3d::path_processor {
public:
    path_processor(source_solver& s, bool own_thread ) 
        : parent(s), my_own_thread(own_thread) {}

    void operator()();

    source_solver&  parent;
private:
    bool my_own_thread;
};

void
TomographicInversion3d::compute_full_reflection_model() {
    assert (do_full_refl);
    double p_low = 1.0/0.33;
    my_full_reflection_velocity_model = (my_velocity_model);
    Array2d<int> irefl(nx,ny);
    for (int i=1; i<=nx; i++){
        for (int j=1; j<=ny; j++){
            Point3d p = slowness_mesh().nodePos(slowness_mesh().nodeIndex(i,j,1));
            double reflz = reflp->z(p.x(),p.y());
            irefl(i,j) = nz;
            for (int k=1; k<=nz; k++){
                if (slowness_mesh().nodePos(slowness_mesh().nodeIndex(i,j,k)).z()>reflz){
                    irefl(i,j) = k;
                    break;
                }
            }
        }
    }
    for (int i=1; i<=nx; i++){
        for (int j=1; j<=ny; j++){
            for (int k=2; k<=nz; k++){
                int idx = slowness_mesh().nodeIndex(i,j,k);
                int below = slowness_mesh().nodeIndex(i,j,k-1);
                if (k==irefl(i,j)){
                    // this is to make the velocity at the node just below the reflector to
                    // be the same as the velocity just above. If I used pwater for this node,
                    // there would be an unwanted low velocity gradient surrounding the reflector.
                    my_full_reflection_velocity_model(idx) = my_full_reflection_velocity_model(below);
                }else if (k>irefl(i,j)){
                    if (i==1){
                        if (j==1){
                            if ((irefl(2,2)>irefl(1,1) && k<=irefl(2,2)) || (irefl(1,2)>irefl(1,1) && k<=irefl(1,2)) || (irefl(2,1)>irefl(1,1) && k<=irefl(2,1))){
                                my_full_reflection_velocity_model(idx) = my_full_reflection_velocity_model(below);;
                            }else{
                                my_full_reflection_velocity_model(idx) = p_low;
                            }
                        }else if (j==ny){
                            if ((irefl(2,ny)>irefl(1,ny) && k<=irefl(2,ny)) || (irefl(1,ny-1)>irefl(1,ny) && k<=irefl(1,ny-1)) || (irefl(2,ny-1)>irefl(1,ny) && k<=irefl(2,ny-1))){
                                my_full_reflection_velocity_model(idx) = my_full_reflection_velocity_model(below);;
                            }else{
                                my_full_reflection_velocity_model(idx) = p_low;
                            }
                        }else{
                            if ((irefl(1,j-1)>irefl(1,j) && k<=irefl(1,j-1)) || (irefl(1,j+1)>irefl(1,j) && k<=irefl(1,j+1)) || (irefl(2,j)>irefl(1,j) && k<=irefl(2,j)) || (irefl(2,j-1)>irefl(1,j) && k<=irefl(2,j-1)) || (irefl(2,j+1)>irefl(1,j) && k<=irefl(2,j+1))){
                                my_full_reflection_velocity_model(idx) = my_full_reflection_velocity_model(below);;
                            }else{
                                my_full_reflection_velocity_model(idx) = p_low;
                            }
                        }
                    }else if (i==nx){
                        if (j==1){
                            if ((irefl(nx,2)>irefl(nx,1) && k<=irefl(nx,2)) || (irefl(nx-1,2)>irefl(nx,1) && k<=irefl(nx-1,2)) || (irefl(nx-1,1)>irefl(nx,1) && k<=irefl(nx-1,1))){
                                my_full_reflection_velocity_model(idx) = my_full_reflection_velocity_model(below);;
                            }else{
                                my_full_reflection_velocity_model(idx) = p_low;
                            }
                        }else if (j==ny){
                            if ((irefl(nx,ny-1)>irefl(nx,ny) && k<=irefl(nx,ny-1)) || (irefl(nx-1,ny-1)>irefl(nx,ny) && k<=irefl(nx-1,ny-1)) || (irefl(nx-1,ny)>irefl(nx,ny) && k<=irefl(nx-1,ny))){
                                my_full_reflection_velocity_model(idx) = my_full_reflection_velocity_model(below);;
                            }else{
                                my_full_reflection_velocity_model(idx) = p_low;
                            }
                        }else{
                            if ((irefl(nx,j-1)>irefl(nx,j) && k<=irefl(nx,j-1)) || (irefl(nx,j+1)>irefl(nx,j) && k<=irefl(nx,j+1)) || (irefl(nx-1,j)>irefl(nx,j) && k<=irefl(nx-1,j)) || (irefl(nx-1,j-1)>irefl(nx,j) && k<=irefl(nx-1,j-1)) || (irefl(nx-1,j+1)>irefl(nx,j) && k<=irefl(nx-1,j+1))){
                                my_full_reflection_velocity_model(idx) = my_full_reflection_velocity_model(below);;
                            }else{
                                my_full_reflection_velocity_model(idx) = p_low;
                            }
                        }
                    }else{
                        if (j==1){
                            if ((irefl(i-1,1)>irefl(i,1) && k<=irefl(i-1,1)) || (irefl(i+1,1)>irefl(i,1) && k<=irefl(i+1,1)) || (irefl(i-1,2)>irefl(i,1) && k<=irefl(i-1,2)) || (irefl(i+1,2)>irefl(i,1) && k<=irefl(i+1,2)) || (irefl(i,2)>irefl(i,1) && k<=irefl(i,2))){
                                my_full_reflection_velocity_model(idx) = my_full_reflection_velocity_model(below);;
                            }else{
                                my_full_reflection_velocity_model(idx) = p_low;
                            }
                        }else if (j==ny){
                            if ((irefl(i-1,ny)>irefl(i,ny) && k<=irefl(i-1,ny)) || (irefl(i+1,ny)>irefl(i,ny) && k<=irefl(i+1,ny)) || (irefl(i-1,ny-1)>irefl(i,ny) && k<=irefl(i-1,ny-1)) || (irefl(i+1,ny-1)>irefl(i,ny) && k<=irefl(i+1,ny-1)) || (irefl(i,ny-1)>irefl(i,ny) && k<=irefl(i,ny-1))){
                                my_full_reflection_velocity_model(idx) = my_full_reflection_velocity_model(below);;
                            }else{
                                my_full_reflection_velocity_model(idx) = p_low;
                            }
                        }else{
                            if((irefl(i-1,j-1)>irefl(i,j) && k<=irefl(i-1,j-1)) || (irefl(i+1,j+1)>irefl(i,j) && k<=irefl(i+1,j+1)) || (irefl(i+1,j-1)>irefl(i,j) && k<=irefl(i+1,j-1)) || (irefl(i-1,j+1)>irefl(i,j) && k<=irefl(i-1,j+1)) || (irefl(i-1,j)>irefl(i,j) && k<=irefl(i-1,j)) || (irefl(i+1,j)>irefl(i,j) && k<=irefl(i+1,j)) || (irefl(i,j-1)>irefl(i,j) && k<=irefl(i,j-1)) || (irefl(i,j+1)>irefl(i,j) && k<=irefl(i,j+1))){
                                my_full_reflection_velocity_model(idx) = my_full_reflection_velocity_model(below);;
                            }else{
                                my_full_reflection_velocity_model(idx) = p_low;
                            }
                        }
                    }
                }
            }
        }
    }
}

void
TomographicInversion3d::compute_paths(solver& globals) {
    if (do_full_refl) {
        compute_full_reflection_model();
    }
    if (use_multithreading()) {
        globals.reset_max_used_threads();
        boost::thread_group workers;
        for (int n = 0; n < src.size() && globals.available_threads() > 0; ++n) {
            workers.create_thread(source_processor(globals));
        }
        workers.join_all();
        if (is_verbose(0) && pick_one_process()) {
            std::cerr << "MPI process " << com().rank() << " used " << globals.max_used_threads() 
                      << " sfw threads to compute paths\n";
        }
    } else {
        source_processor proc(globals);
        proc();
    }
    {
        std::vector<sparse_matrix>& matrix_bits = globals.matrix_bits();
        int flat = 1;
        for( int m = 0; m < src.size(); ++m) {
	  int proc_of_origin = m%com().size();
	  if(!only_forward || first_it){
            if (com().rank() == 0) {
                std::cout << "Broadcasting paths from source " << m << '\n';
            }
	    sparse_matrix& bit = matrix_bits[m];
            boost::mpi::broadcast(com(), bit.stl(), proc_of_origin);
            //re agregate the sparse matrix bits
            for (int j = 1; j <= bit.size(); ++j ) {
                A(flat++).swap(bit(j));
            }
	  }
	  boost::mpi::broadcast(com(), res_ttime.stl()[m].stl(), proc_of_origin);
	  boost::mpi::broadcast(com(), syn_ttime.stl()[m].stl(), proc_of_origin);
        }
    }
}

void
TomographicInversion3d::filter_velocity_perturbations() {
    assert(my_apply_filter);
    if (is_verbose(0)) {
        cerr << "filtering velocity perturbation...";
    }
    Array2d<int> kstart;
    Array1d<double> dx_vec, dxr_vec, dxn_vec;
    dx_vec.resize(dmodel_total.size()); // model perturbation vector (dm)
    dxr_vec.resize(dmodel_total.size());
    dxn_vec.resize(dmodel_total.size());
    kstart.resize(nx,ny);
    for (int i=1; i<=nx; i++){
        for(int j=1; j<=ny; j++){
            Point3d p = slowness_mesh().nodePos(slowness_mesh().nodeIndex(i,j,1));
            double boundz = uboundp->z(p.x(),p.y());
            kstart(i,j) = 1;
            for (int k=1; k<=nz; k++){
                if (slowness_mesh().nodePos(slowness_mesh().nodeIndex(i,j,k)).z()>boundz) {
                    break;
                }
                kstart(i,j) = k;
            }
        }
    }
    Array1d<double> filtered_model(dmodel_total);
    for (int i=1; i<=nx; i++){
        for (int j=1; j<=ny; j++){
            for (int k=kstart(i,j); k<=nz; k++){
                int anode = slowness_mesh().nodeIndex(i,j,k);
                Point3d p = slowness_mesh().nodePos(anode);
                double Lx, Ly, Lz;
                corr_vel_p->at(p,Lx,Ly,Lz);
                
                double Lx2 = Lx*Lx;
                double Ly2 = Ly*Ly;
                double Lz2 = Lz*Lz;
                double sum=0.0, beta_sum=0.0;
                for (int ii=1; ii<=nx; ii++){
                    boost::optional<double> xexp = slowness_mesh().xexp(p,ii,j,k,Lx);
                    if (xexp) {
                        for (int jj=1; jj<=ny; jj++){
                            boost::optional<double> yexp = slowness_mesh().yexp(p,ii,jj,k,Ly);
                            if (yexp){
                                for (int kk=kstart(ii,jj); kk<=nz; kk++){
                                    int dnode = slowness_mesh().nodeIndex(ii,jj,kk);
                                    boost::optional<double> zexp = slowness_mesh().zexp(p,dnode,Ly);
                                    if (zexp){
                                        double beta = (*xexp)*(*yexp)*(*zexp);
                                        sum += beta*dmodel_total(dnode);
                                        beta_sum += beta;
                                    }
                                }
                            }
                        }
                    } else {
                        assert(!slowness_mesh().xexp(p,ii,j,k,Lx));
                    }
                }
                filtered_model(anode) = sum/beta_sum;
            }
        }
    }
    if (is_verbose(0)) cerr << "done.\n";
    
    // conservative filtering (Deal and Nolet, GJI, 1996)
    if (is_verbose(0)) {
        cerr << "calculating a nullspace shuttle...";
    }
    dx_vec=0.0;
    for (int i=1; i<=nnodev; i++) {
        dx_vec(i) = filtered_model(i)-dmodel_total(i);
    }
    Array1d<double> Adm(ndata_valid);
    SparseRectangular sparseA(A,[&](int i){return bool(row_index[i-1]);},identity,nnode_total);
    sparseA.Ax(dx_vec,Adm);
    
    int lsqr_itermax=itermax_LSQR;
    dxr_vec=0.0;
    axb_solver::instance()(sparseA,Adm,dxr_vec,LSQR_ATOL,lsqr_itermax);
    
    // note: correction_vec for depth nodes should be zero, so I don't use them.
    double orig_norm=0, new_norm=0, iprod=0;
    for (int i=1; i<=nnodev; i++){
        orig_norm += filtered_model(i)*filtered_model(i);
        dmodel_total(i) = filtered_model(i)-dxr_vec(i);
        new_norm += dmodel_total(i)*dmodel_total(i);
        dxn_vec(i) = dx_vec(i)-dxr_vec(i);
        iprod += dxr_vec(i)*dxn_vec(i);
    }
    if (is_verbose(0)) {
        cerr << "done.\n";
    }
    {
        ostream_ptr log = log_stream();
        if (log){
            *log << "# a posteriori filter check: "
                 << sqrt(orig_norm/nnodev)
                 << " " << sqrt(new_norm/nnodev) 
                 << " " << iprod << endl;
        }
    }
}

void
TomographicInversion3d::solve(int niter) 
{
    typedef map<int,double>::iterator mapIterator;
    typedef map<int,double>::const_iterator mapBrowser;
    {
        ostream_ptr log = log_stream();
        if (log){
            *log << "# strategy " << jumping << " " << robust << " " << crit_chi << '\n'
                 << "# ray_trace " << my_graph_solver.xOrder() << " " << my_graph_solver.yOrder() << " " << my_graph_solver.zOrder() << " "
                 << my_graph_solver.critLength() << " "
                 << betasp.numIntp() << " " << cg_tolerance << '\n'
                 << "# smooth_vel " << smooth_velocity << " "
                 << wsv_min << " " << wsv_max << " " << dwsv << " "
                 << my_apply_filter << '\n'
                 << "# smooth_dep " << smooth_depth << " "
                 << wsd_min << " " << wsd_max << " " << dwsd << '\n'
		 << "# smooth_anid " << smooth_anid << " "
		 << wsad_min << " " << wsad_max << " " << dwsad << '\n'
		 << "# smooth_anie " << smooth_anie << " "
		 << wsae_min << " " << wsae_max << " " << dwsae << '\n';
            if (damping_is_fixed) {
                *log << "# fixed_damping " << weight_d_v << " " << weight_d_d 
		     << " " << weight_d_ad << " " << weight_d_ae << " " << '\n';
            } else {
                *log << "# damp_vel " << damp_velocity << " " << target_dv << '\n';
                *log << "# damp_dep " << damp_depth << " " << target_dd << '\n';
		*log << "# damp_anid " << damp_anid << " " << target_dad << '\n';
		*log << "# damp_anie " << damp_anie << " " << target_dae << '\n';
            }
            *log << "# ndata " << nb_data() << '\n';
            *log << "# nnodev " << nnodev << " # nnoded " << nb_noded() << " # w " << refl_weight << '\n';
	    *log << "# nnodead " << nnoded << " # nnodeae " << nnodee << '\n'
                 << "# LSQR " << LSQR_ATOL << endl;
        }
    }
    // construct index mappers
    if (reflp){
        for (int i=1; i<=nb_noded(); i++){
            int i_total = i+nnodev;
            tmp_nodedc[i-1] = i_total;
            tmp_nodedr[i-1] = i;
        }
    }

    if (ani){
      for (int j=1; j<=nnoded; j++){
	int j_total = j+nnode_subtotal;
	tmp_nodeadc[j-1] = j_total;
	tmp_nodeadr[j-1] = j;
      }
      for (int k=1; k<=nnodee; k++){
	int k_total = k+nnode_subtotal+nnoded;
	tmp_nodeaec[k-1] = k_total;
	tmp_nodeaer[k-1] = k;
      }
    }

    if (jumping) {
        dmodel_total_sum.resize(dmodel_total.size());
    }

    {
        solver solver_global(*this,niter);

        while(solver_global.solving()) {
            solver_global.init_step();
            
            if (is_verbose(0)) {
		if (!only_forward || first_it){
		    cerr << "TomographicInversion3d::iter="
			 << solver_global.step() << "(" << solver_global.nb_step() << ")\n";
		}
            }
            slowness_mesh().get(my_velocity_model);
            if (reflp) {
                reflp->get(modeld);
            }
	    if (ani) {
	      anid->get(modelad);
	      anie->get(modelae);
	    }
            if (solver_global.first_step() || !jumping){ // set model scaling vectors
                mvscale = velocity_model();
                if (reflp) {
                    mdscale = modeld;
                }
		if(ani) {
		  madscale = modelad;
		  maescale = modelae;
		}
            }
            reset_kernel();
            compute_paths(solver_global);

            // construct data vector
            if (!only_forward || first_it) {
                int idata = 1;
                rms_tres[0]=rms_tres[1]=rms_tres[2]=rms_tres[3]=0;
                init_chi[0]=init_chi[1]=init_chi[2]=init_chi[3]=0;
                ndata_in[0]=ndata_in[1]=ndata_in[2]=ndata_in[3]=0;
                ndata_valid=0;
                std::fill(row_index.begin(), row_index.end(), 0);
                for (int isrc=1; isrc<=src.size(); isrc++){
                    for (int ircv=1; ircv<=rcv(isrc).size(); ircv++){
                        double res=res_ttime(isrc)(ircv);
                        data_vec(idata) = res;
                        double res2=res*r_dt_vec(idata);
                        double res22=res2*res2;
                        ray_type icode = get_ray_type(isrc, ircv);
                        rms_tres[icode] += res*res;
                        init_chi[icode] += res22;
                        ++ndata_in[icode];
                        row_index[idata-1] = ++ndata_valid;
                        idata++;
                    }
                }
            
		rms_tres_total = sqrt((rms_tres[0]+rms_tres[1]+rms_tres[2]+rms_tres[3])/ndata_valid);
		init_chi_total = (init_chi[0]+init_chi[1]+init_chi[2]+init_chi[3])/ndata_valid;
		for (int i=0; i<=3; i++){
		    rms_tres[i] = ndata_in[i]>0 ? sqrt(rms_tres[i]/ndata_in[i]) : 0.0;
		    init_chi[i] = ndata_in[i]>0 ? init_chi[i]/ndata_in[i] : 0.0;
		}
		if (init_chi_total<target_chisq) {
		    solver_global.force_last_step();
		}

		// rescale kernel A and data vector
		// note: averaging matrices are scaled upon their construction.
		for (int i=1; i<=nb_data(); i++){
		    double data_wt=r_dt_vec(i);
		    data_vec(i) *= data_wt;
		    for (mapIterator p=A(i).begin(); p!=A(i).end(); p++){
			int j = p->first;
			double m;
			if (j<=nnodev){
			    m = mvscale(j);
			}else if (j>nnodev && j<=nnode_subtotal){
			    m = mdscale(j-nnodev)*refl_weight;
			}else if (j>nnode_subtotal && j<=nnode_subtotal+nnoded){
			    if (madscale(j-nnode_subtotal) != 0.0){
				m = madscale(j-nnode_subtotal);
			    }else{
				m = 1;
			    }
			}else if (j>nnode_subtotal+nnoded && j<=nnode_subtotal+nnoded+nnodee){
			    if (maescale(j-nnode_subtotal-nnoded) != 0.0){
				m = maescale(j-nnode_subtotal-nnoded);
			    }else{
				m = 1;
			    }
			}
			p->second *= m*data_wt;
		    }
		}
		
		// construct total kernel matrix and data vector, and solve Ax=b
		if (smooth_velocity && !fv) calc_averaging_matrix();
		if (smooth_anid && !fad){
		    calc_anid_averaging_matrix();
		}
		if (smooth_anie && !fae){
		    calc_anie_averaging_matrix();
		}
		if (reflp && smooth_depth && !fd) calc_refl_averaging_matrix();
		if (damp_velocity && !fv) calc_damping_matrix();
		if (reflp && damp_depth && !fd) calc_refl_damping_matrix();
		if (damp_anid && !fad){
		    calc_anid_damping_matrix();
		}
		if (damp_anie && !fae){
		    calc_anie_damping_matrix();
		}

		int iset=0;
		for (double tmp_wsv=wsv_min; tmp_wsv<=wsv_max; tmp_wsv+=dwsv){
		    for (double tmp_wsd=wsd_min; tmp_wsd<=wsd_max; tmp_wsd+=dwsd){
			for (double tmp_wsad=wsad_min; tmp_wsad<=wsad_max; tmp_wsad+=dwsad){
			    for (double tmp_wsae=wsae_min; tmp_wsae<=wsae_max; tmp_wsae+=dwsae){
				iset++;
				
				weight_s_v = logscale_vel ? pow(10.0,tmp_wsv) : tmp_wsv;
				weight_s_d = logscale_dep ? pow(10.0,tmp_wsd) : tmp_wsd;
				weight_s_ad = logscale_anid ? pow(10.0,tmp_wsad) : tmp_wsad;
				weight_s_ae = logscale_anie ? pow(10.0,tmp_wsae) : tmp_wsae;
				
				double wdv,wdd,wdad,wdae;
				int nlsqr=0, lsqr_iter=0;
				clock_t start_t = clock();
				if (damping_is_fixed){
				    fixed_damping(lsqr_iter,nlsqr,wdv,wdd,wdad,wdae);
				}else{
				    auto_damping(lsqr_iter,nlsqr,wdv,wdd,wdad,wdae);
				}
				clock_t end_t = clock();
				double lsqr_time = (end_t-start_t);
				lsqr_time /= CLOCKS_PER_SEC;
				
				if (my_apply_filter){
				    filter_velocity_perturbations();
				}
                
				// take stats
				double pred_chi = calc_chi();
				double dv_norm;
				if(fv){
				    dv_norm = -1.0;
				}else{
				    dv_norm = calc_ave_dmv();
				}
				double dd_norm;
				if(fd){
				    dd_norm = -1.0;
				}else{
				    dd_norm = calc_ave_dmd();
				}
				double dad_norm;
				if(fad){
				    dad_norm = -1.0;
				}else{
				    dad_norm = calc_ave_dmad();
				}
				double dae_norm;
				if(fae){
				    dae_norm = -1.0;
				}else{
				    dae_norm = calc_ave_dmae();
				}
				opt_double Lmvx, Lmvy, Lmdx, Lmdy;
				opt_double Lmadx, Lmady, Lmaex, Lmaey;
				double Lmvz = std::numeric_limits<double>::quiet_NaN();
				double Lmadz = std::numeric_limits<double>::quiet_NaN();
				double Lmaez = std::numeric_limits<double>::quiet_NaN();

   cout << "Entering calc_LM " << endl; //ALOUREIRO

				calc_Lm(Lmvx,Lmvy,Lmvz,Lmdx,Lmdy,Lmadx,Lmady,Lmadz,Lmaex,Lmaey,Lmaez);
				
				if (jumping) dmodel_total_sum += dmodel_total;
				
				// scale to slowness perturbation
				if (!fv){
				    Array1d<double> mv(velocity_model());
				    for (int i=1; i<=nnodev; i++){
					mv(i) += mvscale(i)*dmodel_total(i);
				    }
				    slowness_mesh().set(mv);
				}
				if (reflp && !fd){
				    Array1d<double> md(modeld);
				    for (int i=1; i<=nb_noded(); i++){
					int tmp_i = tmp_nodedc[i-1];
					md(i) += mdscale(i)*dmodel_total(tmp_i)*refl_weight;
				    }
				    reflp->set(md);
				}
				if (ani){
				    if(!fad){
					Array1d<double> mad(modelad);
					for (int i=1; i<=nnoded; i++){
					    int tmp_ii = tmp_nodeadc[i-1];
					    mad(i) += madscale(i)*dmodel_total(tmp_ii);
					}
					anid->set(mad);
				    }
				    if(!fae){
					Array1d<double> mae(modelae);
					for (int i=1; i<=nnodee; i++){
					    int tmp_iii = tmp_nodeaec[i-1];
					    mae(i) += maescale(i)*dmodel_total(tmp_iii);
					}
					anie->set(mae);
				    }
				}
				
				if (printTransient || (printFinal && solver_global.last_step())){
				    if (0 % com().size() == com().rank()) {
					if (!fv){
					    slowness_mesh().outMesh(*indexed_trace("smesh",  solver_global.step(), iset));
         std::cout << "##################################################" << std::endl;
         std::cout << "Wrote final Slowness Mesh!" << std::endl;
         cout << "rms =" << rms_tres_total <<"\n"; //ALOUREIRO
         cout << "chi2 =" << init_chi_total <<"\n"; //ALOUREIRO
         std::cout << "##################################################" << std::endl << std::endl;

					}
				    }
				    if (1 % com().size() == com().rank()) {
					if (reflp && !fd){
					    *indexed_trace("refl",  solver_global.step(), iset) << *reflp;
					}
				    }
				    if (2 % com().size() == com().rank()) {
					if (ani && !fad){
					    anid->outMesh(*indexed_trace("anid", solver_global.step(), iset));
					}
				    }
				    if (3 % com().size() == com().rank()) {
					if (ani && !fae){
					    anie->outMesh(*indexed_trace("anie", solver_global.step(), iset));
					}
				    }
				    if (4 % com().size() == com().rank()) {
					if (out_mask){
					    dws.resize(nnodev);
					    dws=0.0;
					    typedef map<int,double>::iterator mapIterator;
					    for (int i=1; i<=A.size(); i++){
						for (mapIterator p=A(i).begin(); p!=A(i).end(); p++){
						    int inode=p->first;
						    if (inode<=dws.size()) dws(inode) += p->second;
						}
					    }
					    slowness_mesh().printMaskGrid(*indexed_trace("dws", solver_global.step(), iset), dws);
					}
				    }
				    if (5 % com().size() == com().rank()) {
					if (reflp && out_mask_refl){
					    dwsr.resize(nb_noded());
					    dwsr=0.0;
					    typedef map<int,double>::iterator mapIterator;
					    for (int i=1; i<=A.size(); i++){
						for (mapIterator p=A(i).begin(); p!=A(i).end(); p++){
						    int inode=p->first;
						    if (inode>nnodev && inode<=nnodev+dwsr.size()) dwsr(inode-nnodev) += p->second;
						}
					    }
					    reflp->printMaskGridRefl(*indexed_trace("dwsr", solver_global.step(), iset), dwsr);
					}
				    }
				    if (6 % com().size() == com().rank()) {
					if (out_mask_anid){
					    dwsda.resize(nnoded);
					    dwsda=0.0;
					    typedef map<int,double>::iterator mapIterator;
					    for (int i=1; i<=A.size(); i++){
						for (mapIterator p=A(i).begin(); p!=A(i).end(); p++){
						    int inode=p->first;
						    if (inode>nnodev+nb_noded() && inode<=nnodev+nb_noded()+nnoded) dwsda(inode-nnodev-nb_noded()) += fabs(p->second);
						}
					    }
					    anid->printMaskGrid(*indexed_trace("dwsad", solver_global.step(), iset), dwsda);
					}
				    }
				    if (7 % com().size() == com().rank()) {
					if (out_mask_anie){
					    dwsea.resize(nnodee);
					    dwsea=0.0;
					    typedef map<int,double>::iterator mapIterator;
					    for (int i=1; i<=A.size(); i++){
						for (mapIterator p=A(i).begin(); p!=A(i).end(); p++){
						    int inode=p->first;
						    if (inode>nnodev+nb_noded()+nnoded && inode<=nnodev+nb_noded()+nnoded+nnodee) dwsea(inode-nnodev-nb_noded()-nnoded) += fabs(p->second);
						}
					    }
					    anie->printMaskGrid(*indexed_trace("dwsae", solver_global.step(), iset), dwsea);
					}
				    }
				}
				{
				    ostream_ptr log = log_stream();
				    if (log){
					double graph_time = solver_global.graph_time() / CLOCKS_PER_SEC;
					double bend_time  = solver_global.bend_time()  / CLOCKS_PER_SEC;
					*log << solver_global.step() << " " << iset << " "
					     << nb_data()-ndata_valid << " "
					     << rms_tres_total << " " << init_chi_total << " "
					     << ndata_in[0] << " " << rms_tres[0] << " " << init_chi[0] << " "
					     << ndata_in[1] << " " << rms_tres[1] << " " << init_chi[1] << " "
					     << ndata_in[2] << " " << rms_tres[2] << " " << init_chi[2] << " "
					     << ndata_in[3] << " " << rms_tres[3] << " " << init_chi[3] << " "
					     << graph_time << " " << bend_time << " " 
					     << weight_s_v << " " << weight_s_d << " "
					     << weight_s_ad << " " << weight_s_ae << " "
					     << wdv << " " << wdd << " " << wdad << " " << wdae << " "
					     << nlsqr << " " << lsqr_iter << " " << lsqr_time << " "
					     << pred_chi << " " << dv_norm << " " << dd_norm << " "
					     << dad_norm << " " << dae_norm << " "
					     << Lmvx << " " << Lmvy << " " << Lmvz << " " << Lmdx << " " << Lmdy << " "
					     << Lmadx << " " << Lmady << " " << Lmadz << " "
					     << Lmaex << " " << Lmaey << " " << Lmaez;
					*log << endl;
				    }
				}
			    }
			}
		    }
		}
	    }
	    solver_global.next();
	}
    }
}

void
TomographicInversion3d::limit_range(GraphSolver3d& g, int isrc) const {
    // limit range first
    double xmin = slowness_mesh().xmax();
    double xmax = slowness_mesh().xmin();
    double ymin = slowness_mesh().ymax();
    double ymax = slowness_mesh().ymin();
    std::vector<Point3d> const& receivers =  rcv(isrc).stl();
    for (std::vector<Point3d>::const_iterator i = receivers.begin();
	 i != receivers.end(); ++i) {
	double x = i->x();
	double y = i->y();
	xmin = std::min(x,xmin);
	xmax = std::max(x,xmax);
	ymin = std::min(y,ymin);
	ymax = std::max(y,ymax);
    }
    {
	double x = src(isrc).x();
	double y = src(isrc).y();
	xmin = std::min(x,xmin);
	xmax = std::max(x,xmax);
	ymin = std::min(y,ymin);
	ymax = std::max(y,ymax);
    }
    g.limitRange(xmin,xmax,ymin,ymax);
    if (is_verbose(1)){
	cerr << "xrange=(" << xmin << "," << xmax << ") "
	     << "yrange=(" << ymin << "," << ymax << ") ";
    }
}

void
TomographicInversion3d::trace_path_processing(char mark, int bend_iter, int src, int rcv, int step) const {
    if (is_verbose(0)) {
	cerr << mark;
	if (is_verbose(1)) {
	    cerr << '(' << bend_iter << ')';
	}
    }
    if (bend_iter < 0) {
	cerr << "TomographicInversion3d::iteration = " << step << '\n'
	     << "TomographicInversion3d::too many iterations required (" << bend_iter << ")\n"  //ALOUREIRO
	     << "TomographicInversion3d::for bending refinement of ray between (s,r)="
             << src << "," << rcv << " - " << shot_point_list(src)(rcv) << '\n';
    }
}

void
TomographicInversion3d::process_receiver(source_solver& locl, int ircv) const {
    int isrc = locl.src;
    Array1d<Point3d> const& receivers =  rcv(isrc);
    Point3d r = receivers(ircv);
    double orig_t, final_t;
    int iterbend;
    Array1d<Point3d> path;
    std::vector<int> start_i, end_i;
    Array1d<const Interface3d*> interp;
    
    switch (get_ray_type(isrc,ircv)) {
    case RAY_REFRACTION:
    case RAY_MULTIPLE_REFRACTION: 
    {
	GraphSolver3d const& solver = *(locl.refraction_solver);
	if (solver.slowness_mesh().inWater(r) && !solver.slowness_mesh().inWater(src(isrc))){
	    int i0, i1;
	    solver.pickPathRcvW(r,path,i0,i1);
	    start_i.push_back(i0);
	    end_i.push_back(i1); 
	    interp.push_back(bathyp);
	    BendingSolver3d bend(solver.slowness_mesh(), betasp, cg_tolerance, br_tolerance);
	    iterbend = bend.refine(path, orig_t, final_t,
				   start_i, end_i, interp, ani);
         //cerr << "iterbend= " << iterbend  <<'\n'; //ALOUREIRO
	    trace_path_processing('*', iterbend, isrc, ircv, locl.global.step());
	} else if (!solver.slowness_mesh().inWater(r) && solver.slowness_mesh().inWater(src(isrc))){
	    int is0, is1;
	    solver.pickPathSrcW(r,path,is0,is1);
	    start_i.push_back(is0);
	    end_i.push_back(is1); 
	    interp.push_back(bathyp);
	    BendingSolver3d bend(solver.slowness_mesh(), betasp, cg_tolerance, br_tolerance);
	    iterbend = bend.refine(path, orig_t, final_t,
				   start_i, end_i, interp, ani);
	    trace_path_processing('^', iterbend, isrc, ircv, locl.global.step());
	} else if (solver.slowness_mesh().inWater(r) && solver.slowness_mesh().inWater(src(isrc))){
	    int i0, i1, is0, is1;
	    solver.pickPathSrcWRcvW(r,path,i0,i1,is0,is1);
	    start_i.push_back(i0);
	    start_i.push_back(is0);
	    end_i.push_back(i1);
	    end_i.push_back(is1); 
	    interp.push_back(bathyp);
	    interp.push_back(bathyp);
	    BendingSolver3d bend(solver.slowness_mesh(), betasp, cg_tolerance, br_tolerance);
	    iterbend = bend.refine(path, orig_t, final_t,
				   start_i, end_i, interp, ani);
	    trace_path_processing('~', iterbend, isrc, ircv, locl.global.step());
	} else {
	    solver.pickPath(r, path);
	    BendingSolver3d bend(solver.slowness_mesh(), betasp, cg_tolerance, br_tolerance);
	    iterbend = bend.refine(path, orig_t, final_t, ani);
	    trace_path_processing('.', iterbend, isrc, ircv, locl.global.step());
	}
	add_kernel(locl.matrix_bit(ircv),path);
    }
    break;
    case RAY_REFLECTION:
    case RAY_MULTIPLE_REFLECTION:
    {
	GraphSolver3d const& solver = *(locl.reflection_solver);
	assert(reflp);
	int ir0, ir1;
	if (solver.slowness_mesh().inWater(r) && !solver.slowness_mesh().inWater(src(isrc))){
	    int i0, i1;
	    solver.pickReflPathRcvW(r,path,i0,i1,ir0,ir1);
	    start_i.push_back(i0);
	    start_i.push_back(ir0);
	    end_i.push_back(i1);
	    end_i.push_back(ir1);
	    interp.push_back(bathyp);
	    interp.push_back(reflp.get());
	    BendingSolver3d bend(solver.slowness_mesh(), betasp, cg_tolerance, br_tolerance);
	    iterbend = bend.refine(path, orig_t, final_t,
				   start_i, end_i, interp, ani);
	    trace_path_processing('#', iterbend, isrc, ircv, locl.global.step());
	} else if (!solver.slowness_mesh().inWater(r) && solver.slowness_mesh().inWater(src(isrc))){
	    int is0, is1;
	    solver.pickReflPathSrcW(r,path,is0,is1,ir0,ir1);
	    start_i.push_back(ir0);
	    start_i.push_back(is0);
	    end_i.push_back(ir1);
	    end_i.push_back(is1);
	    interp.push_back(reflp.get());
	    interp.push_back(bathyp);
	    BendingSolver3d bend(solver.slowness_mesh(),betasp,cg_tolerance,br_tolerance);
	    iterbend = bend.refine(path, orig_t, final_t,
				   start_i, end_i, interp, ani);
	    trace_path_processing('/', iterbend, isrc, ircv, locl.global.step());
	} else if (solver.slowness_mesh().inWater(r) && solver.slowness_mesh().inWater(src(isrc))){
	    int i0, i1, is0, is1;
	    solver.pickReflPathSrcWRcvW(r,path,i0,i1,is0,is1,ir0,ir1);
	    start_i.push_back(i0);
	    start_i.push_back(ir0);
	    start_i.push_back(is0);
	    end_i.push_back(i1);
	    end_i.push_back(ir1);
	    end_i.push_back(is1);
	    interp.push_back(bathyp);
	    interp.push_back(reflp.get());
	    interp.push_back(bathyp);
	    BendingSolver3d bend(solver.slowness_mesh(), betasp, cg_tolerance, br_tolerance);
	    iterbend = bend.refine(path, orig_t, final_t,
				   start_i, end_i, interp, ani);
	    trace_path_processing('|', iterbend, isrc, ircv, locl.global.step());
	} else {
	    solver.pickReflPath(r,path,ir0,ir1);
	    start_i.push_back(ir0);
	    end_i.push_back(ir1);
	    interp.push_back(reflp.get());
	    BendingSolver3d bend(solver.slowness_mesh(),betasp,cg_tolerance,br_tolerance);
	    iterbend = bend.refine(path, orig_t, final_t,
				   start_i, end_i, interp, ani);
	    trace_path_processing('+', iterbend, isrc, ircv, locl.global.step());
	}
	add_kernel_refl(locl.matrix_bit(ircv),path,ir0,ir1);
    }
    break;
    default:
    {
	std::cerr << "TomographicInversion3d: illegal raycode '" << int(get_ray_type(isrc,ircv)) << "' detected.\n";
	com().abort(-1);
    }
    }
    if (iterbend < 0) {
	com().abort(1);
    }
    Point3d src3d = src(isrc);
    if (is_multiple(isrc, ircv)) {// add travel time in water column
	double pwater = 1.0/1.5;
	final_t += 2*bathyp->z(src3d.x(),src3d.y())*pwater;
    }
    if (only_forward){
	locl.syn_ttime(ircv) = final_t;
	locl.res_ttime(ircv) = obs_ttime(isrc)(ircv)-final_t;
    }else{
	locl.res_ttime(ircv) = obs_ttime(isrc)(ircv)-final_t;
    }
    
    if (printTransient || (printFinal && locl.global.last_step())){
        typedef boost::lock_guard<boost::mutex> os_guard;
        if (my_output_level >= 2){
            // There is  a memory/speed tradeoff here.
            // On the onw hand this printing is expenssive and having a guard here
            // can slow down things.
            // Anoher way to proceed would be to compute and store the curves
            // in paralel and then print them in one chunk.
            // This would use more memory
            os_guard g(locl.ray_os_mutex);
            std::ofstream out(indexed_filename("ray",locl.global.step(), isrc).c_str(), std::ofstream::app);
            out << "> ircv:" << ircv
		<< ", ray type: " << get_ray_type(isrc,ircv)
                << ", rcv in water: " << slowness_mesh().inWater(r)
		<< ", rcv in air: " << slowness_mesh().inAir(r)
                << '\n';
            printCurve(out, path, betasp);
            out.flush();
        }
    }
}

void
TomographicInversion3d::process_source(TomographicInversion3d::source_solver& locl) const {
    clock_t start_t=clock();
    locl.graph_time = 0;
    int isrc = locl.src;
    Array1d<Point3d> const& receivers =  rcv(isrc);
    if (is_verbose(0)){
	cerr << "\nstart isrc=" << isrc << " nrec=" << receivers.size() << '\n';
    }
    
    Point3d src3d = src(isrc);
    boost::shared_ptr<SlownessMesh3d> refr_mesh   (new SlownessMesh3d(slowness_mesh()));
    boost::shared_ptr<GraphSolver3d>  refr_solver (new GraphSolver3d(my_graph_solver, *refr_mesh));
    limit_range(*refr_solver, isrc);
    refr_solver->solve(src3d, ani);

    boost::shared_ptr<SlownessMesh3d> refl_mesh;
    boost::shared_ptr<GraphSolver3d>  refl_solver;
    // Optimization: maybe we could start dealing with refractions 
    // while we do this reflection preparation ?
    if (do_full_refl) {
	refl_mesh.reset(new SlownessMesh3d(slowness_mesh()));
	refl_solver.reset(new GraphSolver3d(my_graph_solver, *refl_mesh));
	refl_mesh->set(my_full_reflection_velocity_model);
    } else {
	refl_mesh   = refr_mesh;
	refl_solver = refr_solver;
    }
    if (reflp) {
      refl_solver->solve_refl(src3d, *reflp, ani);
    }
    locl.res_ttime.resize(receivers.size());
    locl.syn_ttime.resize(receivers.size());
    locl.matrix_bit.resize(receivers.size());
    locl.refraction_solver = refr_solver;
    locl.reflection_solver = refl_solver;

    {
	boost::thread_group workers;
	if (use_multithreading()) {
	    // add "my" thread
	    workers.create_thread(path_processor(locl, false));
	    while (!locl.all_receivers_processed()) {
		if (locl.global.available_threads() > 0) {
		    workers.create_thread(path_processor(locl, true));
		}
	    }
	    workers.join_all();
	} else {
	    path_processor(locl, false)();
	}
        
    }
    if (my_output_level >= 1){
	std::ofstream out(indexed_filename("tres",locl.global.step(), isrc).c_str(), std::ofstream::trunc);
	for (int ircv = 1; ircv <= receivers.size(); ++ircv) {
	    Point3d r = receivers(ircv);
	    out << r.x() << '\t' << r.y() << '\t' << r.z() << '\t' << locl.res_ttime(ircv) << '\n';
	}
    }
    if (is_verbose(0)){
	cerr << "\nfinish isrc=" << isrc << " nrec=" << receivers.size() << '\n';
    }
    locl.bend_time =(clock()-start_t) - locl.graph_time;
}

void
TomographicInversion3d::source_processor::operator()() {
    solver::thread_user tu(global);
    int isrc;
    while ((isrc = global.next_source()) >= 0) {
	source_solver solver(global, isrc);
	global.inversion.process_source(solver);
	global.add_time(solver.graph_time, solver.bend_time);
	global.inversion.res_ttime(isrc).swap(solver.res_ttime);
	global.inversion.syn_ttime(isrc).swap(solver.syn_ttime);
	global.my_matrix_bits[isrc-1].swap(solver.matrix_bit);
    }
}

void
TomographicInversion3d::path_processor::operator()() {
    static boost::mutex io;
    static std::string const own("free");
    static std::string const inherited("parent");

    if (parent.global.inversion.is_verbose(1)) {
	lock g(io);
	std::cerr << "assigning " << (my_own_thread ? own : inherited) << " thread " 
		  << boost::this_thread::get_id() << " to solving paths of source " 
		  << parent.src << '\n';
    }
    {
	solver::thread_user tu(parent.global, my_own_thread);
	int ircv;
	while ((ircv = parent.next_receiver()) >= 0) {
	    path_solver solver(parent, ircv);
	    solver();
	}
    }
    if (parent.global.inversion.is_verbose(1)) {
	lock g(io);
	std::cerr << "releasing " << (my_own_thread ? own : inherited) << " thread " 
		  << boost::this_thread::get_id() << " from solving paths of source " 
		  << parent.src << '\n';
    }
}
