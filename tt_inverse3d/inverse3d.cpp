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
#include <cassert>
#include <iostream>
#include <fstream>
#include <numeric>

#include "boost/thread.hpp"
#include "geom3d.hpp"
#include "inverse3d.hpp"
#include "sparse_rect.hpp"
#include "axb_solver.hpp"
#include "d_nr.hpp"

using std::cout;

TomographicInversion3d::TomographicInversion3d(boost::mpi::communicator& com,
                                               SlownessMesh3d& m, const char* datafn,
                                               int xorder, int yorder, int zorder, double crit_len,
                                               int nintp, double cg_tol, double br_tol) // Modified

:   my_com(com), my_graph_solver(m,xorder,yorder,zorder),
    betasp(1,0,nintp),
    cg_tolerance(cg_tol), br_tolerance(br_tol),
    reflp(), anid(), anie(),
    my_nb_noded(-1), do_full_refl(false), my_nb_data(-1),
    nnodev(slowness_mesh().nb_nodes()),
    nx(slowness_mesh().nx()),
    ny(slowness_mesh().ny()),
    nz(slowness_mesh().nz()),
    my_velocity_model(),
    my_full_reflection_velocity_model(),
    itermax_LSQR(nnodev*10), LSQR_ATOL(1e-3), 
    jumping(false),
    smooth_velocity(false), logscale_vel(false), my_apply_filter(false),
    wsv_min(0.0), wsv_max(0.0), dwsv(1.0), weight_s_v(0.0), corr_vel_p(0),
    smooth_depth(false), logscale_dep(false), wsd_min(0.0),
    wsd_max(0.0), dwsd(1.0), weight_s_d(0.0), corr_dep_p(0),
    smooth_anid(false), logscale_anid(false),
    smooth_anie(false), logscale_anie(false),
    wsad_min(0.0), wsad_max(0.0), dwsad(1.0), weight_s_ad(0.0), corr_anid_p(0),
    wsae_min(0.0), wsae_max(0.0), dwsae(1.0), weight_s_ae(0.0), corr_anie_p(0),
    damp_velocity(false), damp_depth(false), damp_anid(false), damp_anie(false),
    target_dv(0.0), target_dd(0.0), target_dad(0.0), target_dae(0.0),
    damping_is_fixed(false), weight_d_v(0.0), weight_d_d(0.0),
    weight_d_ad(0.0), weight_d_ae(0.0),
    damping_wt_p(0), dampingD_wt_p(0), dampingE_wt_p(0),
    do_squeezing(false), do_squeezingD(false), do_squeezingE(false),
    robust(false), crit_chi(1e3),

    A(), Tv(nnodev), Td(), Tad(), Tae(),
    my_rvz(), my_rdz(), my_radz(), my_raez(),
    my_rvx(), my_rvy(), my_rdx(), my_rdy(),
    my_radx(), my_rady(), my_raex(), my_raey(),

    refl_weight(1.0), target_chisq(0.0),
    phi(-1.0), ani(false),
    my_verbosity_level(0),
    printFinal(false), printTransient(false),
    out_mask(false), out_mask_refl(false), out_mask_anid(false), out_mask_anie(false),
    my_output_level(0),
    nb_threads(boost::thread::hardware_concurrency()),
    my_log_filename()
{
    if (nb_threads == 0) {
        nb_threads = 1;
    }
    if (crit_len>0) {
        my_graph_solver.refine_if_long(crit_len);
    }

    bathyp = new Interface3d(slowness_mesh());
    
    my_nb_data = read_file(datafn);
    A.resize(nb_data());
    
    const int max_np = int(3*pow(double(nnodev),1.0/3.0));
    
    if (!xflat()) {
        my_rvx.reset(new sparse_matrix(nnodev));
    }
    if (!yflat()) {
        my_rvy.reset(new sparse_matrix(nnodev));
    }
    my_rvz.resize(nnodev);

    data_vec.resize(nb_data());
    row_index.resize(nb_data());
    total_data_vec.reserve(2*nb_data()+4*nnodev+3*max_np);
    r_dt_vec.resize(nb_data());
    int idata=1;
    for (int i=1; i<=obs_dt.size(); i++){
        for (int j=1; j<=obs_dt(i).size(); j++){
            r_dt_vec(idata) = 1.0/obs_dt(i)(j);
            idata++;
        }
    }

    mvscale.resize(nnodev);
    dmodel_total.resize(nnodev);
    nnode_total = nnodev;
    nnode_subtotal = nnode_total;
    nnoded = 0;
    nnodee = 0;
}

bool
TomographicInversion3d::need_reflections() const {
    typedef std::vector<std::vector<ray_type> > src_ray_types;
    for(src_ray_types::const_iterator sit = my_ray_types.begin(); sit != my_ray_types.end(); ++sit) {
        typedef std::vector<ray_type> rcv_ray_types;
        for(rcv_ray_types::const_iterator rit = sit->begin(); rit != sit->end(); ++ rit) {
            if (is_reflection(*rit)) {
                return true;
            }
        }
    }
    return false;
}

std::string
TomographicInversion3d::indexed_filename(std::string base, int idx1, int idx2) const {
    std::ostringstream fmt;
    fmt << my_output_prefix << '.' << base << '.' << idx1 << '.' << idx2;
    return fmt.str();
}

TomographicInversion3d::ostream_ptr
TomographicInversion3d::indexed_trace(std::string base, int idx1, int idx2, bool cond) const {
    std::ostream* out = cond ? new std::ofstream(indexed_filename(base,idx1,idx2).c_str(), std::ofstream::app) : 0;
    return ostream_ptr(out);
}

std::istream&
operator>>(std::istream& in, TomographicInversion3d::ray_type& tp )  {
    int n;
    typedef TomographicInversion3d inv;
    in >> n;
    if (!in.good()) {
        tp = static_cast<inv::ray_type>(-1);
    }
    switch(static_cast<inv::ray_type>(n)) {
    case inv::RAY_REFRACTION:
    case inv::RAY_REFLECTION:
    case inv::RAY_MULTIPLE_REFRACTION:
    case inv::RAY_MULTIPLE_REFLECTION:
        tp = static_cast<inv::ray_type>(n);
        break;
    default:
        tp = static_cast<inv::ray_type>(-1);
    }
    return in;
}

int
TomographicInversion3d::read_file(std::string const& datafn) {
    ifstream in(datafn.c_str());
    if (!in){
        cerr << "TomographicInversion3d::read_file: cannot open '" << datafn << "'\n";
        com().abort(1);
    }
    
    int iline=0;
    
    string message;
    int nsrc;
    in >> nsrc;
    ++iline;
    if (nsrc<=0) {
        cerr << "TomographicInversion3d::invalid nsrc in '" << datafn<< "'\n";
        com().abort(1);
    }
    src.resize(nsrc); 
    rcv.resize(nsrc);
    my_ray_types.resize(nsrc);
    obs_ttime.resize(nsrc);
    shot_point_list.resize(nsrc);
    obs_dt.resize(nsrc);
    syn_ttime.resize(nsrc);
    res_ttime.resize(nsrc);
    
    int isrc=0;
    while(in){
        char flag;
        double x, y, z;
        int nrcv;

        in >> flag >> x >> y >> z >> nrcv; 
        getline( in, message ); /*ALOUREIRO Ignore comments after data*/
        cout<< "Current instrument is "<<message<<"\n"; /*ALOUREIRO log what data is being read*/
        ++iline;
        if (flag!='s'){
            cerr << "TomographicInversion3d::bad input (s) at "
                 << datafn << ":" << iline << '\n';
            com().abort(1);
        }
        isrc++;
        src(isrc).set(x,y,z);
        
        rcv(isrc).resize(nrcv);
        my_ray_types[isrc-1].resize(nrcv);
        obs_ttime(isrc).resize(nrcv); 
        shot_point_list(isrc).resize(nrcv); 
        obs_dt(isrc).resize(nrcv);
        for (int ircv=1; ircv<=nrcv; ircv++){
            ray_type n;
            double t_val, dt_val;
            in >> flag >> x >> y >> z >> n >> t_val >> dt_val;
            getline( in, message );
            ++iline;
            if (flag!='r'){
                cerr << "TomographicInversion3d::bad input (r) at "
                     << datafn << ":" << iline << '\n';
                com().abort(1);
            }
            if (dt_val < 1e-6){
                cerr << "TomographicInversion3d::bad input (zero obs_dt) at "
                     << datafn << ":" << iline << '\n';
                com().abort(1);
            }
            rcv(isrc)(ircv).set(x,y,z);
            set_ray_type(isrc, ircv, n);
            obs_ttime(isrc)(ircv) = t_val;
            obs_dt(isrc)(ircv) = dt_val;
            shot_point_list(isrc)(ircv)=message;
            if (ircv==1){
                cout<< "First Shot Point is: "<<message<<"\n";
            }
        }
        if (isrc==nsrc) break;
    }
    if (isrc != nsrc) {
        cerr << "TomographicInversion3d::mismatch in nsrc in '"
             << datafn << "'\n";
        com().abort(1);
    }

    int data_sz = 0;
    for (int isrc=1; isrc<=nsrc; isrc++){
        data_sz += rcv(isrc).size();
    }
    return data_sz;
    
}

void TomographicInversion3d::doRobust(double a)
{
    if (a>0){
        robust = true;
        crit_chi = a;
    }
}

void TomographicInversion3d::removeOutliers()
{
    Array1d<double> Adm(ndata_valid);
    SparseRectangular sparseA(A,[&](int i) { return bool(row_index[i-1]);},identity,nnode_total);
    sparseA.Ax(dmodel_total,Adm);
    std::vector<int> row_index2(row_index.size(), 0);

    int ndata_valid2=0;
    for (int i=1; i<=row_index.size(); i++){
        int j=row_index[i-1];
        if (j>0){
            double res=Adm(j)-data_vec(i);
            if (abs(res)<=crit_chi){
                row_index2[i-1] = ++ndata_valid2;
            }
        }
    }

    // redefine data node vector and related stats
    ndata_valid = ndata_valid2;
    row_index = row_index2;
    rms_tres[0]=rms_tres[1]=rms_tres[2]=rms_tres[3]=0;
    init_chi[0]=init_chi[1]=init_chi[2]=init_chi[3]=0;
    ndata_in[0]=ndata_in[1]=ndata_in[2]=ndata_in[3]=0;
    int idata=1;
    for (int isrc=1; isrc<=src.size(); isrc++){
        for (int ircv=1; ircv<=rcv(isrc).size(); ircv++){
            if (row_index[idata-1]>0){
                double res=res_ttime(isrc)(ircv);
                ray_type icode = get_ray_type(isrc,ircv);
                double res2=res*r_dt_vec(idata);
                double res22=res2*res2;
                rms_tres[icode] += res*res;
                init_chi[icode] += res22;
                ++ndata_in[icode];
            }
            idata++;
        }
    }
    rms_tres_total = sqrt((rms_tres[0]+rms_tres[1]+rms_tres[2]+rms_tres[3])/ndata_valid);
    init_chi_total = (init_chi[0]+init_chi[1]+init_chi[2]+init_chi[3])/ndata_valid;
    cout << "rms =" << rms_tres_total <<"\n"; //ALOUREIRO
    cout << "chi2 =" << init_chi_total <<"\n"; //ALOUREIRO
    for (int i=0; i<=3; i++){
        rms_tres[i] = ndata_in[i]>0 ? sqrt(rms_tres[i]/ndata_in[i]) : 0.0;
        init_chi[i] = ndata_in[i]>0 ? init_chi[i]/ndata_in[i] : 0.0;
    }
}

void TomographicInversion3d::setLSQR_TOL(double a)
{
    if (a<=0){
        cerr << "TomographicInversion3d::setLSQR_TOL - invalid TOL (ignored)\n";
        return;
    }
    LSQR_ATOL = a;
}

void
TomographicInversion3d::set_log(std::string const& fn) {
    my_log_filename = fn;
}

bool
TomographicInversion3d::pick_one_process() const {
    static int round_robin = 0;
    return round_robin++ % com().size() == com().rank();
}

TomographicInversion3d::ostream_ptr
TomographicInversion3d::log_stream() const {
    ostream_ptr log;
    if (!my_log_filename.empty()) {
        if (com().rank() == 0) {
            log.reset(new std::ofstream(my_log_filename.c_str(), std::ofstream::app));
            if (!*log) {
                std::cerr << "could not open log file '" << my_log_filename << "'\n";
            }
        }
    }
    return log;
}

bool
TomographicInversion3d::is_verbose(int level, bool any) const {
    bool high_enough = my_verbosity_level >= level;
    return ( high_enough && 
             (any || pick_one_process()));
}

void TomographicInversion3d::fixing(bool a, bool b, bool c, bool d)
{
    fv=a;
    fd=b;
    fad=c;
    fae=d;
}

void TomographicInversion3d::outStepwise(int i)
{
    printTransient = true;
    my_output_level = i;
}

void TomographicInversion3d::outFinal(int i)
{
    printFinal = true;
    my_output_level = i;
}

void TomographicInversion3d::targetChisq(double c)
{
    target_chisq = c;
}

void TomographicInversion3d::outMask()
{
    out_mask = true;
}

void TomographicInversion3d::outMaskRefl()
{
    out_mask_refl = true;
}

void TomographicInversion3d::outMaskAnid()
{
    out_mask_anid = true;
}

void TomographicInversion3d::outMaskAnie()
{
    out_mask_anie = true;
}

void TomographicInversion3d::onlyForward(bool fit){ only_forward = true; first_it = fit;}

void
TomographicInversion3d::add_reflections(boost::shared_ptr<Interface3d> const& intfp)
{
    reflp = intfp;
    my_nb_noded = reflp->numNodes();
    if (!xflat()) {
        my_rdx.reset(new sparse_matrix(nb_noded()));
    }
    if (!yflat()) {
        my_rdy.reset(new sparse_matrix(nb_noded()));
    }
    my_rdz.resize(nb_noded());
    Td.resize(nb_noded());
    modeld.resize(nb_noded());
    mdscale.resize(nb_noded());
    dmodel_total.resize(nnodev+nb_noded());
    nnode_total = nnodev+nb_noded();
    nnode_subtotal = nnode_total;
    tmp_nodedc.resize(nb_noded());
    tmp_nodedr.resize(nb_noded());
}

void TomographicInversion3d::doFullRefl()
{
    do_full_refl = true;
    my_graph_solver.do_refl_downward();
}

void TomographicInversion3d::setReflWeight(double x)
{
    if (x>0){
        refl_weight = x;
    }else{
        cerr << "TomographicInversion3d::setReflWeight - non-positive value ignored.\n";
    }
}

void TomographicInversion3d::tiltAngle(double x)
{
  if(x>=0 && x<acos(-1.0)){
    phi = x;
  }else{
    cerr << "TomographicInversion3d::tiltAngle - tilt must be 0 or positive and smaller than Pi. Check tilt definition in the user guide.\n";
  }
}

void TomographicInversion3d::add_anisotropy(boost::shared_ptr<AnisotropyMesh3d> const& anidp,
					    boost::shared_ptr<AnisotropyMesh3d> const& aniep)
{
    ani = true;

    anid = anidp;
    anie = aniep;

    nnoded = anid->nb_nodes();
    nnodee = anie->nb_nodes();

    if (!xflat()) {
        my_radx.reset(new sparse_matrix(nnoded));
    }
    if (!yflat()) {
        my_rady.reset(new sparse_matrix(nnoded));
    }
    my_radz.resize(nnoded);

    if (!xflat()) {
        my_raex.reset(new sparse_matrix(nnodee));
    }
    if (!yflat()) {
        my_raey.reset(new sparse_matrix(nnodee));
    }
    my_raez.resize(nnodee);

    Tad.resize(nnoded);
    Tae.resize(nnodee);

    modelad.resize(nnoded);
    modelae.resize(nnodee);
    madscale.resize(nnoded);
    maescale.resize(nnodee);

    dmodel_total.resize(nnode_total+nnoded+nnodee);
    nnode_subtotal = nnode_total;
    nnode_total = nnode_total+nnoded+nnodee;

    tmp_nodeadc.resize(nnoded);
    tmp_nodeadr.resize(nnoded);
    tmp_nodeaec.resize(nnodee);
    tmp_nodeaer.resize(nnodee);
}

void TomographicInversion3d::usingVh()
{
    useVh = true;
}

void TomographicInversion3d::SmoothVelocity(const char* fn,
                                            double start, double end, double d,
                                            bool scale)
{
    smooth_velocity=true;
    wsv_min=start; wsv_max=end; dwsv=d;
    logscale_vel=scale;
    corr_vel_p = new CorrelationLength3d(fn);
}

void TomographicInversion3d::applyFilter(const char *fn)
{
    my_apply_filter=true;
    if (fn[0] != '\0'){
        uboundp = new Interface3d(fn);
    }else{
        uboundp = new Interface3d(slowness_mesh()); // use bathymetry as upperbound
    }
}

void TomographicInversion3d::SmoothDepth(const char* fn,
                                         double start, double end, double d,
                                         bool scale)
{
    SmoothDepth(start,end,d,scale);
    corr_dep_p = new CorrelationLength2d(fn);
}

void TomographicInversion3d::SmoothDepth(double start, double end, double d,
                                         bool scale)
{
    smooth_depth=true;
    wsd_min=start; wsd_max=end; dwsd=d;
    logscale_dep=scale;
}

void TomographicInversion3d::SmoothAniD(const char* fn,
                                         double start, double end, double d,
                                         bool scale)
{
    smooth_anid=true;
    wsad_min=start; wsad_max=end; dwsad=d;
    logscale_anid=scale;
    corr_anid_p = new CorrelationLength3d(fn);
}

void TomographicInversion3d::SmoothAniE(const char* fn,
					double start, double end, double d,
					bool scale)
{
    smooth_anie=true;
    wsae_min=start; wsae_max=end; dwsae=d;
    logscale_anie=scale;
    corr_anie_p = new CorrelationLength3d(fn);
}

void TomographicInversion3d::DampVelocity(double a)
{
    damp_velocity = true;
    target_dv = a/100.0; // % -> fraction
}

void TomographicInversion3d::DampDepth(double a)
{
    damp_depth = true;
    target_dd = a/100.0; // % -> fraction 
}

void TomographicInversion3d::DampAniD(double a)
{
    damp_anid = true;
    target_dad = a/100.0; // % -> fraction
}

void TomographicInversion3d::DampAniE(double a)
{
    damp_anie = true;
    target_dae = a/100.0; // % -> fraction
}

void TomographicInversion3d::FixDamping(double v, double d, double ad, double ae)
{
    damping_is_fixed = true;
    damp_velocity = true;
    damp_depth = true;
    damp_anid = true;
    damp_anie = true;
    weight_d_v = v;
    weight_d_d = d;
    weight_d_ad = ad;
    weight_d_ae = ae;
}

void TomographicInversion3d::Squeezing(const char* fn)
{
    damping_wt_p = new DampingWeight3d(fn);
    do_squeezing = true;
}

void TomographicInversion3d::SqueezingD(const char* fn)
{
    dampingD_wt_p = new DampingWeight3d(fn);
    do_squeezingD = true;
}

void TomographicInversion3d::SqueezingE(const char* fn)
{
    dampingE_wt_p = new DampingWeight3d(fn);
    do_squeezingE = true;
}

void TomographicInversion3d::doJumping()
{
    jumping = true;
}

namespace {
    void
    clear_sparse_mat(sparse_matrix& sm) {
        for(auto& s : sm) {
            s.clear();
        }
    }
    void
    clear_sparse_mat(std::unique_ptr<sparse_matrix>& sm) {
        if (sm) {
            clear_sparse_mat(*sm);
        }
    }
}

void TomographicInversion3d::reset_kernel()
{
    clear_sparse_mat(A);
    clear_sparse_mat(my_rvx);
    clear_sparse_mat(my_rvy);
    clear_sparse_mat(my_rvz);
    clear_sparse_mat(Tv);
    if (reflp){
        clear_sparse_mat(my_rdx);
        clear_sparse_mat(my_rdy);
        clear_sparse_mat(Td);
    }
    if (ani){
	clear_sparse_mat(my_radx);
	clear_sparse_mat(my_rady);
	clear_sparse_mat(my_radz);
	clear_sparse_mat(Tad);
	clear_sparse_mat(my_raex);
	clear_sparse_mat(my_raey);
	clear_sparse_mat(my_raez);
	clear_sparse_mat(Tae);
    }
}

namespace {
    template<typename T>
    bool is_nan(T const& v) {
        return v != v;
    }
}
void
TomographicInversion3d::add_kernel(map<int,double>& A_i, const Array1d<Point3d>& cur_path) const
{
    int np = cur_path.size();
    int nintp = betasp.numIntp();
    
    Array1d<const Point3d*> pp;
    pp.reserve(int(3*pow(double(nnodev),1.0/3.0)));
    makeBSpoints(cur_path,pp);
    
    double path_len=0.0;
    Index3d guess_index = slowness_mesh().nodeIndex(slowness_mesh().nearest(*pp(1)));
    Array1d<Point3d> Q(betasp.numIntp());
    for (int i=1; i<=np+1; i++){
        int j1=i;
        int j2=i+1;
        int j3=i+2;
        int j4=i+3;
        betasp.interpolate(*pp(j1),*pp(j2),*pp(j3),*pp(j4),Q);
        
        for (int j=2; j<=nintp; j++){
            Point3d midp = 0.5*(Q(j-1)+Q(j));
            double dist= Q(j).distance(Q(j-1));
            path_len+=dist;
            
            typedef SlownessMesh3d::opt_idx opt_idx;
            opt_idx jLFU, jRFU, jLBU, jLFL, jLBL, jRFL, jRBU, jRBL; // jXYZ: X=L(eft):R(ight), Y=F(ront):B(ack), Z=U(pper):L(ower).
            Point3d r;
            int icell = slowness_mesh().locate_in_cell(midp,guess_index,
						       jLFU,jRFU,jLBU,jLFL,jLBL,jRFL,jRBU,jRBL,r);
            Point3d rr = Point3d(1,1,1) - r;
	    
	    // Declare anisotropy variables
	    double velmidp, cos_theta, theta, gamma, kanid, kanie, delta, epsilon;
	    Point3d ray, k;
	    typedef AnisotropyMesh3d::opt_idx opt_idxd, opt_idxe;
	    opt_idxd dLFU, dRFU, dLBU, dLFL, dLBL, dRFL, dRBU, dRBL;
	    opt_idxe eLFU, eRFU, eLBU, eLFL, eLBL, eRFL, eRBU, eRBL;
	    int dcell = 0;
	    int ecell = 0;
	    Point3d rd(0,0,0), re(0,0,0);
	    Point3d rrd = Point3d(1,1,1)-rd;
	    Point3d rre = Point3d(1,1,1)-re;
	    if (ani){ // ani = true when anisotropy is considered
	      Index3d dguess_index = anid->nodeIndex(anid->nearest(*pp(1)));
	      Index3d eguess_index = anie->nodeIndex(anie->nearest(*pp(1)));
	      dcell = anid->locate_in_cell(midp,dguess_index,
					   dLFU,dRFU,dLBU,dLFL,dLBL,dRFL,dRBU,dRBL,rd);
	      ecell = anie->locate_in_cell(midp,eguess_index,
					   eLFU,eRFU,eLBU,eLFL,eLBL,eRFL,eRBU,eRBL,re);
	      if (dist != 0.0){
		// Call to SlownessMesh3d::at to get velocity at midp
		velmidp = 1 / slowness_mesh().at(midp);
		// Calculation of angle theta between ray path segment and isotropic axis	
		k.set(0.0,0.0,1.0);
		ray = Q(j)-Q(j-1);
		cos_theta = ray.inner_product(k)/ray.norm();
		theta = acos(cos_theta);
		if (theta > asin(1.0)){
		  theta = acos(-1.0) - theta;
		}
		if(theta >= phi){
		  gamma = theta - phi;
		}else if(theta < phi){
		  gamma = phi - theta;
		}
		// phi passed from tt_inverse3d.cpp
		// Calls to AnisotropyMesh3d::at to get delta and epsilon at midp
		delta = anid->at(midp);
		if(useVh){
		    double vh = anie->at(midp);
		    epsilon = (vh/velmidp)-1;
		}else{
		    epsilon = anie->at(midp);
		}
		// delta kernel
		// Computation of kanid = dist*sin²(theta-phi)*cos²(theta-phi)*velocity
		double ps = (1+delta*sin(gamma)*sin(gamma)*cos(gamma)*cos(gamma)
			     +epsilon*sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma));

		kanid = -dist*sin(gamma)*sin(gamma)*cos(gamma)*cos(gamma)/(velmidp*ps*ps);
		// epsilon kernel
		// Computation of kanie = dist*sin⁴(theta-phi)*velocity
		if(useVh){
		    kanie = -dist*sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma)/(velmidp*velmidp*ps*ps);
		}else{
		    kanie = -dist*sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma)/(velmidp*ps*ps);
		}
		// anisotropic velocity kernel
		// Computation of dist=dist*(1+delta*sin²(theta-phi)*cos²(theta-phi)+epsilon*sin⁴(theta-phi))
		dist = dist/(1+delta*sin(gamma)*sin(gamma)*cos(gamma)*cos(gamma)
			     +epsilon*sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma));
	      }else{
		dist = 0.0;
		kanid = 0.0;
		kanie = 0.0;
	      }
	    }
		    
            if (icell>0 && !fv){
                // depending on the geometry, the number of nodes may vary
                if (jLFU) {
                    A_i[*jLFU] += rr.x()*rr.y()*rr.z()*dist;
                    assert(!is_nan(A_i[*jLFU]));
                }
                if (jRFU) {
                    A_i[*jRFU] += r.x()*rr.y()*rr.z()*dist;
                    assert(!is_nan(A_i[*jRFU]));
                }
                if (jLBU) {
                    A_i[*jLBU] += rr.x()*r.y()*rr.z()*dist;
                    assert(!is_nan(A_i[*jLBU]));
                }
                if (jLFL) {
                    A_i[*jLFL] += rr.x()*rr.y()*r.z()*dist;
                    assert(!is_nan(A_i[*jLFL]));
                }
                if (jLBL) {
                    A_i[*jLBL] += rr.x()*r.y()*r.z()*dist;
                    assert(!is_nan(A_i[*jLBL]));
                }
                if (jRFL) {
                    A_i[*jRFL] += r.x()*rr.y()*r.z()*dist;
                    assert(!is_nan(A_i[*jRFL]));
                }
                if (jRBU) {
                    A_i[*jRBU] += r.x()*r.y()*rr.z()*dist;
                    assert(!is_nan(A_i[*jRBU]));
                }
                if (jRBL) {
                    A_i[*jRBL] += r.x()*r.y()*r.z()*dist;
                    assert(!is_nan(A_i[*jRBL]));
                }
            }
	    if(ani && dcell > 0 && !fad){
		if(dLFU){
		    A_i[*dLFU+nnode_subtotal] += rrd.x()*rrd.y()*rrd.z()*kanid;
		    assert(!is_nan(A_i[*dLFU+nnode_subtotal]));
		}
		if(dRFU){
		    A_i[*dRFU+nnode_subtotal] += rd.x()*rrd.y()*rrd.z()*kanid;
		    assert(!is_nan(A_i[*dRFU+nnode_subtotal]));
		}
		if(dLBU){
		    A_i[*dLBU+nnode_subtotal] += rrd.x()*rd.y()*rrd.z()*kanid;
		    assert(!is_nan(A_i[*dLBU+nnode_subtotal]));
		}
		if(dLFL){
		    A_i[*dLFL+nnode_subtotal] += rrd.x()*rrd.y()*rd.z()*kanid;
		    assert(!is_nan(A_i[*dLFL+nnode_subtotal]));
		}
		if(dLBL){
		    A_i[*dLBL+nnode_subtotal] += rrd.x()*rd.y()*rd.z()*kanid;
		    assert(!is_nan(A_i[*dLBL+nnode_subtotal]));
		}
		if(dRFL){
		    A_i[*dRFL+nnode_subtotal] += rd.x()*rrd.y()*rd.z()*kanid;
		    assert(!is_nan(A_i[*dRFL+nnode_subtotal]));
		}
		if(dRBU){
		    A_i[*dRBU+nnode_subtotal] += rd.x()*rd.y()*rrd.z()*kanid;
		    assert(!is_nan(A_i[*dRBU+nnode_subtotal]));
		}
		if(dRBL){
		    A_i[*dRBL+nnode_subtotal] += rd.x()*rd.y()*rd.z()*kanid;
		    assert(!is_nan(A_i[*dRBL+nnode_subtotal]));
		}
	    }
	    if(ani && ecell > 0 && !fae){
		if(eLFU){
		    A_i[*eLFU+nnode_subtotal+nnoded] += rre.x()*rre.y()*rre.z()*kanie;
		    assert(!is_nan(A_i[*eLFU+nnode_subtotal+nnoded]));
		}
		if(eRFU){
		    A_i[*eRFU+nnode_subtotal+nnoded] += re.x()*rre.y()*rre.z()*kanie;
		    assert(!is_nan(A_i[*eRFU+nnode_subtotal+nnoded]));
		}
		if(eLBU){
		    A_i[*eLBU+nnode_subtotal+nnoded] += rre.x()*re.y()*rre.z()*kanie;
		    assert(!is_nan(A_i[*eLBU+nnode_subtotal+nnoded]));
		}
		if(eLFL){
		    A_i[*eLFL+nnode_subtotal+nnoded] += rre.x()*rre.y()*re.z()*kanie;
		    assert(!is_nan(A_i[*eLFL+nnode_subtotal+nnoded]));
		}
		if(eLBL){
		    A_i[*eLBL+nnode_subtotal+nnoded] += rre.x()*re.y()*re.z()*kanie;
		    assert(!is_nan(A_i[*eLBL+nnode_subtotal+nnoded]));
		}
		if(eRFL){
		    A_i[*eRFL+nnode_subtotal+nnoded] += re.x()*rre.y()*re.z()*kanie;
		    assert(!is_nan(A_i[*eRFL+nnode_subtotal+nnoded]));
		}
		if(eRBU){
		    A_i[*eRBU+nnode_subtotal+nnoded] += re.x()*re.y()*rre.z()*kanie;
		    assert(!is_nan(A_i[*eRBU+nnode_subtotal+nnoded]));
		}
		if(eRBL){
		    A_i[*eRBL+nnode_subtotal+nnoded] += re.x()*re.y()*re.z()*kanie;
		    assert(!is_nan(A_i[*eRBL+nnode_subtotal+nnoded]));
		}
	    }
        }
    }
}

void TomographicInversion3d::add_kernel_refl(map<int,double>& A_i, const Array1d<Point3d>& cur_path,
                                             int ir0, int ir1) const
{
    typedef Interface3d::opt_idx opt_idx;
    opt_idx jLF, jRF, jLB, jRB; // jRB is dispensable as its coordinates are given by JRF and JLB.
    reflp->locate_in_plane(cur_path(ir0).x(),cur_path(ir0).y(),jLF,jRF,jLB,jRB);

    if (!jLF){
        cerr << "TomographicInversion3d::add_kernel_refl: bottoming point out of bounds\n"
             << "x=" << cur_path(ir0).x() << " y=" << cur_path(ir0).y() << '\n';
        return; // ignore this ray info
    }

    Point3d a = cur_path(ir0-1)-cur_path(ir0);
    a /= a.norm();
    Point3d b = cur_path(ir1+1)-cur_path(ir1);
    b /= b.norm();
    Point3d c = a+b; 
    c /= c.norm();
    double cos_theta = a.inner_product(c);
    double alt_cos_theta = b.inner_product(c);

    double x1  = reflp->x(*jLF);
    double x2  = reflp->xflat() ? x1 : reflp->x(*jRF);
    double y1  = reflp->y(*jLF);
    double y2  = reflp->yflat() ? y1 : reflp->y(*jLB);
    double x   = cur_path(ir0).x();
    double y   = cur_path(ir0).y();
    double z   = reflp->z(x,y);
    double z1  = reflp->z(x1,y1);
    double cos_alpha; 
    if (reflp->xflat()) {
        double z2y = reflp->z(x1,y2);
        Point3d q(0.,y2-y1,z2y-z1);
        Point3d ey(0.0,1.0,0.0);
        cos_alpha = q.inner_product(ey)/q.norm();
    } else if (reflp->yflat()) {
        double z2x = reflp->z(x2,y1);
        Point3d p(x2-x1,0.0,z2x-z1);
        Point3d ex(1.0,0.0,0.0);
        cos_alpha = p.inner_product(ex)/p.norm();
    } else {
        double z2x = reflp->z(x2,y1);
        double z2y = reflp->z(x1,y2);
        Point3d p(x2-x1,0.,z2x-z1); // note: this doesn't fail even when x==x1 (or x2)
        Point3d q(0.,y2-y1,z2y-z1);
        Point3d n = p.outer_product(q);
        Point3d ez(0.0,0.0,1.0);
        
        cos_alpha = n.inner_product(ez)/n.norm();
    }

    double p1 = slowness_mesh().at(cur_path(ir0));
    double com_factor = 2.0*cos_theta*cos_alpha*p1;
    double rx = reflp->xflat() ? 0 : (x-x1)/(x2-x1);
    double ry = reflp->yflat() ? 0 : (y-y1)/(y2-y1);

    if(rx<0. && fabs(rx)>1e-3){
        cerr<<"too big negative rx\n";
    }
    if(ry<0. && fabs(ry)>1e-3){
        cerr<<"too big negative ry\n";
    }

    if(fabs(rx)<1e-3){
        rx=0.;
    }
    if(fabs(ry)<1e-3){
        ry=0.;
    }
    if (!fd){
	if (jLF) {
	    A_i[nnodev+*jLF] = com_factor*(1-rx)*(1-ry);
	    assert(!is_nan(A_i[nnodev+*jLF]));
	}
	if (jRF) {
	    A_i[nnodev+*jRF] = com_factor*rx*(1-ry);
	    assert(!is_nan(A_i[nnodev+*jRF]));
	}
	if (jLB) {
	    A_i[nnodev+*jLB] = com_factor*(1-rx)*ry;
	    assert(!is_nan(A_i[nnodev+*jLB]));
	}
	if (jRB) {
	    A_i[nnodev+*jRB] = com_factor*rx*ry;
	    assert(!is_nan(A_i[nnodev+*jRB]));
	}
    }

    add_kernel(A_i,cur_path);
}

void
TomographicInversion3d::calc_averaging_matrix()
{
    for (int i=1; i<=nx; i++){
        for (int j=1; j<=ny; j++){
            for (int k=1; k<=nz; k++){
                int anode = slowness_mesh().nodeIndex(i,j,k);
                Point3d p = slowness_mesh().nodePos(anode);
                double Lx, Ly, Lz;
                corr_vel_p->at(p,Lx,Ly,Lz);
                if (!xflat()) {
                    double Lx2 = Lx*Lx;
                    double beta_sum = 0.0;
                    for (int ii=1; ii<=nx; ii++){
                        int bnode = slowness_mesh().nodeIndex(ii,j,k);
                        if (bnode!=anode){
                            double dx = slowness_mesh().nodePos(bnode).x()-p.x();
                            if (abs(dx)<=Lx){
                                double dxL2 = dx*dx/Lx2;
                                double beta = exp(-dxL2);
                                rvx(anode)[bnode] = beta*mvscale(bnode);
                                beta_sum += beta;
                            }
                        }
                    }
                    rvx(anode)[anode] = -beta_sum*mvscale(anode);
                }
                if (!yflat()) {
                    double Ly2 = Ly*Ly;
                    double beta_sum = 0.0;
                    for (int jj=1; jj<=ny; jj++){
                        int bnode = slowness_mesh().nodeIndex(i,jj,k);
                        if (bnode!=anode){
                            double dy = slowness_mesh().nodePos(bnode).y()-p.y();
                            if (abs(dy)<=Ly){
                                double dyL2 = dy*dy/Ly2;
                                double beta = exp(-dyL2);
                                rvy(anode)[bnode] = beta*mvscale(bnode);
                                beta_sum += beta;
                            }
                        }
                    }
                    rvy(anode)[anode] = -beta_sum*mvscale(anode);
                }
                {
                    double Lz2 = Lz*Lz;
                    double beta_sum = 0.0;
                    for (int kk=1; kk<=nz; kk++){
                        int bnode = slowness_mesh().nodeIndex(i,j,kk);
                        if (bnode!=anode){
                            double dz = slowness_mesh().nodePos(bnode).z()-p.z();
                            if (abs(dz)<=Lz){
                                double dzL2 = dz*dz/Lz2;
                                double beta = exp(-dzL2);
                                rvz(anode)[bnode] = beta*mvscale(bnode);
                                beta_sum += beta;
                            }
                        }
                    }
                    rvz(anode)[anode] = -beta_sum*mvscale(anode);
                }
            }
        }
    }
}

void
TomographicInversion3d::calc_anid_averaging_matrix()
{
    for (int i=1; i<=anid->nx(); i++){
        for (int j=1; j<=anid->ny(); j++){
            for (int k=1; k<=anid->nz(); k++){
                int anode = anid->nodeIndex(i,j,k);
                Point3d p = anid->nodePos(anode);
                double Lx, Ly, Lz;
                corr_anid_p->at(p,Lx,Ly,Lz);
                if (!xflat()) {
                    double Lx2 = Lx*Lx;
                    double beta_sum = 0.0;
                    for (int ii=1; ii<=anid->nx(); ii++){
                        int bnode = anid->nodeIndex(ii,j,k);
                        if (bnode!=anode){
                            double dx = anid->nodePos(bnode).x()-p.x();
                            if (abs(dx)<=Lx){
                                double dxL2 = dx*dx/Lx2;
                                double beta = exp(-dxL2);
                                radx(anode)[bnode] = beta*madscale(bnode);
                                beta_sum += beta;
                            }
                        }
                    }
                    radx(anode)[anode] = -beta_sum*madscale(anode);
                }
                if (!yflat()) {
                    double Ly2 = Ly*Ly;
                    double beta_sum = 0.0;
                    for (int jj=1; jj<=anid->ny(); jj++){
                        int bnode = anid->nodeIndex(i,jj,k);
                        if (bnode!=anode){
                            double dy = anid->nodePos(bnode).y()-p.y();
                            if (abs(dy)<=Ly){
                                double dyL2 = dy*dy/Ly2;
                                double beta = exp(-dyL2);
                                rady(anode)[bnode] = beta*madscale(bnode);
                                beta_sum += beta;
                            }
                        }
                    }
                    rady(anode)[anode] = -beta_sum*madscale(anode);
                }
                {
                    double Lz2 = Lz*Lz;
                    double beta_sum = 0.0;
                    for (int kk=1; kk<=anid->nz(); kk++){
                        int bnode = anid->nodeIndex(i,j,kk);
                        if (bnode!=anode){
                            double dz = anid->nodePos(bnode).z()-p.z();
                            if (abs(dz)<=Lz){
                                double dzL2 = dz*dz/Lz2;
                                double beta = exp(-dzL2);
                                radz(anode)[bnode] = beta*madscale(bnode);
                                beta_sum += beta;
                            }
                        }
                    }
                    radz(anode)[anode] = -beta_sum*madscale(anode);
                }
            }
        }
    }
}

void
TomographicInversion3d::calc_anie_averaging_matrix()
{
    for (int i=1; i<=anie->nx(); i++){
        for (int j=1; j<=anie->ny(); j++){
            for (int k=1; k<=anie->nz(); k++){
                int anode = anie->nodeIndex(i,j,k);
                Point3d p = anie->nodePos(anode);
                double Lx, Ly, Lz;
                corr_anie_p->at(p,Lx,Ly,Lz);
                if (!xflat()) {
                    double Lx2 = Lx*Lx;
                    double beta_sum = 0.0;
                    for (int ii=1; ii<=anie->nx(); ii++){
                        int bnode = anie->nodeIndex(ii,j,k);
                        if (bnode!=anode){
                            double dx = anie->nodePos(bnode).x()-p.x();
                            if (abs(dx)<=Lx){
                                double dxL2 = dx*dx/Lx2;
                                double beta = exp(-dxL2);
                                raex(anode)[bnode] = beta*maescale(bnode);
                                beta_sum += beta;
                            }
                        }
                    }
                    raex(anode)[anode] = -beta_sum*maescale(anode);
                }
                if (!yflat()) {
                    double Ly2 = Ly*Ly;
                    double beta_sum = 0.0;
                    for (int jj=1; jj<=anie->ny(); jj++){
                        int bnode = anie->nodeIndex(i,jj,k);
                        if (bnode!=anode){
                            double dy = anie->nodePos(bnode).y()-p.y();
                            if (abs(dy)<=Ly){
                                double dyL2 = dy*dy/Ly2;
                                double beta = exp(-dyL2);
                                raey(anode)[bnode] = beta*maescale(bnode);
                                beta_sum += beta;
                            }
                        }
                    }
                    raey(anode)[anode] = -beta_sum*maescale(anode);
                }
                {
                    double Lz2 = Lz*Lz;
                    double beta_sum = 0.0;
                    for (int kk=1; kk<=anie->nz(); kk++){
                        int bnode = anie->nodeIndex(i,j,kk);
                        if (bnode!=anode){
                            double dz = anie->nodePos(bnode).z()-p.z();
                            if (abs(dz)<=Lz){
                                double dzL2 = dz*dz/Lz2;
                                double beta = exp(-dzL2);
                                raez(anode)[bnode] = beta*maescale(bnode);
                                beta_sum += beta;
                            }
                        }
                    }
                    raez(anode)[anode] = -beta_sum*maescale(anode);
                }
            }
        }
    }
}

void TomographicInversion3d::calc_refl_averaging_matrix()
{
    for (int i=1; i<=nb_noded(); i++){
        double x = reflp->x(i);
        double y = reflp->y(i);
        double Lx, Ly;
        if (corr_dep_p != 0){
            Lx = corr_dep_p->atx(x,y);
            Ly = corr_dep_p->aty(x,y);
        }else if (corr_vel_p != 0){
            double Lz;
            corr_vel_p->at(Point3d(x,y,reflp->z(x,y)),Lx,Ly,Lz);
        }else{
            error("TomographicInversion3d::calc_refl_averaging_matrix - no correlation length available");
        }
        if (!xflat()) {
            double Lx2 = std::pow(Lx,2);
            double beta_sum = 0.0;
            for (int ii=1; ii<=nb_noded(); ii++){
                if (ii!=i){
                    double dx = reflp->x(ii)-x;
                    if (abs(dx)<=Lx){
                        double beta = exp(-dx*dx/Lx2);
                        rdx(i)[ii] = beta*mdscale(ii)*refl_weight;
                        beta_sum += beta;
                    }
                }
            }
            rdx(i)[i] = -beta_sum*mdscale(i)*refl_weight;
        }
        if (!yflat()) {
            double Ly2 = std::pow(Ly,2);
            double beta_sum = 0.0;
            for (int jj=1; jj<=nb_noded(); jj++){
                if (jj!=i){
                    double dy = reflp->y(jj)-y;
                    if (abs(dy)<=Ly){
                        double beta = exp(-dy*dy/Ly2);
                        rdy(i)[jj] = beta*mdscale(jj)*refl_weight;
                        beta_sum += beta;
                    }
                }
            }
            rdy(i)[i] = -beta_sum*mdscale(i)*refl_weight;
        }
    }
}
    
void TomographicInversion3d::calc_damping_matrix()
{
    Array2d<double> local_T;
    std::vector<int>    j;
    if (slowness_mesh().flat()) {
        local_T.resize(4,4);
        j.resize(4);
    } else {
        local_T.resize(8,8);
        j.resize(8);
    }
    
    for (int i=1; i<=slowness_mesh().numCells(); i++){
        slowness_mesh().cellNodes(i, j);
        slowness_mesh().cellNormKernel(i, local_T);
        
        if (do_squeezing){
            Point3d p(0,0,0);
            for (int k = 1; k <= j.size(); ++k) {
                p += slowness_mesh().nodePos(j[k-1]);
            }
            p /= j.size();
            double w;
            damping_wt_p->at(p,w);
            local_T *= w;
        }
        add_global(local_T, j, Tv);
    }
}

void TomographicInversion3d::calc_refl_damping_matrix() // variable interval version (beta).
{
    Td_common.resize(4,4);
    Td_common(1,1) = 4; Td_common(1,2) = 2; Td_common(1,3) = 1; Td_common(1,4) = 2; // from the 2D velocity damping matrix in tomo2d
    Td_common(2,1) = 2; Td_common(2,2) = 4; Td_common(2,3) = 2; Td_common(2,4) = 1;
    Td_common(3,1) = 1; Td_common(3,2) = 2; Td_common(3,3) = 4; Td_common(3,4) = 2;
    Td_common(4,1) = 2; Td_common(4,2) = 1; Td_common(4,3) = 1; Td_common(4,4) = 4;
    Td_common *= (1.0/36.0);

    Array1d<double> p(4);
    d_choldc(Td_common.toRecipe(), 4, p.toRecipe());
    for (int i=1; i<=4; i++){
        Td_common(i,i) = p(i);
        for (int j=i+1; j<=4; j++){
            Td_common(i,j) = Td_common(j,i);
        }
    }

    std::vector<int> r(4);
    Array2d<double> Td_local(4,4);
    double fac;
    for (int i=1; i<=reflp->numCells(); i++){
        reflp->cellNodes(i,r[0],r[1],r[2],r[3]);
        reflp->cellNormKernel(i,fac);
        Td_local = Td_common * fac;
        add_global(Td_local, r, Td);
    }
}

void TomographicInversion3d::calc_anid_damping_matrix()
{
    Array2d<double> local_T;
    std::vector<int>    j;
    if (anid->flat()) {
        local_T.resize(4,4);
        j.resize(4);
    } else {
        local_T.resize(8,8);
        j.resize(8);
    }
    
    for (int i=1; i<=anid->numCells(); i++){
        anid->cellNodes(i, j);
        anid->cellNormKernel(i, local_T);
        
        if (do_squeezingD){
            Point3d p(0,0,0);
            for (int k = 1; k <= j.size(); ++k) {
                p += anid->nodePos(j[k-1]);
            }
            p /= j.size();
            double w;
            dampingD_wt_p->at(p,w);
            local_T *= w;
        }
        add_global(local_T, j, Tad);
    }
}

void TomographicInversion3d::calc_anie_damping_matrix()
{
    Array2d<double> local_T;
    std::vector<int>    j;
    if (anie->flat()) {
        local_T.resize(4,4);
        j.resize(4);
    } else {
        local_T.resize(8,8);
        j.resize(8);
    }
    
    for (int i=1; i<=anie->numCells(); i++){
        anie->cellNodes(i, j);
        anie->cellNormKernel(i, local_T);
        
        if (do_squeezingE){
            Point3d p(0,0,0);
            for (int k = 1; k <= j.size(); ++k) {
                p += anie->nodePos(j[k-1]);
            }
            p /= j.size();
            double w;
            dampingE_wt_p->at(p,w);
            local_T *= w;
        }
        add_global(local_T, j, Tae);
    }
}

// from local damping matrix to global damping matrix
void
TomographicInversion3d::add_global(Array2d<double> const& a,
                                   std::vector<int> const& j, sparse_matrix& global) 
{
    assert(j.size() == a.nCol());
    assert(j.size() == a.nRow());
    for (int m=1; m <= j.size(); ++m){
        for ( int n=1; n <= j.size(); ++n){
            global(j[m-1])[j[n-1]] += a(m,n);
        }
    }
}

double TomographicInversion3d::calc_ave_dmv()
{
    double dm_norm=0.0;
    for (int i=1; i<=nnodev; i++){
        dm_norm += dmodel_total(i)*dmodel_total(i);
    }
    return sqrt(dm_norm/nnodev);
}

double TomographicInversion3d::calc_ave_dmd()
{
    if (reflp) {
        double dm_norm=0.0;
        for (int i=1+nnodev; i<=nb_noded()+nnodev; i++){
            dm_norm += dmodel_total(i)*dmodel_total(i);
        }
        return sqrt(dm_norm/nb_noded())*refl_weight;
    } else {
        return 0;
    }
}

double TomographicInversion3d::calc_ave_dmad()
{
    if (ani) {
        double dm_norm=0.0;
        for (int i=1+nnode_subtotal; i<=nnoded+nnode_subtotal; i++){
            dm_norm += dmodel_total(i)*dmodel_total(i);
        }
        return sqrt(dm_norm/nnoded);
    } else {
        return 0;
    }
}

double TomographicInversion3d::calc_ave_dmae()
{
    if (ani) {
        double dm_norm=0.0;
        for (int i=1+nnoded+nnode_subtotal; i<=nnodee+nnoded+nnode_subtotal; i++){
            dm_norm += dmodel_total(i)*dmodel_total(i);
        }
        return sqrt(dm_norm/nnodee);
    } else {
        return 0;
    }
}

double
square_mean(std::vector<double> const& v) {
    return (v.empty()
            ? 0 
            : std::sqrt(std::accumulate(v.begin(), v.end(), 0.0,
                                        [](double d1, double d2) { return d1 + std::pow(d2, 2); })
                        / v.size()));
}

double
square_mean(SparseRectangular const& A, std::vector<double> const& x) {
    std::vector<double> b(A.nRow());
    A.Ax(x, b);
    return square_mean(b);
}

void 
TomographicInversion3d::calc_Lm(opt_double& lmx, opt_double& lmy, double& lmz, opt_double& lmdx, opt_double& lmdy,
				opt_double& lmadx, opt_double& lmady, double& lmadz, opt_double& lmaex, opt_double& lmaey, double& lmaez)
{
    if (smooth_velocity && !fv){
        std::vector<double> dmv(dmodel_total.begin(), dmodel_total.begin()+nnodev);
        lmx = xflat() ? opt_double() : opt_double(square_mean(SparseRectangular(rvx(),nnodev), dmv));
        lmy = yflat() ? opt_double() : opt_double(square_mean(SparseRectangular(rvy(),nnodev), dmv));
        lmz = square_mean(SparseRectangular(rvz(),nnodev), dmv);
    }else{
	lmx = lmy = lmz = -1;
    }
    if (reflp && smooth_depth && !fd){
        std::vector<double> dmd(dmodel_total.begin()+nnodev, dmodel_total.begin()+nnodev+nb_noded());        
        lmdx = xflat() ? opt_double() : opt_double(square_mean(SparseRectangular(rdx(),tmp_nodedr,tmp_nodedr), dmd));
        lmdy = yflat() ? opt_double() : opt_double(square_mean(SparseRectangular(rdy(),tmp_nodedr,tmp_nodedr), dmd));
    }else{
	lmdx = lmdy = -1;
    }
    if (ani && smooth_anid && !fad){
	std::vector<double> dmad(dmodel_total.begin()+nnode_subtotal, dmodel_total.begin()+nnode_subtotal+nnoded);
        lmadx = xflat() ? opt_double() : opt_double(square_mean(SparseRectangular(radx(),nnoded), dmad));
        lmady = yflat() ? opt_double() : opt_double(square_mean(SparseRectangular(rady(),nnoded), dmad));
        lmadz = square_mean(SparseRectangular(radz(),nnoded), dmad);
    }else{
	lmadx = lmady = lmadz = -1;
    }
    if (ani && smooth_anie && !fae){
	std::vector<double> dmae(dmodel_total.begin()+nnode_subtotal+nnoded, dmodel_total.begin()+nnode_subtotal+nnoded+nnodee);
        lmaex = xflat() ? opt_double() : opt_double(square_mean(SparseRectangular(raex(),nnodee), dmae));
        lmaey = yflat() ? opt_double() : opt_double(square_mean(SparseRectangular(raey(),nnodee), dmae));
        lmaez = square_mean(SparseRectangular(raez(),nnodee), dmae);
    }else{
	lmaex = lmaey = lmaez = -1;
    }
}
 
double TomographicInversion3d::calc_chi()
{
    Array1d<double> Adm(ndata_valid);
    SparseRectangular sparseA(A,[&](int i) { return bool(row_index[i-1]);},identity,nnode_total);
    sparseA.Ax(dmodel_total,Adm);
    double val=0.0;
    for (int i=1; i<=nb_data(); i++){
        int j=row_index[i-1];
        if (j>0){
            double res=Adm(j)-data_vec(i);
            val += res*res;
        }
    }
    return val/ndata_valid;
}

int TomographicInversion3d::_solve(bool sv, double wsv, bool sd, double wsd,
				   bool sad, double wsad, bool sae, double wsae,
                                   bool dv, double wdv, bool dd, double wdd,
				   bool dad, double wdad, bool dae, double wdae)
{
    if (is_verbose(0)) {
        std::cerr << "enter TomographicInversion3d::_solve\n";
    }

    // construct total kernel
    Array1d<const sparse_matrix*> As;
    Array1d<SparseMatAux> matspec(0);
    As.push_back(&A);
    matspec.push_back(SparseMatAux(1.0,[&](int i){return bool(row_index[i-1]);},identity));
    int ndata_plus=ndata_valid;
    
    if (sv && !fv){
        if (wsv<0) {
            error("TomographicInversion3d::_solve - negative wsv");
        }
        if (!xflat()) {
            As.push_back(&(rvx()));
            matspec.push_back(SparseMatAux(wsv, unfiltered, identity));
            ndata_plus += nnodev;
        }
        if (!yflat()) {
            As.push_back(&(rvy()));
            matspec.push_back(SparseMatAux(wsv, unfiltered, identity));
            ndata_plus += nnodev;
        }
        {
            As.push_back(&(rvz()));
            matspec.push_back(SparseMatAux(wsv, unfiltered, identity));
            ndata_plus += nnodev;
        }
    }
    if (reflp && sd && !fd){
        if (wsd<0) {
            error("TomographicInversion3d::_solve - negative wsd");
        }
        if (!xflat()) {
            As.push_back(&(rdx()));
            matspec.push_back(SparseMatAux(wsd,tmp_nodedr,tmp_nodedc));
            ndata_plus += nb_noded();
        }
        if (!yflat()) {
            As.push_back(&(rdy()));
            matspec.push_back(SparseMatAux(wsd,tmp_nodedr,tmp_nodedc));
            ndata_plus += nb_noded();
        }
    }
    if (ani && sad && !fad){
        if (wsad<0) {
            error("TomographicInversion3d::_solve - negative wsad");
        }
        if (!xflat()) {
            As.push_back(&(radx()));
            matspec.push_back(SparseMatAux(wsad,tmp_nodeadr,tmp_nodeadc));
            ndata_plus += nnoded;
        }
        if (!yflat()) {
            As.push_back(&(rady()));
            matspec.push_back(SparseMatAux(wsad,tmp_nodeadr,tmp_nodeadc));
            ndata_plus += nnoded;
        }
        {
            As.push_back(&(radz()));
            matspec.push_back(SparseMatAux(wsad,tmp_nodeadr,tmp_nodeadc));
            ndata_plus += nnoded;
        }	
    }
    if (ani && sae && !fae){
        if (wsae<0) {
            error("TomographicInversion3d::_solve - negative wsae");
        }
        if (!xflat()) {
            As.push_back(&(raex()));
            matspec.push_back(SparseMatAux(wsae,tmp_nodeaer,tmp_nodeaec));
            ndata_plus += nnodee;
        }
        if (!yflat()) {
            As.push_back(&(raey()));
            matspec.push_back(SparseMatAux(wsae,tmp_nodeaer,tmp_nodeaec));
            ndata_plus += nnodee;
        }
        {
            As.push_back(&(raez()));
            matspec.push_back(SparseMatAux(wsae,tmp_nodeaer,tmp_nodeaec));
            ndata_plus += nnodee;
        }	
    }
    if (dv && !fv){
        if (wdv<0) {
            error("TomographicInversion3d::_solve - negative wdv");
        }
        As.push_back(&Tv);
        matspec.push_back(SparseMatAux(wdv, unfiltered, identity));
        ndata_plus += nnodev;
    }
    if (reflp && dd && !fd){
        if (wdd<0) {
            error("TomographicInversion3d::_solve - negative wdd");
        }
        As.push_back(&Td);
        matspec.push_back(SparseMatAux(wdd,tmp_nodedr,tmp_nodedc));
        ndata_plus += nb_noded();
    }
    if (ani && dad && !fad){
        if (wdad<0) {
            error("TomographicInversion3d::_solve - negative wdad");
        }
        As.push_back(&Tad);
        matspec.push_back(SparseMatAux(wdad, unfiltered, identity));
        ndata_plus += nnoded;
    }
    if (ani && dae && !fae){
        if (wdae<0) {
            error("TomographicInversion3d::_solve - negative wdae");
        }
        As.push_back(&Tae);
        matspec.push_back(SparseMatAux(wdae, unfiltered, identity));
        ndata_plus += nnodee;
    }

    SparseRectangular B(As,matspec,nnode_total);

    // construct total data vector
    total_data_vec.resize(ndata_plus);
    int idata=1;
    for (int i=1; i<=nb_data(); i++){
        if (row_index[i-1]>0){
            total_data_vec(idata++) = data_vec(i);
        }
    }
    if (sv && !fv){
        if (jumping){
            Array1d<double> Rvm(nnodev), dmv_vec(nnodev);
            for (int i=1; i<=nnodev; i++) {
                dmv_vec(i) = dmodel_total_sum(i);
            }
            if (!xflat()) {
                SparseRectangular sparseRv_x(rvx(),nnodev);
                sparseRv_x.Ax(dmv_vec,Rvm);
                for (int i=1; i<=nnodev; i++) {
                    total_data_vec(idata++) = -wsv*Rvm(i); // Lx
                }
            }
            if (!yflat()) {
                SparseRectangular sparseRv_y(rvy(),nnodev);
                sparseRv_y.Ax(dmv_vec,Rvm);
                for (int i=1; i<=nnodev; i++) {
                    total_data_vec(idata++) = -wsv*Rvm(i); // Ly
                }
            }
            SparseRectangular sparseRv_z(rvz(),nnodev);
            sparseRv_z.Ax(dmv_vec,Rvm);
            for (int i=1; i<=nnodev; i++) {
                total_data_vec(idata++) = -wsv*Rvm(i); // Lz
            }
        }else{
            if (!xflat()) {
                for (int i=1; i<=nnodev; i++) {
                    total_data_vec(idata++) = 0.0; // Lx
                }
            }
            if (!yflat()) {
                for (int i=1; i<=nnodev; i++) {
                    total_data_vec(idata++) = 0.0; // Ly
                }
            }
            for (int i=1; i<=nnodev; i++) {
                total_data_vec(idata++) = 0.0; // Lz
            }
        }
    }
    if (reflp && sd && !fd){
        if (jumping){
            Array1d<double> Rdm(nb_noded()), dmd_vec(nb_noded());
            for (int i=1+nnodev, j=1; i<=nb_noded()+nnodev; i++) {
                dmd_vec(j++) = dmodel_total_sum(i);
            }
            if (!xflat()) {
                SparseRectangular sparseRd_x(rdx(),tmp_nodedr,tmp_nodedr);
                sparseRd_x.Ax(dmd_vec,Rdm);
                for (int i=1; i<=nb_noded(); i++) {
                    total_data_vec(idata++) = -wsd*Rdm(i); // Lx
                }
            }
            if (!yflat()) {
                SparseRectangular sparseRd_y(rdy(),tmp_nodedr,tmp_nodedr);
                sparseRd_y.Ax(dmd_vec,Rdm);
                for (int i=1; i<=nb_noded(); i++) {
                    total_data_vec(idata++) = -wsd*Rdm(i); // Ly
                }
            }
        }else{
            if (!xflat()) {
                for (int i=1; i<=nb_noded(); i++) {
                    total_data_vec(idata++) = 0.0; // Lx
                }
            }
            if (!yflat()) {
                for (int i=1; i<=nb_noded(); i++) {
                    total_data_vec(idata++) = 0.0; // Ly
                }
            }
        }
    }
    if (ani && sad && !fad){
        if (jumping){
            Array1d<double> Radm(nnoded), dmad_vec(nnoded);
            for (int i=1+nnode_subtotal, j=1; i<=nnode_subtotal+nnoded; i++) {
                dmad_vec(j++) = dmodel_total_sum(i);
            }
            if (!xflat()) {
                SparseRectangular sparseRad_x(radx(),nnoded);
                sparseRad_x.Ax(dmad_vec,Radm);
                for (int i=1; i<=nnoded; i++) {
                    total_data_vec(idata++) = -wsad*Radm(i); // Lx
                }
            }
            if (!yflat()) {
                SparseRectangular sparseRad_y(rady(),nnoded);
                sparseRad_y.Ax(dmad_vec,Radm);
                for (int i=1; i<=nnoded; i++) {
                    total_data_vec(idata++) = -wsad*Radm(i); // Ly
                }
            }
            SparseRectangular sparseRad_z(radz(),nnoded);
            sparseRad_z.Ax(dmad_vec,Radm);
            for (int i=1; i<=nnoded; i++) {
                total_data_vec(idata++) = -wsad*Radm(i); // Lz
            }
        }else{
            if (!xflat()) {
                for (int i=1; i<=nnoded; i++) {
                    total_data_vec(idata++) = 0.0; // Lx
                }
            }
            if (!yflat()) {
                for (int i=1; i<=nnoded; i++) {
                    total_data_vec(idata++) = 0.0; // Ly
                }
            }
            for (int i=1; i<=nnoded; i++) {
                total_data_vec(idata++) = 0.0; // Lz
            }
        }
    }
    if (ani && sae && !fae){
        if (jumping){
            Array1d<double> Raem(nnodee), dmae_vec(nnodee);
            for (int i=1+nnode_subtotal+nnoded, j=1; i<=nnode_subtotal+nnoded+nnodee; i++) {
                dmae_vec(j++) = dmodel_total_sum(i);
            }
            if (!xflat()) {
                SparseRectangular sparseRae_x(raex(),nnodee);
                sparseRae_x.Ax(dmae_vec,Raem);
                for (int i=1; i<=nnodee; i++) {
                    total_data_vec(idata++) = -wsae*Raem(i); // Lx
                }
            }
            if (!yflat()) {
                SparseRectangular sparseRae_y(raey(),nnodee);
                sparseRae_y.Ax(dmae_vec,Raem);
                for (int i=1; i<=nnodee; i++) {
                    total_data_vec(idata++) = -wsae*Raem(i); // Ly
                }
            }
            SparseRectangular sparseRae_z(raez(),nnodee);
            sparseRae_z.Ax(dmae_vec,Raem);
            for (int i=1; i<=nnodee; i++) {
                total_data_vec(idata++) = -wsae*Raem(i); // Lz
            }
        }else{
            if (!xflat()) {
                for (int i=1; i<=nnodee; i++) {
                    total_data_vec(idata++) = 0.0; // Lx
                }
            }
            if (!yflat()) {
                for (int i=1; i<=nnodee; i++) {
                    total_data_vec(idata++) = 0.0; // Ly
                }
            }
            for (int i=1; i<=nnodee; i++) {
                total_data_vec(idata++) = 0.0; // Lz
            }
        }
    }
    if (dv && !fv){
        if (jumping){
            Array1d<double> Tvm(nnodev), dmv_vec(nnodev);
            SparseRectangular sparseTv(Tv,nnodev);
            for (int i=1; i<=nnodev; i++) dmv_vec(i) = dmodel_total_sum(i);
            sparseTv.Ax(dmv_vec,Tvm);
            for (int i=1; i<=nnodev; i++) {
                total_data_vec(idata++) = -wdv*Tvm(i); 
            }
        }else{
            for (int i=1; i<=nnodev; i++) {
                total_data_vec(idata++) = 0.0;
            }
        }
    }
    if (reflp && dd && !fd){
        if (jumping){
            Array1d<double> Tdm(nb_noded()), dmd_vec(nb_noded());
            SparseRectangular sparseTd(Td,tmp_nodedr,tmp_nodedr);
            for (int i=1+nnodev, j=1; i<=nb_noded()+nnodev; i++) {
                dmd_vec(j++) = dmodel_total_sum(i);
            }
            sparseTd.Ax(dmd_vec,Tdm);
            for (int i=1; i<=nb_noded(); i++) {
                total_data_vec(idata++) = -wdd*Tdm(i); 
            }
        }else{
            for (int i=1; i<=nb_noded(); i++) {
                total_data_vec(idata++) = 0.0;
            }
        }
    }
    if (ani && dad && !fad){
        if (jumping){
            Array1d<double> Tadm(nnoded), dmad_vec(nnoded);
            SparseRectangular sparseTad(Tad,nnoded);
            for (int i=1+nnode_subtotal, j=1; i<=nnode_subtotal+nnoded; i++){
		dmad_vec(j++) = dmodel_total_sum(i);
	    }
            sparseTad.Ax(dmad_vec,Tadm);
            for (int i=1; i<=nnoded; i++) {
                total_data_vec(idata++) = -wdad*Tadm(i); 
            }
        }else{
            for (int i=1; i<=nnoded; i++) {
                total_data_vec(idata++) = 0.0;
            }
        }
    }
    if (ani && dae && !fae){
        if (jumping){
            Array1d<double> Taem(nnodee), dmae_vec(nnodee);
            SparseRectangular sparseTae(Tae,nnodee);
            for (int i=1+nnode_subtotal+nnoded, j=1; i<=nnode_subtotal+nnoded+nnodee; i++){
		dmae_vec(j++) = dmodel_total_sum(i);
	    }
            sparseTae.Ax(dmae_vec,Taem);
            for (int i=1; i<=nnodee; i++) {
                total_data_vec(idata++) = -wdae*Taem(i); 
            }
        }else{
            for (int i=1; i<=nnodee; i++) {
                total_data_vec(idata++) = 0.0;
            }
        }
    }

    int lsqr_itermax=itermax_LSQR;
    int istop = axb_solver::instance()(B,total_data_vec,dmodel_total,LSQR_ATOL,lsqr_itermax);
    if (is_verbose(0)) {
        std::cerr << "exit TomographicInversion3d::_solve\n";
    }
    return lsqr_itermax;
}

void TomographicInversion3d::fixed_damping(int& iter, int& n, double& wdv, double& wdd, double& wdad, double& wdae)
{
    wdv = weight_d_v;
    wdd = weight_d_d;
    wdad = weight_d_ad;
    wdae = weight_d_ae;
    
    iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		   smooth_anid,weight_s_ad,smooth_anie,weight_s_ae,
                   true,weight_d_v,true,weight_d_d,true,weight_d_ad,true,weight_d_ae);
    n++;

    if (robust){
        if (is_verbose(0)){
            cerr << "TomographicInversion3d:: check outliers... ";
        }
        int ndata_valid_orig=ndata_valid;
        removeOutliers();
        if (is_verbose(0)){
            cerr << ndata_valid_orig-ndata_valid << " found\n";
        }
        if (ndata_valid < ndata_valid_orig){
            if (is_verbose(0)){
                cerr << "TomographicInversion3d:: re-inverting without outliers...\n";
            }
            iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
			   smooth_anid,weight_s_ad,smooth_anie,weight_s_ae,
                           true,weight_d_v,true,weight_d_d,true,weight_d_ad,true,weight_d_ae);
            n++;
        }
    }
}

void TomographicInversion3d::auto_damping(int& iter, int& n, double& wdv, double& wdd, double& wdad, double& wdae)
{
    // check if damping is necessary
    iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		   smooth_anid,weight_s_ad,smooth_anie,weight_s_ae,
                   false,0.0,false,0.0,false,0.0,false,0.0);
    n++;

    if (robust){
        if (is_verbose(0)){
            cerr << "TomographicInversion3d:: check outliers... ";
        }
        int ndata_valid_orig=ndata_valid;
        removeOutliers();
        if (is_verbose(0)){
            cerr << ndata_valid_orig-ndata_valid << " found\n";
        }
        if (ndata_valid < ndata_valid_orig){
            if (is_verbose(0)){
                cerr << "TomographicInversion3d:: re-inverting without outliers...\n";
            }
            iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
			   smooth_anid,weight_s_ad,smooth_anie,weight_s_ae,
                           false,0.0,false,0.0,false,0.0,false,0.0);
            n++;
        }
    }
    double ave_dmv0;
    if(fv){
	ave_dmv0 = -1.0;
    }else{
	ave_dmv0 = calc_ave_dmv();
    }
    double ave_dmd0;
    if(fd){
	ave_dmd0 = -1.0;
    }else{
	ave_dmd0 = calc_ave_dmd();
    }
    double ave_dmad0;
    if(fad){
	ave_dmad0 = -1.0;
    }else{
	ave_dmad0 = calc_ave_dmad();
    }
    double ave_dmae0;
    if(fae){
	ave_dmae0 = -1.0;
    }else{
	ave_dmae0 = calc_ave_dmae();
    }
    if (is_verbose(0)){
        cerr << "\t\tave_dm = " << ave_dmv0*100
             << "%, " << ave_dmd0*100 << "%, "
	     << ave_dmad0*100 << "%, " << ave_dmae0*100 << "%, with no damping.\n";
    }
    wdv = ave_dmv0<=target_dv || !damp_velocity || fv ? 0 : 1;
    wdd = ave_dmd0<=target_dd || !damp_depth || fd ? 0 : 1;
    wdad = ave_dmad0<=target_dad || !damp_anid || fad ? 0 : 1;
    wdae = ave_dmae0<=target_dae || !damp_anie || fae ? 0 : 1;
    if (wdd>0.0) wdd = auto_damping_depth(wdv,wdad,wdae,iter,n);
    if (wdae>0.0) wdae = auto_damping_anie(wdv,wdad,wdd,iter,n);
    if (wdad>0.0) wdad = auto_damping_anid(wdv,wdae,wdd,iter,n);
    if (wdv>0.0) wdv = auto_damping_vel(wdad,wdae,wdd,iter,n);
    
}

double TomographicInversion3d::auto_damping_depth(double wdv, double wdad, double wdae, int& iter, int& n)
{
    // secant and bisection search
    if (is_verbose(0)) cerr << "\t\tsearching weight_d_depth...\n";
    const double wd_max = 1e7; // absolute bound
    double wd1 = 1.0; // initial guess
    double wd2 = 1e2;
    const double ddm = target_dd*0.3; // target accuracy = 30 %
    const int secant_itermax=10;
    double rts, xl, swap, dx, fl, f;

    iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		   smooth_anid,weight_s_ad,smooth_anie,weight_s_ae,
                   damp_velocity,wdv,true,wd1,damp_anid,wdad,damp_anie,wdae); n++;
    fl = calc_ave_dmd()-target_dd;

    if (is_verbose(0)){
        cerr << "\t\tave_dmd = " << (fl+target_dd)*100 << "% at wd1("
             << wd1 << ")\n";
    }
    iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		   smooth_anid,weight_s_ad,smooth_anie,weight_s_ae,
                   damp_velocity,wdv,true,wd2,damp_anid,wdad,damp_anie,wdae); n++;
    f = calc_ave_dmd()-target_dd;

    if (is_verbose(0)){
        cerr << "\t\tave_dmd = " << (f+target_dd)*100 << "% at wd2("
             << wd2 << ")\n";
    }

    if (abs(fl) < abs(f)){
        rts = wd1;
        xl = wd2;
        swap=fl; fl=f; f=swap;
    }else{
        xl = wd1; rts = wd2;
    }
    for (int j=1; j<=secant_itermax; j++){
        dx = (xl-rts)*f/(f-fl);
        if (rts+dx<=0){ // switch to bisection
            xl=rts; fl=f;
            rts *= 0.5;
        }else if(rts+dx>wd_max){
            xl=rts; fl=f;
            rts += (wd_max-rts)*0.5;
        }else{
            xl=rts; fl=f;
            rts+=dx;
        }

        iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		       smooth_anid,weight_s_ad,smooth_anie,weight_s_ae,
                       damp_velocity,wdv,true,rts,damp_anid,wdad,damp_anie,wdae); n++;
        f = calc_ave_dmd()-target_dd;

        if (is_verbose(0)){
            cerr << "\t\tave_dmd = " << (f+target_dd)*100
                 << "% at w_d of " << rts << '\n';
        }
        if (rts > wd_max) break;
        if (abs(f) < ddm) break;
    }
    return rts;
}

double TomographicInversion3d::auto_damping_anie(double wdv, double wdad, double wdd, int& iter, int& n)
{
    // secant and bisection search
    if (is_verbose(0)) cerr << "\t\tsearching weight_d_anie...\n";
    const double wd_max = 1e5; // absolute bound
    double wd1 = 1.0; // initial guess
    double wd2 = 1e2;
    const double ddm = target_dae*0.3; // target accuracy = 30 %
    const int secant_itermax=10;
    double rts, xl, swap, dx, fl, f;

    iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		   smooth_anid,weight_s_ad,smooth_anie,weight_s_ae,
                   damp_velocity,wdv,damp_depth,wdd,damp_anid,wdad,true,wd1); n++;
    fl = calc_ave_dmae()-target_dae;
    if (is_verbose(0)){
        cerr << "\t\tave_dmae = " << (fl+target_dae)*100 << "% at wd1("
             << wd1 << ")\n";
    }
    iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		   smooth_anid,weight_s_ad,smooth_anie,weight_s_ae,
                   damp_velocity,wdv,damp_depth,wdd,damp_anid,wdad,true,wd2); n++;
    f = calc_ave_dmae()-target_dae;
    if (is_verbose(0)){
        cerr << "\t\tave_dmae = " << (f+target_dae)*100 << "% at wd2("
             << wd2 << ")\n";
    }

    if (abs(fl) < abs(f)){
        rts = wd1;
        xl = wd2;
        swap=fl; fl=f; f=swap;
    }else{
        xl = wd1; rts = wd2;
    }
    for (int j=1; j<=secant_itermax; j++){
        dx = (xl-rts)*f/(f-fl);
        if (rts+dx<=0){ // switch to bisection
            xl=rts; fl=f;
            rts *= 0.5;
        }else{
            xl=rts; fl=f;
            rts+=dx;
        }
        iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		       smooth_anid,weight_s_ad,smooth_anie,weight_s_ae,
                       damp_velocity,wdv,damp_depth,wdd,damp_anid,wdad,true,rts); n++;
        f = calc_ave_dmae()-target_dae;
        if (is_verbose(0)){
            cerr << "\t\tave_dmae = " << (f+target_dae)*100
                 << "% at w_d of " << rts << '\n';
        }
        if (rts > wd_max) break;
        if (abs(f) < ddm) break;
    }
    return rts;
}

double TomographicInversion3d::auto_damping_anid(double wdv, double wdae, double wdd, int& iter, int& n)
{
    // secant and bisection search
    if (is_verbose(0)) cerr << "\t\tsearching weight_d_anid...\n";
    const double wd_max = 1e5; // absolute bound
    double wd1 = 1.0; // initial guess
    double wd2 = 1e2;
    const double ddm = target_dad*0.3; // target accuracy = 30 %
    const int secant_itermax=10;
    double rts, xl, swap, dx, fl, f;

    iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		   smooth_anid,weight_s_ad,smooth_anie,weight_s_ae,
                   damp_velocity,wdv,damp_depth,wdd,true,wd1,damp_anie,wdae); n++;
    fl = calc_ave_dmad()-target_dad;
    if (is_verbose(0)){
        cerr << "\t\tave_dmad = " << (fl+target_dad)*100 << "% at wd1("
             << wd1 << ")\n";
    }
    iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		   smooth_anid,weight_s_ad,smooth_anie,weight_s_ae,
                   damp_velocity,wdv,damp_depth,wdd,true,wd2,damp_anie,wdae); n++;
    f = calc_ave_dmad()-target_dad;
    if (is_verbose(0)){
        cerr << "\t\tave_dmad = " << (f+target_dad)*100 << "% at wd2("
             << wd2 << ")\n";
    }

    if (abs(fl) < abs(f)){
        rts = wd1;
        xl = wd2;
        swap=fl; fl=f; f=swap;
    }else{
        xl = wd1; rts = wd2;
    }
    for (int j=1; j<=secant_itermax; j++){
        dx = (xl-rts)*f/(f-fl);
        if (rts+dx<=0){ // switch to bisection
            xl=rts; fl=f;
            rts *= 0.5;
        }else{
            xl=rts; fl=f;
            rts+=dx;
        }
        iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		       smooth_anid,weight_s_ad,smooth_anie,weight_s_ae,
                       damp_velocity,wdv,damp_depth,wdd,true,rts,damp_anie,wdae); n++;
        f = calc_ave_dmad()-target_dad;
        if (is_verbose(0)){
            cerr << "\t\tave_dmad = " << (f+target_dad)*100
                 << "% at w_d of " << rts << '\n';
        }
        if (rts > wd_max) break;
        if (abs(f) < ddm) break;
    }
    return rts;
}

double TomographicInversion3d::auto_damping_vel(double wdad, double wdae, double wdd, int& iter, int& n)
{
    // secant and bisection search
    if (is_verbose(0)) cerr << "\t\tsearching weight_d_vel...\n";
    const double wd_max = 1e5; // absolute bound
    double wd1 = 1.0; // initial guess
    double wd2 = 1e2;
    const double ddm = target_dv*0.3; // target accuracy = 30 %
    const int secant_itermax=10;
    double rts, xl, swap, dx, fl, f;

    iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		   smooth_anid,weight_s_ad,smooth_anie,weight_s_ae,
                   true,wd1,damp_depth,wdd,damp_anid,wdad,damp_anie,wdae); n++;
    fl = calc_ave_dmv()-target_dv;
    if (is_verbose(0)){
        cerr << "\t\tave_dmv = " << (fl+target_dv)*100 << "% at wd1("
             << wd1 << ")\n";
    }
    iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		   smooth_anid,weight_s_ad,smooth_anie,weight_s_ae,
                   true,wd2,damp_depth,wdd,damp_anid,wdad,damp_anie,wdae); n++;
    f = calc_ave_dmv()-target_dv;
    if (is_verbose(0)){
        cerr << "\t\tave_dmv = " << (f+target_dv)*100 << "% at wd2("
             << wd2 << ")\n";
    }

    if (abs(fl) < abs(f)){
        rts = wd1;
        xl = wd2;
        swap=fl; fl=f; f=swap;
    }else{
        xl = wd1; rts = wd2;
    }
    for (int j=1; j<=secant_itermax; j++){
        dx = (xl-rts)*f/(f-fl);
        if (rts+dx<=0){ // switch to bisection
            xl=rts; fl=f;
            rts *= 0.5;
        }else{
            xl=rts; fl=f;
            rts+=dx;
        }
        iter += _solve(smooth_velocity,weight_s_v,smooth_depth,weight_s_d,
		       smooth_anid,weight_s_ad,smooth_anie,weight_s_ae,
                       true,rts,damp_depth,wdd,damp_anid,wdad,damp_anie,wdae); n++;
        f = calc_ave_dmv()-target_dv;
        if (is_verbose(0)){
            cerr << "\t\tave_dmv = " << (f+target_dv)*100
                 << "% at w_d of " << rts << '\n';
        }
        if (rts > wd_max) break;
        if (abs(f) < ddm) break;
    }
    return rts;
}

Array1d<double> const&
TomographicInversion3d::velocity_model() const {
    assert(!my_velocity_model.empty()); 
    return my_velocity_model;
}

void
TomographicInversion3d::printSynTime(ostream& os, double vred) const //modified
{
    const double inv_dt = 0.01; // assume 10 ms error
    
    os << src.size() << '\n'; 
    for (int isrc=1; isrc<=src.size(); isrc++){
	os << 's' << " "
	   << src(isrc).x() << " "
	   << src(isrc).y() << " "
	   << src(isrc).z() << " "
	   << rcv(isrc).size() << '\n';
	for (int ircv=1; ircv<=rcv(isrc).size(); ircv++){
	    double redtime = syn_ttime(isrc)(ircv);
	    if (vred!=0){
		double dist_x = rcv(isrc)(ircv).x()-src(isrc).x();
		double dist_y = rcv(isrc)(ircv).y()-src(isrc).y();
		double dist = sqrt(dist_x*dist_x+dist_y*dist_y);
		redtime -= dist/vred;
	    }
	    os << 'r' << " "
	       << rcv(isrc)(ircv).x() << " "
	       << rcv(isrc)(ircv).y() << " "
	       << rcv(isrc)(ircv).z() << " "
	       << get_ray_type(isrc,ircv) << " "
	       << redtime << " "
	       << obs_dt(isrc)(ircv) << '\n';
	}
    }
}
