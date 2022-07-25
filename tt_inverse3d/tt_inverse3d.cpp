/*
 * tt_inverse3d.cc - traveltime inversion 
 * based on tt_inverse.cc
 *
 * Adria Melendez, Jun Korenaga and Alain Miniussi
 * Fall 2011
 * adria.melendez.catalan@gmail.com
 */

#include <iostream>
#include <fstream>
#include <cstdio>

#include <sysexits.h>
#include <unistd.h>
#include <omp.h>

#include "boost/mpi.hpp"
#include "boost/program_options.hpp"

#include "axb_solver.hpp"
#include "in_house_solver.hpp"
#include "in_house_omp_solver.hpp"
#include "array.hpp"
#include "inverse3d.hpp"

//includes for capturing cout
#include <cassert>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
//end modification

void
wait_debugger() {
  int i = 0;
  char hostname[256];
  gethostname(hostname, sizeof(hostname));
  printf("PID %d on %s ready for attach\n", getpid(), hostname);
  fflush(stdout);
  while (0 == i) {
	sleep(5);
  }
}

int main(int argc, char** argv)
{

 /* Capture cout to file my_file.log
int fdum = open("tt_inverse3d_cout.log", O_WRONLY|O_CREAT|O_TRUNC, 0660);
   assert(fdum >= 0);
   int ret = dup2(fdum, 1);
   assert(ret >= 0);
   printf("cout for tt_inverse3d\n");
   std::cout << "Starting log!" << std::endl;
   close(fdum);
// end modification */

    namespace po = boost::program_options;

    bool getMesh=false, getData=false;
    bool forward=false;
    bool first_iter=false;
    double vred=0.0;
    char *meshfn, *datafn; //, *synfn;
//    char *anidfn, *aniefn;
    int xorder=3, yorder=3, zorder=3, nintp=8;
    double crit_len=-1, bend_cg_tol=1e-4, bend_br_tol=1e-7;
    double crit_chi=-1, lsqr_atol=1e-3;
    double refl_weight=1;
    double tilt=0;
    bool getRefl=false, outMask=false;
    bool getAniD=false, getAniE=false;
    bool getVh=false;
    bool outMaskRefl=false;
    bool outMaskAnid=false;
    bool outMaskAnie=false;
    bool doFullRefl=false;
    std::string synfn;
    std::string reflfn;
    std::string anidfn, aniefn;
    std::string logfile;
    std::string output_prefix = "inv_out";
    int outlevel=0;
    bool gotError=false;
    int verbosity_level = -1;

    int niter=1;
    bool smooth_vel=false, smooth_dep=false;
    bool smooth_anid=false, smooth_anie=false;
    bool apply3Dfilter=false;
    char *boundfn;
    double wsv_min=-1, wsv_max=-1, dwsv=-1;
    double wsd_min=-1, wsd_max=-1, dwsd=-1;
    double wsad_min=-1, wsad_max=-1, dwsad=-1;
    double wsae_min=-1, wsae_max=-1, dwsae=-1;
    double max_dv=-1, max_dd=-1, max_dad=-1, max_dae=-1;
    double wdv=0, wdd=0, wdad=0, wdae=0;
    double target_chisq=0;
    bool vlogscale=false, dlogscale=false;
    bool adlogscale=false, aelogscale=false;
    bool getCorrVel=false, getCorrDep=false;
    bool getCorrAniD=false, getCorrAniE=false;
    char *corr_velfn, *corr_depfn;
    char *corr_anidfn, *corr_aniefn;
    bool printFinalOnly=false;
    bool auto_damping=false, fixed_damping=false, jumping=false;
    bool getDamp3d=false, getDamp3dD=false, getDamp3dE=false;
    char *damp_velfn, *damp_anidfn, *damp_aniefn;

    bool fv=false, fd=false, fad=false, fae=false;

    int nb_threads = 1;
    bool dump_placement = false;
    bool debug_wait_asked = false;

    std::string solver_name      = "default";
    std::string solver_trace_dir = ".";
    bool        solver_stats     = false;
    int         solver_trace_limit = -1;
    
    po::options_description opts("Among supported options:");
    opts.add_options() 
        ("help,h", "Print this help message.")
        ("placement", "Print placement of the processes at startup.")
        ("wait-db", "wait for the debugger.")
        ("nb-threads", po::value<int>(&nb_threads), "The number of thread to use.")
        ("solver", po::value<std::string>(&solver_name)->default_value("omp"), "select the solver name. Supported values are: \"omp\", \"default\".")
        ("solver-stats", "print solver time along with problem sizes.")
        ("solver-trace-dir", po::value<std::string>(&solver_trace_dir), "specify the directoy name where to dump the solver data dump.")
        ("solver-trace-time", po::value<int>(&solver_trace_limit), "only dum solver data when solving time is above that number of seconds.");
    {
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv)
                  .options(opts).allow_unregistered().run(),
                  vm);
        po::notify(vm);
        dump_placement   = vm.count("placement")   > 0;
        debug_wait_asked = vm.count("wait-db")     > 0;
        solver_stats     = vm.count("solver-stats") > 0;
        if (vm.count("help") > 0) {
            std::cerr << opts;
            return 1;
        }
    }
    for (int i=1; i<argc; i++){
        if (argv[i][0] == '-'){
            switch(argv[i][1]){
	    case 'M':
		meshfn = &argv[i][2];
		getMesh = true;
		break;
	    case 'G':
		datafn = &argv[i][2];
		getData = true;
		break;
	    case 'p':
		int a, b, c, d;
		if(sscanf(&argv[i][2], "%d/%d/%d/%d", &a, &b, &c, &d)==4){
		    if((a==0 || a==1) && (b==0 || b==1) && (c==0 || c==1) && (d==0 || d==1)){
			fv=a;
			fd=b;
			fad=c;
			fae=d;
		    }else{
			error("the four values for -p option must be either 0 or 1.\n");
		    }
		}else{
		    error("-p option requires four 0 or 1 values.\n");
		}
		break;
	    case 'f':
   		std::cout << "##################################################" << std::endl;
   		std::cout << "Generating Synthetic Traveltimes!" << std::endl;
   		std::cout << "##################################################" << std::endl << std::endl;
		forward = true;
		synfn = &argv[i][2];
		break;
	    case 'i':
		first_iter = true;
		break;
	    case 'r':
		vred = atof(&argv[i][2]);
		break;
	    case 'N':
	    {
		int tmp1a, tmp1b, tmp1c, tmp3;
		double tmp2, tmp4, tmp5;
		if (sscanf(&argv[i][2], "%d/%d/%d/%lf/%d/%lf/%lf",
			   &tmp1a, &tmp1b, &tmp1c, &tmp2, &tmp3, &tmp4, &tmp5)==7){
		    if (tmp1a>0 && tmp1b>0 && tmp1c>0 && tmp2>0 && tmp3 > 0 && tmp4>0 && tmp5>0){
			xorder=tmp1a; 
                        yorder=tmp1b; 
                        zorder=tmp1c;
			crit_len=tmp2;
			nintp=tmp3;
			bend_cg_tol=tmp4;
			bend_br_tol=tmp5;
		    }else{
			error("invalid -N option, ignored.\n");
		    }
		}
		break;
	    }
	    case 'P':
		jumping = true; // pure jumping strategy
		break;
	    case 'R':
		crit_chi = atof(&argv[i][2]); // robust inversion
		break;
	    case 'Q':
		lsqr_atol = atof(&argv[i][2]);
		break;
	    case 'A':
		doFullRefl = true;
		break;
	    case 'F':
		reflfn = &argv[i][2];
		break;
	    case 'd':
		anidfn = &argv[i][2];
		getAniD = true;
		break;
	    case 'e':
		aniefn = &argv[i][2];
		getAniE = true;
		break;
	    case 'H':
		getVh = true;
		break;
	    case 't':
		tilt = atof(&argv[i][2]);
		break;
	    case 'W':
		refl_weight = atof(&argv[i][2]);
		break;
	    case 'L':
                logfile = &argv[i][2];
		break;
	    case 'O':
		output_prefix = &argv[i][2];
		break;
	    case 'o':
		outlevel = atoi(&argv[i][2]);
		break;
	    case 'l':
		printFinalOnly = true;
		break;
	    case 'K':
		outMask = true;
		break;
	    case 'k':
		outMaskRefl = true;
		break;
	    case 'x':
 	        outMaskAnid = true;
		break;
	    case 'y':
	        outMaskAnie = true;
		break;
	    case 'I':
		niter = atoi(&argv[i][2]);
		break;
	    case 'J':
		target_chisq = atof(&argv[i][2]);
		break;
	    case 'S':
		switch (argv[i][2]){
		case 'V':
		{
		    double a, b, c;
		    int nitem=sscanf(&argv[i][3], "%lf/%lf/%lf", &a, &b, &c);
		    if (nitem==1){ wsv_min=a; wsv_max=a; dwsv=a+1.0; }
		    else if (nitem==3){ wsv_min=a; wsv_max=b; dwsv=c; }
		    else{
			cerr << "invalid -SV option.\n";
			gotError=true;
		    }
		    smooth_vel = true;
		    break;
		}
		case 'D':
		{
		    double a, b, c;
		    int nitem=sscanf(&argv[i][3], "%lf/%lf/%lf", &a, &b, &c);
		    if (nitem==1){ wsd_min=a; wsd_max=a; dwsd=a+1.0; }
		    else if (nitem==3){ wsd_min=a; wsd_max=b; dwsd=c; }
		    else{
			cerr << "invalid -SD option.\n";
			gotError=true;
		    }
		    smooth_dep = true;
		    break;
		}
		case 'd':
		{
		    double a, b, c;
		    int nitem=sscanf(&argv[i][3], "%lf/%lf/%lf", &a, &b, &c);
		    if (nitem==1){ wsad_min=a; wsad_max=a; dwsad=a+1.0; }
		    else if (nitem==3){ wsad_min=a; wsad_max=b; dwsad=c; }
		    else{
			cerr << "invalid -Sd option.\n";
			gotError=true;
		    }
		    smooth_anid = true;
		    break;
		}
		case 'e':
		{
		    double a, b, c;
		    int nitem=sscanf(&argv[i][3], "%lf/%lf/%lf", &a, &b, &c);
		    if (nitem==1){ wsae_min=a; wsae_max=a; dwsae=a+1.0; }
		    else if (nitem==3){ wsae_min=a; wsae_max=b; dwsae=c; }
		    else{
			cerr << "invalid -Se option.\n";
			gotError=true;
		    }
		    smooth_anie = true;
		    break;
		}
		default:
		    cerr << "invalid -S option.\n";
		    gotError=true;
		    break;
		}
		break;
	    case 's':
		apply3Dfilter = true;
		boundfn = &argv[i][2];
		break;
	    case 'X':
		switch (argv[i][2]){
		case 'V':
		    vlogscale=true;
		    break;
		case 'D':
		    dlogscale=true;
		    break;
		case 'd':
		    adlogscale=true;
		    break;
		case 'e':
		    aelogscale=true;
		    break;	    
		default:
		    cerr << "invalid -X option.\n";
		    gotError = true;
		    break;
		}
		break;
	    case 'T': // auto damping
		auto_damping = true;
		switch (argv[i][2]){
		case 'V':
		    max_dv = atof(&argv[i][3]);
		    break;
		case 'D':
		    max_dd = atof(&argv[i][3]);
		    break;
		case 'd':
		    max_dad = atof(&argv[i][3]);
		    break;
		case 'e':
		    max_dae = atof(&argv[i][3]);
		    break;
		default:
		    cerr << "invalid -T option.\n";
		    gotError = true;
		    break;
		}
		break; 
	    case 'D': // fixed damping
		fixed_damping = true;
		switch (argv[i][2]){
		case 'V':
		    wdv = atof(&argv[i][3]);
		    break;
		case 'D':
		    wdd = atof(&argv[i][3]);
		    break;
		case 'd':
		    wdad = atof(&argv[i][3]);
		    break;
		case 'e':
		    wdae = atof(&argv[i][3]);
		    break;
		case 'Q':
		    getDamp3d=true;
		    damp_velfn = &argv[i][3];
		    break;
		case 'R':
		    getDamp3dD=true;
		    damp_anidfn = &argv[i][3];
		    break;
		case 'S':
		    getDamp3dE=true;
		    damp_aniefn = &argv[i][3];
		    break;
		default:
		    cerr << "invalid -D option.\n";
		    gotError = true;
		    break;
		}
		break; 
	    case 'C':
		switch(argv[i][2]){
		case 'V':
		    getCorrVel=true;
		    corr_velfn = &argv[i][3];
		    break;
		case 'D':
		    getCorrDep=true;
		    corr_depfn = &argv[i][3];
		    break;
		case 'd':
		    getCorrAniD=true;
		    corr_anidfn = &argv[i][3];
		    break;
		case 'e':
		    getCorrAniE=true;
		    corr_aniefn = &argv[i][3];
		    break;
		default:
		    cerr << "invalid -C option.\n";
		    gotError = true;
		    break;
		}
		break;
	    case 'V':
		verbosity_level = 0;
		if (isdigit(argv[i][2])){
		    verbosity_level = atoi(&argv[i][2]);
		}
		break;
	    default:
		// Maybe dealt with the option parser
		break;
	    }
	} else {
            // Maybe dealt with the option parser
            break;
	}
    }

    if (solver_name == "omp")  {
        select_in_house_omp_solver();
    } else if (solver_name == "nothread" || solver_name == "default" )  {
        select_in_house_solver();
    } else {
        std::cerr << "unknown solver " << solver_name << ". Bye\n";
        return -1;
    }

    cout << "TOMO 3D v0.3 - modified by A. Loureiro JUL2022\n";
    if (!getMesh || !getData){
	cerr << "A velocity model and a data file are required.\n";
	gotError=true;
    }
    if ((getAniD && !getAniE) || (!getAniD && getAniE)){
	cerr << "Anisotropy inversion requires two input files.\n";
	gotError=true;
    }
    if (smooth_vel && !getCorrVel){
	cerr << "Velocity smoothing needs correlation length info.\n";
	gotError=true;
    }
    if (smooth_dep && (!getCorrDep && !getCorrVel)){
	cerr << "Depth smoothing needs correlation length info.\n";
	gotError=true;
    }
    if (smooth_anid && !getCorrAniD){
	if (getCorrAniE){
	    corr_anidfn = corr_aniefn;
	    cerr << "Using epsilon correlation lengths for delta smoothing.\n";
	    getCorrAniD = true;
	}else if (!getCorrAniE && getCorrVel){
	    corr_anidfn = corr_velfn;
	    cerr << "Using velocity correlation lengths for delta smoothing.\n";
	    getCorrAniD = true;
	}else{
	    cerr << "Anisotropy smoothing requires at least one file for correlation length info.\n";
	    gotError=true;
	}
    }
    if (smooth_anie && !getCorrAniE){
	if (getCorrAniD){
	    corr_aniefn = corr_anidfn;
	    cerr << "Using delta correlation lengths for epsilon smoothing.\n";
	    getCorrAniE = true;
	}else if (!getCorrAniD && getCorrVel){
	    corr_aniefn = corr_velfn;
	    cerr << "Using velocity correlation lengths for epsilon smoothing.\n";
	    getCorrAniE = true;
	}else{
	    cerr << "Anisotropy smoothing requires at least one file for correlation length info.\n";
	    gotError=true;
	}
    }    
    if (auto_damping && fixed_damping){
	cerr << "-T and -D options are mutually exclusive.\n";
	gotError=true;
    }
    if (!forward && first_iter){
	cerr << "-i must be used with -f.\n";
	gotError=true;
    }
    if (gotError) error("usage: tt_inverse3d ...");
    
    if (debug_wait_asked) {
        wait_debugger();
    }
    namespace mpi = boost::mpi;
    mpi::environment env(argc, argv, mpi::threading::funneled);
    mpi::communicator world;
    
    if (world.rank() == 0) {
        std::cout << "MPI threading level is " << env.thread_level() << '\n';
        axb_solver::instance().computation_statistics(solver_stats);
        axb_solver::instance().trace_computation_above(solver_trace_limit);
        axb_solver::instance().computation_trace_dir(solver_trace_dir);
    }
    if (dump_placement) {
        for(int p = 0; p < world.size(); ++p) {
            world.barrier();
            if (p == world.rank()) {
                std::cout << "Process " << p << " on " << env.processor_name() << " with " << nb_threads << "threads" << std::endl;
            }
        }
    }

    
    SlownessMesh3d smesh(meshfn);
    TomographicInversion3d inv(world,smesh,datafn,xorder,yorder,zorder,crit_len,
			       nintp,bend_cg_tol,bend_br_tol);
    inv.set_output_prefix(output_prefix);
    inv.set_nb_threads(nb_threads);
#pragma omp parallel
    omp_set_num_threads(nb_threads);
    if (crit_chi>0) inv.doRobust(crit_chi);
    inv.setLSQR_TOL(lsqr_atol);
    if (printFinalOnly){
	inv.outFinal(outlevel);
    }else{
	inv.outStepwise(outlevel);
    }

    if (!logfile.empty()) {
        inv.set_log(logfile);
    }
    if (outMask) inv.outMask();
    if (inv.need_reflections() && outMaskRefl) inv.outMaskRefl();
    if (getAniD && outMaskAnid) inv.outMaskAnid();
    if (getAniE && outMaskAnie) inv.outMaskAnie();
    inv.set_verbosity_level(verbosity_level);

    if (reflfn.empty()) {
        if (inv.need_reflections()) {
	  std::cerr << "Missing reflector input file.\n";
	  world.abort(EX_DATAERR);
        }
    } else {
        boost::shared_ptr<Interface3d> reflp(new Interface3d(reflfn));
	inv.add_reflections(reflp);
	inv.setReflWeight(refl_weight);
	if (doFullRefl) {
            inv.doFullRefl();
        }
    }

    if(getAniD && getAniE){
	if (tilt < 0 || tilt >= acos(-1.0)){
	    cerr << "Tilt angle must be 0 or positive and smaller than Pi. Check tilt definition in the user guide.\n";
	    world.abort(EX_DATAERR);
	}else{
	    inv.tiltAngle(tilt);
	    smesh.tiltAngle(tilt);
	}
    }

    if (anidfn.empty() && getAniD && aniefn.empty() && getAniE){ // Here pass anidfn & aniefn to inv.
	std::cerr << "Missing anisotropy input file.\n";
	world.abort(EX_DATAERR);
    }else if (!anidfn.empty() && getAniD && !aniefn.empty() && getAniE){
	boost::shared_ptr<AnisotropyMesh3d> anidp(new AnisotropyMesh3d(anidfn));
	boost::shared_ptr<AnisotropyMesh3d> aniep(new AnisotropyMesh3d(aniefn));
	inv.add_anisotropy(anidp,aniep);
	smesh.add_anisotropy(anidp,aniep);
	if (getVh){
	    inv.usingVh();
	    smesh.usingVh();
	}
    }

    if (smooth_vel){
	inv.SmoothVelocity(corr_velfn,wsv_min,wsv_max,dwsv,vlogscale);
	if (apply3Dfilter) inv.applyFilter(boundfn);
    }
    if (smooth_dep){
	if (getCorrDep){
	    inv.SmoothDepth(corr_depfn,wsd_min,wsd_max,dwsd,dlogscale);
	}else{
	    // use velocity's correlation length
	    inv.SmoothDepth(wsd_min,wsd_max,dwsd,dlogscale);
	}
    }
    if (smooth_anid && getCorrAniD){ // Here pass corr_anidfn to inv.
	inv.SmoothAniD(corr_anidfn,wsad_min,wsad_max,dwsad,adlogscale);
    }
    if (smooth_anie && getCorrAniE){ // Here pass corr_aniefn to inv.
	inv.SmoothAniE(corr_aniefn,wsae_min,wsae_max,dwsae,aelogscale);
    }

    if (target_chisq>0) inv.targetChisq(target_chisq);

    // In this if-condition include anisotropy damping
    if (auto_damping){
	if (max_dv>0) inv.DampVelocity(max_dv);
	if (max_dd>0) inv.DampDepth(max_dd);
	if (max_dad>0) inv.DampAniD(max_dad);
	if (max_dae>0) inv.DampAniE(max_dae);
    }else if (fixed_damping){
	if (wdv>=0 && wdd>=0 && wdad>=0 && wdae>=0){
	    inv.FixDamping(wdv,wdd,wdad,wdae);
	    if (getDamp3d){
		inv.Squeezing(damp_velfn);
	    }
	    if (getDamp3dD){
		inv.SqueezingD(damp_anidfn);
	    }
	    if (getDamp3dE){
		inv.SqueezingE(damp_aniefn);
	    }
	}
    }
    if (jumping) inv.doJumping();

    if (forward){
	niter = 1;
	inv.onlyForward(first_iter);
    }

    inv.fixing(fv,fd,fad,fae);
    inv.solve(niter);
    
    if (forward){
	ofstream os(synfn);
	inv.printSynTime(os,vred);
    }

    return 0;
}

/*
 x [0] = -0.220238
 x [1] = -0.029762
 x [2] =  0.395089
 x [3] = -0.290179
 x [4] =  0.116071
*/
