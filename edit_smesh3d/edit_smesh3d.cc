/*
 * edit_smesh3d.cc based on edit_smesh.cc
 *
 * usage: edit_smesh3d smesh_file -C<cmd> [-Lvcorr_file -Uupper_file]
 *
 *        -C<cmd>
 *          cmd = 'a' - set 1-D average
 *                'p<smesh>' - paste <smesh> on the original smesh
 *                'P<1dprof>' - paste 1d-prof 
 *                'Z<s1dprof>' - paste 1d-prof hanging from upperbound file
 *                'sx/y/z' - apply Gaussian smoothing operator with a box of (x,y,z) (km)
 *                'rmx/my/mz' - refine mesh by mx for x-direction, by my for y-direction
 *                           and by mz for z-direction
 *                'cA/x/y/z' - add checkerboard pattern
 *                           (A:amplitude [%], with x, y, and z cycles in km)
 *                'dA/xmin/xmax/ymin/ymax/zmin/zmax' - add rectangular anomaly
 *                           (A:amplitude [%])
 *                'gA/x0/y0/z0/Lx/Ly/Lz' - add gaussian anomaly
 *                'l' - remove low velocity zone
 *                'R<seed>/A/nrand' - randomize the velocity field
 *                'S<seed>/A/xmin/xmax/dx/ymin/ymax/dy/zmin/zmax/dz' 
 *                'G<seed>/A/N/xmin/xmax/ymin/ymax/zmin/zmax'
 *                'm<v>/mohofile' - set sub-Moho velocity as <v>
 *                'v' - keep velocity vertically constant (use with -U)
 *
 *        -Lvcorr.file
 *        -Uupperbound.file
 *
 * Adria Melendez and Jun Korenaga
 * Fall 2011
 */

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "array.h"
#include "util3d.h"
#include "smesh3d.h"
#include "corrlen3d.h"
#include "interface3d.h"

int main(int argc, char **argv)
{
    bool getMesh=false, getCmd=false, setAverage=false;
    bool pasteSmesh=false, smooth=false, checkerboard=false;
    bool paste1dprof=false, sheared1dprof=false, constVert=false;
    bool refineMesh=false, rectangular=false, removeLVZ=false;
    bool gaussian=false, randomize=false, randomize2=false, randomize3=false;
    bool verbose=false, err=false;
    bool getCorr=false, getBound=false;
    char *sfn, *pfn, *corrfn, *prof1dfn, *sprof1dfn;
    double Lx, Ly, Lz, A, cx, cy, cz;
    double dxmin, dxmax, dymin, dymax, dzmin, dzmax;
    double x0, y0, z0;
    int mx, my, mz, seed, nrand, N;
    double rand_xmax, rand_xmin, rand_dx, rand_ymin, rand_ymax, rand_dy, rand_zmax, rand_zmin, rand_dz;
    bool subMohoVel=false;
    double vsubmoho;
    char uboundfn[MaxStr];
    
    for (int i=1; i<argc; i++){
	if (argv[i][0] == '-'){
	    switch(argv[i][1]){
	    case 'C':
		getCmd = true;
		switch(argv[i][2]){
		case 'a':
		    setAverage = true;
		    break;
		case 'p':
		    pasteSmesh = true;
		    pfn = &argv[i][3];
		    break;
		case 'P':
		    paste1dprof = true;
		    prof1dfn = &argv[i][3];
		    break;
		case 'Z':
		    sheared1dprof = true;
		    sprof1dfn = &argv[i][3];
		    break;
		case 's':
		    if (sscanf(&argv[i][3], "%lf/%lf/%lf", &Lx, &Ly, &Lz)==3){
			smooth = true;
		    }else{
			cerr << "invalid -Cs option\n";
			err = true;
		    }
		    break;
		case 'c':
		    if (sscanf(&argv[i][3], "%lf/%lf/%lf/%lf", &A, &cx, &cy, &cz)==4){
			checkerboard = true;
		    }else{
			cerr << "invalid -Cc option\n";
			err = true;
		    }
		    break;
		case 'd':
		    if (sscanf(&argv[i][3], "%lf/%lf/%lf/%lf/%lf/%lf/%lf",
			       &A, &dxmin, &dxmax, &dymin, &dymax, &dzmin, &dzmax)==7){
			rectangular=true;
		    }else{
			cerr << "invalid -Cd option\n";
			err = true;
		    }
		    break;
		case 'g':
		    if (sscanf(&argv[i][3], "%lf/%lf/%lf/%lf/%lf/%lf/%lf",
			       &A, &x0, &y0, &z0, &Lx, &Ly, &Lz)==7){
			gaussian=true;
		    }else{
			cerr << "invalid -Cg option\n";
			err = true;
		    }
		    break;
		case 'r':
		    if (sscanf(&argv[i][3], "%d/%d/%d", &mx, &my, &mz)==3){
			refineMesh = true;
		    }else{
			cerr << "invalid -Cr option\n";
			err = true;
		    }
		    break;
		case 'l':
		    removeLVZ = true;
		    break;
		case 'R':
		    if (sscanf(&argv[i][3], "%d/%lf/%d",
			       &seed, &A, &nrand)==3){
			randomize = true;
		    }else{
			cerr << "invalid -CR option\n";
			err = true;
		    }
		    break;
		case 'S':
		    if (sscanf(&argv[i][3], "%d/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf/%lf",
			       &seed, &A, &rand_xmin, &rand_xmax, &rand_dx,
			       &rand_ymin, &rand_ymax, &rand_dy,
			       &rand_zmin, &rand_zmax, &rand_dz)==11){
			randomize2=true;
		    }else{
			cerr << "invalid -CS option\n";
			err = true;
		    }
		    break;
		case 'G':
		    if (sscanf(&argv[i][3], "%d/%lf/%d/%lf/%lf/%lf/%lf/%lf/%lf",
			       &seed, &A, &N, &rand_xmin, &rand_xmax,
			       &rand_ymin, &rand_ymax,
			       &rand_zmin, &rand_zmax)==9){
			randomize3=true;
		    }else{
			cerr << "invalid -CG option\n";
			err = true;
		    }
		    break;
		case 'm':
		    if (sscanf(&argv[i][3], "%lf/%[^/]", &vsubmoho, uboundfn)==2){
			subMohoVel = true;
		    }else{
			cerr << "invalid -Cm option\n";
			err = true;
		    }
		    break;
		case 'v':
		    constVert = true;
		    break;
		default:
		    cerr << "invalid -C option\n";
		    err = true;
		    break;
		}
		break;
	    case 'L':
		getCorr = true;
		corrfn = &argv[i][2];
		break;
	    case 'U':
		getBound = true;
		sscanf(&argv[i][2], "%s", uboundfn);
		break;
	    default:
		err = true;
		break;
	    }
	}else{
	    getMesh = true;
	    sfn = argv[i];
	}
    }

    if (!getMesh || !getCmd) err=true;
    if (err) error("usage: edit_smesh3d [ -options ]");

    // read smesh
    ifstream s(sfn);
    if (!s){
	cerr << "edit_smesh3d::cannot open " << sfn << "\n";
	exit(1);
    }

    int nx, ny, nz;
    double v_water, v_air;
    s >> nx >> ny >> nz >> v_water >> v_air;
    int nnodes = nx*ny*nz;
    
    Array1d<double> xpos(nx), ypos(ny), zpos(nz);
    Array2d<double> topo(nx,ny);
    Array3d<double> vgrid(nx,ny,nz);
    for (int i=1; i<=nx; i++) s >> xpos(i);
    for (int i=1; i<=ny; i++) s >> ypos(i);
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    s >> topo(i,j);
	}
    }
    for (int k=1; k<=nz; k++) s >> zpos(k);
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    for (int k=1; k<=nz; k++){
		double vel;
		s >> vel;
		vgrid(i,j,k) = vel;
	    }
	}
    }

    // 3-D correlation length
    CorrelationLength3d *corr_p;
    if (getCorr){
	corr_p = new CorrelationLength3d(corrfn);
    }

    // upper bound
    Interface3d *ubound_p;
    Array2d<int> kstart(nx,ny);
    kstart = 1;
    if (getBound || subMohoVel){
	ubound_p = new Interface3d(uboundfn);
	for (int i=1; i<=nx; i++){
	    for (int j=1; j<=ny; j++){
		double boundz = ubound_p->z(xpos(i),ypos(j));
		for (int k=1; k<=nz; k++){
		    if (zpos(k)+topo(i,j)>boundz) break;
		    kstart(i,j) = k;
		}
		if (kstart(i,j)>1) kstart(i,j)++;
	    }
	}
    }
    
    // do something
    if (setAverage){
	for (int k=1; k<=nz; k++){
	    double ave=0;
	    for (int i=1; i<=nx; i++){
		for (int j=1; j<=ny; j++){
		    ave += vgrid(i,j,k);
		}
	    }
	    ave /= nx*ny;
	    for (int i=1; i<=nx; i++){
		for (int j=1; j<=ny; j++){
		    if(k>=kstart(i,j)){
			vgrid(i,j,k) = ave;
		    }
		}
	    }
	}
    }else if (pasteSmesh){
	SlownessMesh3d smesh1(pfn);
	for (int i=1; i<=nx; i++){
	    double x = xpos(i);
	    for (int j=1; j<=ny; j++){
		double y = ypos(j);
		double t = topo(i,j);
		for (int k=kstart(i,j); k<=nz; k++){
		    double z = t+zpos(k);
		    Point3d p(x,y,z);
		    if (!smesh1.inWater(p) && !smesh1.inAir(p)){
			vgrid(i,j,k) = 1.0/smesh1.at(p); // override with smesh1
		    }
		}
	    }
	}
    }else if (constVert){
	if (!getBound) error("use -Cv together with -U");
	for (int i=1; i<=nx; i++){
	    for (int j=1; j<=ny; j++){
		for (int k=kstart(i,j); k<=nz; k++){
		    vgrid(i,j,k)=vgrid(i,j,k-1);
		}
	    }
	}
    }else if (paste1dprof){
	ifstream p1d(prof1dfn);

	double z1d, v1d;
	for (int k=1; k<=nz; k++){
	    p1d >> z1d >> v1d;
	    for (int i=1; i<=nx; i++){
		for (int j=1; j<=ny; j++){
		    if (z1d > zpos(kstart(i,j))){
			vgrid(i,j,k) = v1d; 
		    }
		}
	    }
	}   
    }else if (sheared1dprof){
	ifstream sp1d(sprof1dfn);
	int snz;
	Array1d<double> sv1d;
	sp1d >> snz;
	sv1d.resize(snz);
	for (int k=1; k<=snz; k++){
	    sp1d >> sv1d(k);
	}
	int ii=1;
	for (int i=1; i<=nx; i++){
	    for (int j=1; j<=ny; j++){
		for(int k=kstart(i,j); k<=nz; k++){
		    vgrid(i,j,k) = sv1d(ii++);
		}
		ii=1;
	    }
	}
    }else if (smooth){
	double Lx2 = Lx*Lx;
	double Ly2 = Ly*Ly;
	double Lz2 = Lz*Lz;
	Array3d<double> new_vgrid(nx,ny,nz);
	new_vgrid = vgrid;
	for (int i=1; i<=nx; i++){
	    double x = xpos(i);
	    for (int j=1; j<=ny; j++){
		double y = ypos(j);
		double t = topo(i,j);
		for (int k=kstart(i,j); k<=nz; k++){
		    double z = t+zpos(k);

		    if (getCorr){
			Point3d p(x,y,z);
			corr_p->at(p,Lx,Ly,Lz);
			Lx2 = Lx*Lx;
			Ly2 = Ly*Ly;
			Lz2 = Lz*Lz;
		    }
		
		    double sum = 0.0, beta_sum = 0.0;
		    for (int ii=1; ii<=nx; ii++){
			double xx = xpos(ii);
			double dx = x-xx;
			for (int jj=1; jj<=ny; jj++){
			    double yy = ypos(jj);
			    double tt = topo(ii,jj);
			    double dy = y-yy;
			    if (abs(dx)<=Lx){
				double xexp = exp(-dx*dx/Lx2);
				if (abs(dy)<=Ly){
				    double yexp = exp(-dy*dy/Ly2);
				    for (int kk=kstart(ii,jj); kk<=nz; kk++){
					double zz = tt+zpos(kk);
					double dz = z-zz;
					if (abs(dz)<=Lz){
					    double beta = xexp*yexp*exp(-dz*dz/Lz2);
					    sum += beta*vgrid(ii,jj,kk);
					    beta_sum += beta;
					}
				    }
				}
			    }
			}
		    }
		    new_vgrid(i,j,k) = sum/beta_sum;
		}
	    }
	}
	vgrid = new_vgrid;	
    }else if (checkerboard){
	double pi = 3.1415927;
	double coeffx = 2.0*pi/cx;
	double coeffy = 2.0*pi/cy;
	double coeffz = 2.0*pi/cz;
	for (int i=1; i<=nx; i++){
	    double x = xpos(i);
	    for (int j=1; j<=ny; j++){
		double y = ypos(j);
		double t = topo(i,j);
		for (int k=kstart(i,j); k<=nz; k++){
		    double z = t+zpos(k);
		    double origvel = vgrid(i,j,k);
		    double newvel = origvel*(1.0+0.01*A*sin(coeffx*x)*sin(coeffy*y)*sin(coeffz*z));
		    vgrid(i,j,k) = newvel;
		}
	    }
	}
    }else if (rectangular){
	for (int i=1; i<=nx; i++){
	    double x = xpos(i);
	    for (int j=1; j<=ny; j++){
		double y = ypos(j);
		double t = topo(i,j);
		for (int k=1; k<=nz; k++){
		    double z = t+zpos(k);
		    if (x>=dxmin && x<=dxmax && y>=dymin && y<=dymax && z>=dzmin && z<=dzmax){
			double origvel = vgrid(i,j,k);
			double newvel = origvel*(1.0+0.01*A);
			vgrid(i,j,k) = newvel;
		    }
		}
	    }
	}
    }else if (gaussian){
	double rLx2 = 1.0/(Lx*Lx);
	double rLy2 = 1.0/(Ly*Ly);
	double rLz2 = 1.0/(Lz*Lz);
	A *= 0.01; // percent to fraction
	for (int i=1; i<=nx; i++){
	    double x = xpos(i);
	    double dx = x-x0;
	    double dx2 = dx*dx;
	    for (int j=1; j<=ny; j++){
		double y = ypos(j);
		double dy = y-y0;
		double dy2 = dy*dy;
		double t = topo(i,j);
		for (int k=kstart(i,j); k<=nz; k++){
		    double z = t+zpos(k);
		    double dz = z-z0;
		    double dz2 = dz*dz;
		    double coeff = exp(-dx2*rLx2-dy2*rLy2-dz2*rLz2);
		    vgrid(i,j,k) *= (1.0+A*coeff);
		}
	    }
	}
    }else if (removeLVZ){
	for (int i=1; i<=nx; i++){
	    for (int j=1; j<=ny; j++){
		double vmax = 0.0;
		for (int k=kstart(i,j); k<=nz; k++){
		    double origvel = vgrid(i,j,k);
		    if (origvel > vmax) vmax = origvel;
		    if (origvel < vmax) vgrid(i,j,k) = vmax;
		}
	    }
	}
    }else if (randomize){
	double coeff = 1.0/RAND_MAX;
	A *= 0.01; // percent to fraction
	Array3d<double> pur(nx,ny,nz), new_pur(nx,ny,nz);
	
	srand(seed);
	int ik=0;
	for (int i=1; i<=nx; i++){
//		double x = xpos(i);
	    for (int j=1; j<=ny; j++){
//		    double y = ypos(j);
//		    double t = topo(i,j);
		for (int k=kstart(i,j); k<=nz; k++){
//			double z = t+zpos(k);
		    double amp=0;
		    if (ik%nrand==0) amp = A*2.0*(coeff*rand()-0.5);

//		if (getCorr){
//		    Point2d p(x,z);
//		    double h,v;
//		    corr_p->at(p,h,v);
//		    amp *= h*v;
//		}

		    pur(i,j,k) = amp;
		    ik++;
		}
	    }
	}
	if (getCorr){
	    for (int i=1; i<=nx; i++){
		double x = xpos(i);
		for (int j=1; j<=ny; j++){
		    double y = ypos(j);
		    double t = topo(i,j);
		    for (int k=kstart(i,j); k<=nz; k++){
			double z = t+zpos(k);
			Point3d p(x,y,z);
			corr_p->at(p,Lx,Ly,Lz);
			double Lx2 = Lx*Lx;
			double Ly2 = Ly*Ly;
			double Lz2 = Lz*Lz;
		
			double sum = 0.0, beta_sum = 0.0;
			for (int ii=1; ii<=nx; ii++){
			    double xx = xpos(ii);
			    double dx = x-xx;
			    for (int jj=1; jj<=ny; jj++){
				double yy = ypos(jj);
				double tt = topo(ii,jj);
				double dy = y-yy;
				if (abs(dx)<=Lx){
				    double xexp = exp(-dx*dx/Lx2);
				    if (abs(dy)<=Ly){
					double yexp = exp(-dy*dy/Ly2);
					for (int kk=kstart(ii,jj); kk<=nz; kk++){
					    double zz = tt+zpos(kk);
					    double dz = z-zz;
					    if (abs(dz)<=Lz){
						double beta = xexp*yexp*exp(-dz*dz/Lz2);
						sum += beta*pur(ii,jj,kk);
						beta_sum += beta;
					    }
					}
				    }
				}
			    }
			}
			new_pur(i,j,k) = sum/beta_sum;
		    }
		}
	    }
	    pur = new_pur;
	}
	for (int i=1; i<=nx; i++){
	    for (int j=1; j<=ny; j++){
		for (int k=kstart(i,j); k<=nz; k++){
		    vgrid(i,j,k) *= (1.0+pur(i,j,k));
		}
	    }
	}
    }else if (randomize2){
	double coeff = 1.0/RAND_MAX;
	A *= 0.01*2; // percent to fraction
	A = sqrt(A); 
	srand(seed);

	double pi2 = 3.1415927*2.0;
	Array1d<double> xamp, yamp, zamp, xph, yph, zph, kx, ky, kz;
	for (double x=rand_xmin; x<=rand_xmax; x+=rand_dx){
	    xamp.push_back(A*(coeff*rand()-0.5));
	    kx.push_back(pi2/x);
	    xph.push_back(pi2*(coeff*rand()-0.5));
//cerr << xamp.back() << ' ' << xph.back() << '\n';
	}
	for (double y=rand_ymin; y<=rand_ymax; y+=rand_dy){
	    yamp.push_back(A*(coeff*rand()-0.5));
	    ky.push_back(pi2/y);
	    yph.push_back(pi2*(coeff*rand()-0.5));
	}
	for (double z=rand_zmin; z<=rand_zmax; z+=rand_dz){
	    zamp.push_back(A*(coeff*rand()-0.5));
	    kz.push_back(pi2/z);
	    zph.push_back(pi2*(coeff*rand()-0.5));
//cerr << zamp.back() << ' ' << zph.back() << '\n';

	}

	for (int i=1; i<=nx; i++){
	    double x = xpos(i);
	    for (int j=1; j<=ny; j++){
		double y = ypos(j);
		double t = topo(i,j);
		for (int k=kstart(i,j); k<=nz; k++){
		    double z = t+zpos(k);
		    double purx=0.0, pury=0.0, purz=0.0;
		    for (int ii=1; ii<=xamp.size(); ii++){
			purx += xamp(ii)*sin(kx(ii)*x+xph(ii));
		    }
		    for (int jj=1; jj<=yamp.size(); jj++){
			pury += yamp(jj)*sin(ky(jj)*y+yph(jj));
		    }
		    for (int kk=1; kk<=zamp.size(); kk++){
			purz += zamp(kk)*sin(kz(kk)*z+zph(kk));
		    }
		
		    vgrid(i,j,k) *= (1.0+purx*pury*purz);
		}
	    }
	}
    }else if (randomize3){
	double coeff = 1.0/RAND_MAX;
	A *= 0.01*2; // percent to fraction
	srand(seed);

	for (int nn=1; nn<=N; nn++){
	    double Lx = rand_xmin+coeff*rand()*(rand_xmax-rand_xmin);
	    double Ly = rand_ymin+coeff*rand()*(rand_ymax-rand_ymin);		
	    double Lz = rand_zmin+coeff*rand()*(rand_zmax-rand_zmin);
	    double rLx2 = 1.0/(Lx*Lx);
	    double rLy2 = 1.0/(Ly*Ly);
	    double rLz2 = 1.0/(Lz*Lz);
	    double amp = A*(coeff*rand()-0.5);
	    int ci = int(nx*coeff*rand()); if (ci==0) ci=1;
	    int cj = int(ny*coeff*rand()); if (cj==0) cj=1;
	    int ck = int(nz*coeff*rand()); if (ck==0) ck=1;

	    for (int i=1; i<=nx; i++){
		double x = xpos(i);
		double dx = x-xpos(ci);
		double dx2 = dx*dx;
		for (int j=1; j<=ny; j++){
		    double y = ypos(j);
		    double dy = y-ypos(cj);
		    double dy2 = dy*dy;
		    double t = topo(i,j);
		    for (int k=kstart(i,j); k<=nz; k++){
			double z = t+zpos(k);
			double dz = z-zpos(ck);
			double dz2 = dz*dz;
			double gauval = amp*exp(-dx2*rLx2-dy2*rLy2-dz2*rLz2);
			vgrid(i,j,k) *= (1.0+gauval);
		    }
		}
	    }
	}
    }else if (refineMesh){
	int rnx = (nx-1)*mx+1;
	int rny = (ny-1)*my+1;
	int rnz = (nz-1)*mz+1;
	double dr = 1.0/mx;
	double ds = 1.0/my;
	double dt = 1.0/mz;
	Array1d<double> rxpos(rnx), rypos(rny), rzpos(rnz);
	Array2d<double> rtopo(rnx,rny);
	Array3d<double> rvgrid(rnx,rny,rnz);

	for (int i=1; i<nx; i++){
	    int ri=i+(i-1)*(mx-1);
	    double x1=xpos(i), x2=xpos(i+1);
	    for (int ii=0; ii<=mx; ii++){
		double r=ii*dr;
		rxpos(ri+ii) = (1-r)*x1+r*x2;
	    }
	}

	for (int j=1; j<ny; j++){
	    int rj=j+(j-1)*(my-1);
	    double y1=ypos(j), y2=ypos(j+1);
	    for (int jj=0; jj<=my; jj++){
		double s=jj*ds;
		rypos(rj+jj) = (1-s)*y1+s*y2;
	    }
	}

	for (int i=1; i<nx; i++){
	    int ri=i+(i-1)*(mx-1);
	    for (int j=1; j<ny; j++){
		int rj=j+(j-1)*(my-1);
		double t1=topo(i,j), t2=topo(i+1,j), t3=topo(i,j+1), t4=topo(i+1,j+1);
		for (int ii=0; ii<=mx; ii++){
		    double r=ii*dr;
		    for (int jj=0; jj<=my; jj++){
			double s=jj*ds;
			rtopo(ri+ii,rj+jj) = (1-r)*(1-s)*t1+r*(1-s)*t2+(1-r)*s*t3+r*s*t4;
		    }
		}
	    }
	}

	for (int k=1; k<nz; k++){
	    int rk=k+(k-1)*(mz-1);
	    double z1=zpos(k), z2=zpos(k+1);
	    for (int kk=0; kk<=mz; kk++){
		double t=kk*dt;
		rzpos(rk+kk) = (1-t)*z1+t*z2;
	    }
	}

	rvgrid = -1;
	for (int i=1; i<nx; i++){
	    int ri=i+(i-1)*(mx-1);
	    for (int j=1; j<ny; j++){
		int rj=j+(j-1)*(my-1);
		for (int k=1; k<nz; k++){
		    int rk=k+(k-1)*(mz-1);
		    double v1=vgrid(i,j,k), v2=vgrid(i+1,j,k), v3=vgrid(i,j+1,k), v4=vgrid(i,j,k+1),
			v5=vgrid(i+1,j+1,k), v6=vgrid(i+1,j,k+1), v7=vgrid(i,j+1,k+1), v8=vgrid(i+1,j+1,k+1);
		    for (int ii=0; ii<=mx; ii++){
			double r=ii*dr;
			for (int jj=0; jj<=my; jj++){
			    double s=jj*ds;
			    for (int kk=0; kk<=mz; kk++){
				double t=kk*dt;
				if (rvgrid(ri+ii,rj+jj,rk+kk)<0){
				    rvgrid(ri+ii,rj+jj,rk+kk)
					= v1*(1-r)*(1-s)*(1-t)+v2*r*(1-s)*(1-t)+v3*(1-r)*s*(1-t)+v4*(1-r)*(1-s)*t
					+v5*r*s*(1-t)+v6*r*(1-s)*t+v7*(1-r)*s*t+v8*r*s*t;
				}
			    }
			}
		    }
		}
	    }
	}
	// output smesh
	cout << rnx << " " << rny << " " << rnz << " "
	     << v_water << " " << v_air << '\n';
	for (int i=1; i<=rnx; i++) cout << rxpos(i) << " ";
	cout << '\n';
	for (int j=1; j<=rny; j++) cout << rypos(j) << " ";
	cout << '\n';
	for (int i=1; i<=rnx; i++){
	    for (int j=1; j<=rny; j++){
		cout << rtopo(i,j) << " ";
	    }
	    cout << '\n';
	}
	for (int k=1; k<=rnz; k++) cout << rzpos(k) << " ";
	cout << '\n';
	for (int i=1; i<=rnx; i++){
	    for (int j=1; j<=rny; j++){
		for (int k=1; k<=rnz; k++){
		    cout << rvgrid(i,j,k) << " ";
		}
		cout << '\n';
	    }
	}
	return 0; // end the program here.
    }else if (subMohoVel){
	for (int i=1; i<=nx; i++){
	    for (int j=1; j<=ny; j++){
		for (int k=kstart(i,j); k<=nz; k++){
		    vgrid(i,j,k) = vsubmoho;
		}
	    }
	}
    }

    // output smesh
    cout << nx << " " << ny << " " << nz << " "
	 << v_water << " " << v_air << '\n';
    for (int i=1; i<=nx; i++) cout << xpos(i) << " ";
    cout << '\n';
    for (int j=1; j<=ny; j++) cout << ypos(j) << " ";
    cout << '\n';
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    cout << topo(i,j) << " ";
	}
	cout << '\n';
    }
    for (int k=1; k<=nz; k++) cout << zpos(k) << " ";
    cout << '\n';
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    for (int k=1; k<=nz; k++){
		cout << vgrid(i,j,k) << " ";
	    }
	    cout << '\n';
	}
    }
}
