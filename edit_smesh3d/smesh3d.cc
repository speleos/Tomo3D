/*
 * smesh3d.cc - slowness mesh implementation
 * based on smesh.cc
 *
 * Adria Melendez and Jun Korenaga
 * Spring 2011
 */
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "smesh3d.h"
#include "interface3d.h"
#include "d_nr.h"

SlownessMesh3d::SlownessMesh3d(const char* fn) // Modified
    : eps(1e-6)
{
    ifstream s(fn);
    if (!s){
	cerr << "SlownessMesh3d::cannot open " << fn << "\n";
	exit(1);
    }

    s >> nx >> ny >> nz >> p_water >> p_air;
    if (p_water <= 0.0)
	cerr << "SlownessMesh3d::invalid water velocity\n";
    if (p_air <= 0.0)
	cerr << "SlownessMesh3d::invalid air velocity\n";
    p_water = 1.0/p_water;
    p_air = 1.0/p_air;
    nnodes = nx*ny*nz;
    ncells = (nx-1)*(ny-1)*(nz-1);

    pgrid.resize(nx,ny,nz);
    vgrid.resize(nx,ny,nz);
    ser_index.resize(nx,ny,nz);
    node_index.resize(nnodes);
    cell_index.resize(ncells);
    index2cell.resize(nx-1,ny-1,nz-1);
    xpos.resize(nx);
    ypos.resize(ny);
    topo.resize(nx,ny);
    zpos.resize(nz);

    for (int i=1; i<=nx; i++) s >> xpos(i);
    for (int j=1; j<=ny; j++) s >> ypos(j);
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++) s >> topo(i,j);
    }
    for (int k=1; k<=nz; k++) s >> zpos(k);
    
    int N=1;
    vgrid = 0.0;
    for (int i=1; i<=nx; i++){
	for(int j=1; j<=ny; j++){
	    for (int k=1; k<=nz; k++){
		double vel;
		s >> vel;
		vgrid(i,j,k) = vel;
		pgrid(i,j,k) = 1.0/vel;
		ser_index(i,j,k) = N;
		node_index(N).set(i,j,k);
		N++;
	    }
	}
    }
    int icell=1;
    for (int i=1; i<nx; i++){
	for(int j=1; j<ny; j++){
	    for (int k=1; k<nz; k++){
		index2cell(i,j,k) = icell;
		cell_index(icell++) = ser_index(i,j,k);
	    
	    }
	}
    }

    // sanity check
    for (int i=2; i<=nx; i++){
	if (xpos(i)<=xpos(i-1)){
	    cerr << "SlownessMesh3d: illegal ordering of x nodes at xpos(" << i << ")="  << xpos(i)
		 << "<=" << "xpos(" << i-1 << ")=" << xpos(i-1) << '\n';
	    exit(1);
	}
    }
    for (int j=2; j<=ny; j++){
	if (ypos(j)<=ypos(j-1)){
	    cerr << "SlownessMesh3d: illegal ordering of y nodes at ypos(" << j << ")="  << ypos(j)
		 << "<=" << "ypos(" << j-1 << ")=" << ypos(j-1) << '\n';
	    exit(1);
	}
    }
    for (int i=2; i<=nz; i++){
	if (zpos(i)<=zpos(i-1)){
	    cerr << "SlownessMesh3d: illegal ordering of z nodes zpos(" << i << ")="  << zpos(i)
		 << "<=" << "zpos(" << i-1 << ")=" << zpos(i-1) << '\n';
	    exit(1);
	}
    }
    for (int i=1; i<=nx; i++){
	for(int j=1; j<=ny; j++){
	    for (int k=1; k<=nz; k++){
		if (vgrid(i,j,k)<=0){
		    cerr << "SlownessMesh3d: non-positive velocity v("
			 << i << "," << j << "," << k << ")=" << vgrid(i,j,k) << '\n';
		    exit(1);
		}
	    }
	}
    }

    dx_vec.resize(nx-1); rdx_vec.resize(nx-1);
    dy_vec.resize(ny-1); rdy_vec.resize(ny-1);
    dz_vec.resize(nz-1); rdz_vec.resize(nz-1);
    bx_vec.resize(nx-1,ny);
    by_vec.resize(nx,ny-1);

    for (int i=1; i<nx; i++){
	dx_vec(i) = xpos(i+1)-xpos(i);
	rdx_vec(i) = 1.0/dx_vec(i);
    }
    for (int j=1; j<ny; j++){
	dy_vec(j) = ypos(j+1)-ypos(j);
	rdy_vec(j) = 1.0/dy_vec(j);
    }
    for (int k=1; k<nz; k++){
	dz_vec(k) = zpos(k+1)-zpos(k);
	rdz_vec(k) = 1.0/dz_vec(k);
    }
    for (int j=1; j<=ny; j++){
	for (int i=1; i<nx; i++){
	    bx_vec(i,j) = topo(i+1,j)-topo(i,j);
	}
    }
    for (int i=1; i<=nx; i++){
	for (int j=1; j<ny; j++){
	    by_vec(i,j) = topo(i,j+1)-topo(i,j);
	}
    }

    T_common.resize(8,8);
    commonNormKernel();
}

void SlownessMesh3d::set(const Array1d<double>& u) // Modified
{
    if (u.size() != nnodes) error("SlownessMesh3d::set - size mismatch");

    int N=1;
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    for (int k=1; k<=nz; k++){
		pgrid(i,j,k) = u(N);
		vgrid(i,j,k) = 1.0/pgrid(i,j,k);
		N++;
	    }
	}
    }
}

void SlownessMesh3d::get(Array1d<double>& u) const // Modified
{
    if (u.size() != nnodes) error("SlownessMesh3d::get - size mismatch");

    int N=1;
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    for (int k=1; k<=nz; k++){
		u(N) = pgrid(i,j,k);
		N++;
	    }
	}
    }
}

void SlownessMesh3d::vget(Array1d<double>& v) const // Modified
{
    if (v.size() != nnodes) error("SlownessMesh3d::vget - size mismatch");

    int N=1;
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    for (int k=1; k<=nz; k++){
		v(N) = vgrid(i,j,k);
		N++;
	    }
	}
    }
}

double SlownessMesh3d::xmin() const { return xpos.front(); }

double SlownessMesh3d::xmax() const { return xpos.back(); }

double SlownessMesh3d::ymin() const { return ypos.front(); } // Modified

double SlownessMesh3d::ymax() const { return ypos.back(); } // Modified

double SlownessMesh3d::zmin() const // Modified
{
    double tmin=1e30;
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    if (topo(i,j)<tmin) tmin = topo(i,j);
	}
    }
    return tmin+zpos.front();
}

double SlownessMesh3d::zmax() const // Modified
{
    double tmax=-1e30;
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    if (topo(i,j)>tmax) tmax = topo(i,j);
	}
    }
    return tmax+zpos.back();
}

const Index3d&
SlownessMesh3d::nodeIndex(int i) const { return node_index(i); }

int SlownessMesh3d::nodeIndex(int i, int j, int k) const { return ser_index(i,j,k); } // Modified

Point3d SlownessMesh3d::nodePos(int i) const // Modified
{
    int node_i = node_index(i).i();
    int node_j = node_index(i).j();
    int node_k = node_index(i).k();
    Point3d p(xpos(node_i),ypos(node_j),zpos(node_k)+topo(node_i,node_j));

    return p;
}

double SlownessMesh3d::calc_ttime(int nsrc, int nrcv) const // Modified
{
    if (nsrc<1 || nsrc>nnodes || nrcv<1 || nrcv>nnodes)	error("SlownessMesh3d::calc_ttime - invalid input");

    Index3d src=node_index(nsrc);
    Index3d rcv=node_index(nrcv);
    int isrc=src.i(), jsrc=src.j(), ksrc=src.k();
    int ircv=rcv.i(), jrcv=rcv.j(), krcv=rcv.k();
    double xs=xpos(isrc), xr=xpos(ircv), ys=ypos(jsrc), yr=ypos(jrcv);
    double zs=topo(isrc,jsrc)+zpos(ksrc), zr=topo(ircv,jrcv)+zpos(krcv);
    Point3d sloc(xs,ys,zs);
    Point3d rloc(xr,yr,zr);
    int di = ircv-isrc;
    int dj = jrcv-jsrc;
    int dk = krcv-ksrc;

    /*
     *The following if-condition selects the plane in which we will interpolate velocity.
     *In this version 1.9, we maximize the number of interpolations.
     *We interpolate velocity at the planes perpendicular to the direction with max(abs(di),abs(dj)).
     *We may prefer to minimize interpolations by selecting min(abs(di),abs(dj)).
     *
     *We discard the z-direction because it is not possible to proceed as for both horizontal
     *directions (using proportionality) due to relief.
     */

    double ttime=0.0;
   
    if(abs(di)>abs(dj)){

	if(di<0){ //exchange src and rcv to get the first for-loop to work properly.
	    swap(xs,xr);
	    swap(ys,yr);
	    swap(zs,zr);
	    swap(isrc,ircv);
	    swap(jsrc,jrcv);
	    swap(ksrc,krcv);
	    swap(sloc,rloc);
	    di=-di;
	    dj=-dj;
	    dk=-dk;
	}

	double distx = xr-xs;
	double disty = yr-ys;
	double distz = zr-zs;
	double ratioyx = disty/distx; //Calculate proportions.
	double ratiozx = distz/distx;
	Point3d precut = sloc;
	double pave1 = pgrid(isrc,jsrc,ksrc);

	for(int i=isrc+1; i<=ircv; i++){
	    double incx = (xpos(i)-xs);
	    double ycut = ys+ratioyx*incx; //Find intersections.
	    double zcut = zs+ratiozx*incx;
	    Point3d cut(xpos(i),ycut,zcut);
	    double lng = cut.distance(precut); //Distance between previous and present intersection points.

	    int jstore=0;
	    if(dj>0){
		for(int j=jsrc; j<=jrcv; j++){
		    if(fabs(ycut-ypos(j))<eps){
			if(j!=ny){
			    jstore=j;
			    break;
			}else{
			    jstore=ny-1;
			    break;
			}
		    }
		}
		if(jstore==0){
		    for(int j=jsrc; j<=jrcv; j++){
			if(j!=ny){
			    if(ycut>ypos(j) && ycut<ypos(j+1)){
				jstore=j;
				break;
			    }
			}else{
			    cerr<<"ycut="<<ycut<<" out of bounds.\n";
			    break;
			}
		    }
		}
	    }else if(dj<0){
		for(int j=jrcv; j<=jsrc; j++){ 
		    if(fabs(ycut-ypos(j))<eps){
			if(j!=ny){
			    jstore=j;
			    break;
			}else{
			    jstore=ny-1;
			    break;
			}
		    }
		}
		if(jstore==0){
		    for(int j=jrcv; j<=jsrc; j++){
			if(j!=ny){
			    if(ycut>ypos(j) && ycut<ypos(j+1)){
				jstore=j;
				break;
			    }
			}else{
			    cerr<<"ycut="<<ycut<<" out of bounds.\n";
			    break;
			}
		    }
		}
	    }else{
		if(jsrc!=ny){
		    jstore=jsrc;
		}else{
		    jstore=ny-1;
		}
	    }

	    //Interpolate bathymetry at ycut to remove it from zcut and be able to locate it between vertical nodes.

	    double ratiob = (ycut-ypos(jstore))/(ypos(jstore+1)-ypos(jstore));
	    double bat = (1-ratiob)*topo(i,jstore)+ratiob*topo(i,jstore+1);
	    double znb = zcut-bat;
    
	    //Velocity interpolation

	    double pave2;
	    if(znb<zpos(1)){
		if(zcut>=0){
		    pave2=p_water;
		}else if(zcut<0){
		    pave2=p_air;
		}
	    }else if(znb>zpos(nz)){
		double wy = (ycut-ypos(jstore))/(ypos(jstore+1)-ypos(jstore)); 
		pave2 = wy*pgrid(i,jstore+1,nz)+(1-wy)*pgrid(i,jstore,nz);
	    }else{
		int kstore=0;
		if(dk!=0){
		    for(int k=1; k<=nz; k++){
			if(fabs(znb-zpos(k))<eps){
			    if(k!=nz){
				kstore=k;
				break;
			    }else{
				kstore=nz-1;
				break;
			    }
			}
		    }
		    if(kstore==0){
			for(int k=1; k<nz; k++){
			    if(znb>zpos(k) && znb<zpos(k+1)){
				kstore=k;
				break;
			    }
			}
		    }
		    if(kstore==0){
			cerr<<"znb="<<znb<<" out of bounds.\n";
			break;
		    }
		}else{
		    if(ksrc!=nz){
			kstore=ksrc;
		    }else{
			kstore=nz-1;
		    }
		}

		double wy = (ycut-ypos(jstore))/(ypos(jstore+1)-ypos(jstore)); 
		double wz = (znb-zpos(kstore))/(zpos(kstore+1)-zpos(kstore));
		pave2 = (1-wy-wz+wy*wz)*pgrid(i,jstore,kstore)+(wy-wy*wz)*pgrid(i,jstore+1,kstore)+
		    (wz-wy*wz)*pgrid(i,jstore,kstore+1)+wy*wz*pgrid(i,jstore+1,kstore+1);
	    }
	    double pave = (pave1+pave2)*0.5;
	    ttime += lng*pave;
	    pave1=pave2;
	    precut=cut;
	}

    }else if(abs(di)<abs(dj)){

	if(dj<0){ //exchange src and rcv to get the first for loop to work properly.
	    swap(xs,xr);
	    swap(ys,yr);
	    swap(zs,zr);
	    swap(isrc,ircv);
	    swap(jsrc,jrcv);
	    swap(ksrc,krcv);
	    swap(sloc,rloc);
	    di=-di;
	    dj=-dj;
	    dk=-dk;
	}

	double distx = xr-xs;
	double disty = yr-ys;
	double distz = zr-zs;
	double ratioxy = distx/disty;
	double ratiozy = distz/disty;
	Point3d precut = sloc;
	double pave1 = pgrid(isrc,jsrc,ksrc);

	for(int j=jsrc+1; j<=jrcv; j++){
	    double incy = (ypos(j)-ys);
	    double xcut = xs+ratioxy*incy;
	    double zcut = zs+ratiozy*incy;
	    Point3d cut(xcut,ypos(j),zcut);
	    double lng = cut.distance(precut);
     
	    int istore=0;
	    if(di>0){
		for(int i=isrc; i<=ircv; i++){
		    if(fabs(xcut-xpos(i))<eps){
			if(i!=nx){
			    istore=i;
			    break;
			}else{
			    istore=nx-1;
			    break;
			}
		    }
		}
		if(istore==0){
		    for(int i=isrc; i<=ircv; i++){
			if(i!=nx){
			    if(xcut>xpos(i) && xcut<xpos(i+1)){
				istore=i;
				break;
			    }
			}else{
			    cerr<<"xcut="<<xcut<<" out of bounds.\n";
			    break;
			}
		    }
		}
	    }else if(di<0){
		for(int i=ircv; i<=isrc; i++){
		    if(fabs(xcut-xpos(i))<eps){
			if(i!=nx){
			    istore=i;
			    break;
			}else{
			    istore=nx-1;
			    break;
			}
		    }
		}
		if(istore==0){
		    for(int i=ircv; i<=isrc; i++){
			if(i!=nx){
			    if(xcut>xpos(i) && xcut<xpos(i+1)){
				istore=i;
				break;
			    }
			}else{
			    cerr<<"xcut="<<xcut<<" out of bounds.\n";
			    break;
			}
		    }
		}
	    }else{
		if(isrc!=nx){
		    istore=isrc;
		}else{
		    istore=nx-1;
		}
	    }

	    //Interpolate bathymetry at xcut to remove it from zcut and be able to locate it between vertical nodes.

	    double ratiob = (xcut-xpos(istore))/(xpos(istore+1)-xpos(istore));
	    double bat = (1-ratiob)*topo(istore,j)+ratiob*topo(istore+1,j);
	    double znb = zcut-bat;

	    //Velocity interpolation

	    double pave2;
	    if(znb<zpos(1)){
		if(zcut>=0){
		    pave2=p_water;
		}else if(zcut<0){
		    pave2=p_air;
		}
	    }else if(znb>zpos(nz)){
		double wx = (xcut-xpos(istore))/(xpos(istore+1)-xpos(istore)); 
		pave2 = wx*pgrid(istore+1,j,nz)+(1-wx)*pgrid(istore,j,nz);
	    }else{
		int kstore=0;
		if(dk!=0){
		    for(int k=1; k<=nz; k++){
			if(fabs(znb-zpos(k))<eps){
			    if(k!=nz){
				kstore=k;
				break;
			    }else{
				kstore=nz-1;
				break;
			    }
			}
		    }
		    if(kstore==0){
			for(int k=1; k<nz; k++){
			    if(znb>zpos(k) && znb<zpos(k+1)){
				kstore=k;
				break;
			    }
			}
		    }
		    if(kstore==0){
			cerr<<"znb="<<znb<<" out of bounds.\n";
			break;
		    }
		}else{
		    if(ksrc!=nz){
			kstore=ksrc;
		    }else{
			kstore=nz-1;
		    }
		}

		double wx = (xcut-xpos(istore))/(xpos(istore+1)-xpos(istore)); 
		double wz = (znb-zpos(kstore))/(zpos(kstore+1)-zpos(kstore));
		pave2 = (1-wx-wz+wx*wz)*pgrid(istore,j,kstore)+(wx-wx*wz)*pgrid(istore+1,j,kstore)+
		    (wz-wx*wz)*pgrid(istore,j,kstore+1)+wx*wz*pgrid(istore+1,j,kstore+1);
	    }
	    double pave = (pave1+pave2)*0.5;
	    ttime += lng*pave;
	    pave1=pave2;
	    precut=cut;
	}  
    
    }else if(abs(di)==abs(dj)){

	if(di==0){ //vertical path particular case
	    if(dk<0){ //exchange src and rcv for the first loop to work.
		swap(isrc,ircv);
		swap(jsrc,jrcv);
		swap(ksrc,krcv);
	    }
      
	    double zprev = zpos(ksrc);
	    double pave1 = pgrid(isrc,jsrc,ksrc);
	    for(int k=ksrc+1; k<=krcv; k++){
		double lng = fabs(zpos(k)-zprev); //Distance between previous and present nodes.
		double pave = (pave1+pgrid(isrc,jsrc,k))*0.5;
		ttime += lng*pave;
		zprev=zpos(k);
		pave1=pgrid(isrc,jsrc,k);
	    }
	}else{
	    if(di<0){ //exchange src and rcv to get the first for loop to work properly.
		swap(xs,xr);
		swap(ys,yr);
		swap(zs,zr);
		swap(isrc,ircv);
		swap(jsrc,jrcv);
		swap(ksrc,krcv);
		swap(sloc,rloc);
		di=-di;
		dj=-dj;
		dk=-dk;
	    }
      
	    double distx = xr-xs;
	    double disty = yr-ys;
	    double distz = zr-zs;
	    double ratioyx = disty/distx; //Calculate proportions.
	    double ratiozx = distz/distx;
	    Point3d precut = sloc;
	    double pave1 = pgrid(isrc,jsrc,ksrc);
      
	    for(int i=isrc+1; i<=ircv; i++){
		double incx = (xpos(i)-xs);
		double ycut = ys+ratioyx*incx; //Find intersections.
		double zcut = zs+ratiozx*incx;
		Point3d cut(xpos(i),ycut,zcut);
		double lng = cut.distance(precut); //Distance between previous and present intersection points.
	
		int jstore=0;
		if(dj>0){
		    for(int j=jsrc; j<=jrcv; j++){
			if(fabs(ycut-ypos(j))<eps){
			    if(j!=ny){
				jstore=j;
				break;
			    }else{
				jstore=ny-1;
				break;
			    }
			}
		    }
		    if(jstore==0){
			for(int j=jsrc; j<=jrcv; j++){
			    if(j!=ny){
				if(ycut>ypos(j) && ycut<ypos(j+1)){
				    jstore=j;
				    break;
				}
			    }else{
				cerr<<"ycut="<<ycut<<" out of bounds.\n";
				break;
			    }
			}
		    }
		}else if(dj<0){
		    for(int j=jrcv; j<=jsrc; j++){ 
			if(fabs(ycut-ypos(j))<eps){
			    if(j!=ny){
				jstore=j;
				break;
			    }else{
				jstore=ny-1;
				break;
			    }
			}
		    }
		    if(jstore==0){
			for(int j=jrcv; j<=jsrc; j++){
			    if(j!=ny){
				if(ycut>ypos(j) && ycut<ypos(j+1)){
				    jstore=j;
				    break;
				}
			    }else{
				cerr<<"ycut="<<ycut<<" out of bounds.\n";
				break;
			    }
			}
		    }
		}else{
		    if(jsrc!=ny){
			jstore=jsrc;
		    }else{
			jstore=ny-1;
		    }
		}

		//Interpolate bathymetry at ycut to remove it from zcut and be able to locate it between vertical nodes.
	
		double ratiob = (ycut-ypos(jstore))/(ypos(jstore+1)-ypos(jstore));
		double bat = (1-ratiob)*topo(i,jstore)+ratiob*topo(i,jstore+1);
		double znb = zcut-bat;

		//Velocity interpolation

		double pave2;
		if(znb<zpos(1)){
		    if(zcut>=0.0){
			pave2=p_water;
		    }else if(zcut<0.0){
			pave2=p_air;
		    }
		}else if(znb>zpos(nz)){
		    double wy = (ycut-ypos(jstore))/(ypos(jstore+1)-ypos(jstore)); 
		    pave2 = wy*pgrid(i,jstore+1,nz)+(1-wy)*pgrid(i,jstore,nz);
		}else{
		    int kstore=0;
		    if(dk!=0){
			for(int k=1; k<=nz; k++){
			    if(fabs(znb-zpos(k))<eps){
				if(k!=nz){
				    kstore=k;
				    break;
				}else{
				    kstore=nz-1;
				    break;
				}
			    }
			}
			if(kstore==0){
			    for(int k=1; k<nz; k++){
				if(znb>zpos(k) && znb<zpos(k+1)){
				    kstore=k;
				    break;
				}
			    }
			}
			if(kstore==0){
			    cerr<<"znb="<<znb<<" out of bounds.\n";
			    break;
			}
		    }else{
			if(ksrc!=nz){
			    kstore=ksrc;
			}else{
			    kstore=nz-1;
			}
		    }

		    double wy = (ycut-ypos(jstore))/(ypos(jstore+1)-ypos(jstore)); 
		    double wz = (znb-zpos(kstore))/(zpos(kstore+1)-zpos(kstore));
		    pave2 = (1-wy-wz+wy*wz)*pgrid(i,jstore,kstore)+(wy-wy*wz)*pgrid(i,jstore+1,kstore)+
			(wz-wy*wz)*pgrid(i,jstore,kstore+1)+wy*wz*pgrid(i,jstore+1,kstore+1);
		}

		double pave = (pave1+pave2)*0.5;
		ttime += lng*pave;
		pave1=pave2;
		precut=cut;

	    }
	}
    }

    if (ttime<0){
	cerr << "SlownessMesh3d::calc_ttime - negative traveltime encountered for ("
	     << nsrc << ", " << nrcv << ")\n";
	cerr << "pgrid(" << isrc << "," << jsrc << "," << ksrc << ")=" << pgrid(isrc,jsrc,ksrc)
	     << ", pgrid(" << ircv << "," << jrcv << "," << krcv << ")=" << pgrid(ircv,jrcv,krcv)
	     << '\n';
	exit(1);
    }

    return ttime;
}

void SlownessMesh3d::upperleft(const Point3d& p, Index3d& guess) const // Modified
{
    // note: this code guarantees that the final guess index is bounded by
    //       valid range (1...nx-1)(1...ny-1)(1...nz-1).
    int i_guess=guess.i(), j_guess=guess.j(), k_guess=guess.k();
    if (xpos(i_guess)<=p.x()){
	if (i_guess < nx) i_guess++;
	while (xpos(i_guess)<=p.x() && i_guess<nx) i_guess++;
	i_guess--;
    }else{
	while (xpos(i_guess)>p.x() && i_guess>1) i_guess--;
    }

    if (ypos(j_guess)<=p.y()){
	if (j_guess < ny) j_guess++;
	while (ypos(j_guess)<=p.y() && j_guess<ny) j_guess++;
	j_guess--;
    }else{
	while (ypos(j_guess)>p.y() && j_guess>1) j_guess--;
    }

    double rx = (p.x()-xpos(i_guess))*rdx_vec(i_guess);
    double ry = (p.y()-ypos(j_guess))*rdy_vec(j_guess);
    double dz = rx*bx_vec(i_guess,j_guess)+ry*by_vec(i_guess,j_guess)+topo(i_guess,j_guess);

    if (zpos(k_guess)+dz<=p.z()){
	if (k_guess < nz) k_guess++;
	while (zpos(k_guess)+dz<=p.z() && k_guess<nz) k_guess++;
	k_guess--;
    }else{
	while (zpos(k_guess)+dz>p.z() && k_guess>1) k_guess--;
    }

    guess.set(i_guess,j_guess,k_guess);
}

void SlownessMesh3d::calc_local(const Point3d& pos, int i, int j, int k,
				double& r, double& s, double& t, double& rr, double& ss, double& tt) const // Modified
{
    double rDX = rdx_vec(i);
    double rDY = rdy_vec(j);
    double rDZ = rdz_vec(k);
    //  double BX = bx_vec(i,j);
    //  double BY = by_vec(i,j);

    // global coordinates with an origin at the upperleft node
    double x = pos.x()-xpos(i);
    double y = pos.y()-ypos(j);
    double z = pos.z()-zpos(k);
	
    // local coordinates
    r = x*rDX;
    s = y*rDY;
    rr=1-r;
    ss=1-s;
    double b = topo(i,j)*rr*ss+topo(i+1,j)*r*ss+topo(i,j+1)*rr*s+topo(i+1,j+1)*r*s;
    t = (z-b)*rDZ;
    tt=1-t;
}

double SlownessMesh3d::at(const Point3d& pos, Index3d& guess) const // Modified
{

    upperleft(pos,guess);
  
    if (in_water(pos,guess)) return p_water;
    if (in_air(pos,guess)) return p_air;

    int i=guess.i(), j=guess.j(), k=guess.k();
    double r,s,t,rr,ss,tt;
    calc_local(pos,i,j,k,r,s,t,rr,ss,tt);

    double u=
	rr*ss*tt*vgrid(i,j,k)+r*ss*tt*vgrid(i+1,j,k)+rr*s*tt*vgrid(i,j+1,k)+rr*ss*t*vgrid(i,j,k+1)
	+rr*s*t*vgrid(i,j+1,k+1)+r*ss*t*vgrid(i+1,j,k+1)+r*s*tt*vgrid(i+1,j+1,k)+r*s*t*vgrid(i+1,j+1,k+1);
    //     if (u<0){
    // 	cerr << "at i=" << i << " k=" << k << " r=" << r
    // 	     << " s=" << s << " x=" << pos.x() << " z=" << pos.y()
    // 	     << " grid x = " << xpos(i) << "," << xpos(i+1)
    // 	     << " grid z = " << topo(i)+zpos(k) << "," << topo(i)+zpos(k+1)
    // 	     << " " << topo(i+1)+zpos(k) << "," << topo(i+1)+zpos(k+1) << '\n';
    //     }

    return 1.0/u;
}

double SlownessMesh3d::at(const Point3d& pos) const
{
    Index3d guess = nodeIndex(nearest(pos));
    return at(pos,guess);
}

double SlownessMesh3d::at(const Point3d& pos, Index3d& guess,
			  double& dudx, double& dudy, double& dudz) const // Modified
{

    upperleft(pos,guess);
    if (in_water(pos,guess)){
	dudx=dudy=dudz=0.0;
	return p_water;
    }
    if (in_air(pos,guess)){
	dudx=dudy=dudz=0.0;
	return p_air;
    }
    
    int i=guess.i(), j=guess.j(), k=guess.k();
    double r,s,t,rr,ss,tt;
    calc_local(pos,i,j,k,r,s,t,rr,ss,tt);

    double u1 = vgrid(i,j,k);
    double u2 = vgrid(i+1,j,k);
    double u3 = vgrid(i,j+1,k);
    double u4 = vgrid(i,j,k+1);
    double u5 = vgrid(i,j+1,k+1);
    double u6 = vgrid(i+1,j,k+1);
    double u7 = vgrid(i+1,j+1,k);
    double u8 = vgrid(i+1,j+1,k+1);
    double u = rr*ss*tt*u1+r*ss*tt*u2+rr*s*tt*u3+rr*ss*t*u4+rr*s*t*u5+r*ss*t*u6+r*s*tt*u7+r*s*t*u8;

    double rDX = rdx_vec(i);
    double rDY = rdy_vec(j);
    double rDZ = rdz_vec(k);
    double BX1 = bx_vec(i,j);
    double BX2 = bx_vec(i,j+1);
    double BY1 = by_vec(i,j);
    double BY2 = by_vec(i+1,j);

    double dudr = ss*tt*u2+ss*t*u6+s*tt*u7+s*t*u8;
    double duds = rr*tt*u3+rr*t*u5+r*tt*u7+r*t*u8;
    double dudt = rr*ss*u4+rr*s*u5+r*ss*u6+r*s*u8;
    double dudrr = ss*tt*u1+s*tt*u3+ss*t*u4+s*t*u5;
    double dudss = rr*tt*u1+r*tt*u2+rr*t*u4+r*t*u6;
    double dudtt = rr*ss*u1+r*ss*u2+rr*s*u3+r*s*u7;
 
    // before here we've worked with velocity but we need the spatial partial derivatives of slowness: factor -1/(u*u)
    dudz = -rDZ*(dudt-dudtt)/(u*u);
    dudx = -rDX*((dudr-dudrr)/(u*u)-(s*BX1+ss*BX2)*dudz);
    dudy = -rDY*((duds-dudss)/(u*u)-(r*BY1+rr*BY2)*dudz);

    return 1.0/u;
}

void SlownessMesh3d::cellNodes(int icell, int& j1, int& j2, int& j3, int& j4,
			       int& j5, int& j6, int& j7, int& j8) const // Modified
{
    j1 = cell_index(icell);
    Index3d index=node_index(j1);
    int i=index.i(), j=index.j(), k=index.k();
    j2 = ser_index(i+1,j,k);
    j3 = ser_index(i,j+1,k);
    j4 = ser_index(i,j,k+1);
    j5 = ser_index(i,j+1,k+1);
    j6 = ser_index(i+1,j,k+1);
    j7 = ser_index(i+1,j+1,k);
    j8 = ser_index(i+1,j+1,k+1);
}

void SlownessMesh3d::cellNormKernel(int icell, Array2d<double>& T) const // Modified
{
    Index3d index=node_index(cell_index(icell));
    int i=index.i(), j=index.j(), k=index.k();
    double dx=dx_vec(i), dy=dy_vec(j), dz=dz_vec(k);
    T = T_common*sqrt(dx*dy*dz);
}

void SlownessMesh3d::commonNormKernel() // Modified
{
    T_common(1,1) = 24; T_common(1,2) = 4; T_common(1,3) = 4; T_common(1,4) = 2;
    T_common(2,1) = 4; T_common(2,2) = 24; T_common(2,3) = 2; T_common(2,4) = 4;
    T_common(3,1) = 4; T_common(3,2) = 2; T_common(3,3) = 24; T_common(3,4) = 4;
    T_common(4,1) = 2; T_common(4,2) = 4; T_common(4,3) = 4; T_common(4,4) = 24;
    T_common(5,1) = 4; T_common(5,2) = 2; T_common(5,3) = 2; T_common(5,4) = 1;
    T_common(6,1) = 2; T_common(6,2) = 4; T_common(6,3) = 1; T_common(6,4) = 2;
    T_common(7,1) = 2; T_common(7,2) = 1; T_common(7,3) = 4; T_common(7,4) = 2;
    T_common(8,1) = 1; T_common(8,2) = 2; T_common(8,3) = 2; T_common(8,4) = 4;
    T_common(1,5) = 4; T_common(1,6) = 2; T_common(1,7) = 2; T_common(1,8) = 1;
    T_common(2,5) = 2; T_common(2,6) = 4; T_common(2,7) = 1; T_common(2,8) = 2;
    T_common(3,5) = 2; T_common(3,6) = 1; T_common(3,7) = 4; T_common(3,8) = 2;
    T_common(4,5) = 1; T_common(4,6) = 2; T_common(4,7) = 2; T_common(4,8) = 4;
    T_common(5,5) = 24; T_common(5,6) = 4; T_common(5,7) = 4; T_common(5,8) = 2;
    T_common(6,5) = 4; T_common(6,6) = 24; T_common(6,7) = 2; T_common(6,8) = 4;
    T_common(7,5) = 4; T_common(7,6) = 2; T_common(7,7) = 24; T_common(7,8) = 4;
    T_common(8,5) = 2; T_common(8,6) = 4; T_common(8,7) = 4; T_common(8,8) = 24;
    T_common *= (1.0/216.0);

    Array1d<double> p(8);
    d_choldc(T_common.toRecipe(), 8, p.toRecipe());
    for (int i=1; i<=8; i++){
	T_common(i,i) = p(i);
	for (int j=i+1; j<=8; j++){
	    T_common(i,j) = T_common(j,i);
	}
    }
}

int SlownessMesh3d::locateInCell(const Point3d& p, Index3d& guess,
				 int& j1, int& j2, int& j3 , int& j4, int& j5, int& j6, int& j7, int& j8,
				 double& r, double& s, double& t, double& rr, double& ss, double& tt) const // Modified
{
    upperleft(p,guess);
    if (in_water(p,guess)) return -1;
    if (in_air(p,guess)) return -2;
    
    int i=guess.i(), j=guess.j(), k=guess.k();
    j1 = ser_index(i,j,k);
    j2 = ser_index(i+1,j,k);
    j3 = ser_index(i,j+1,k);
    j4 = ser_index(i,j,k+1);
    j5 = ser_index(i,j+1,k+1);
    j6 = ser_index(i+1,j,k+1);
    j7 = ser_index(i+1,j+1,k);
    j8 = ser_index(i+1,j+1,k+1);
    calc_local(p,i,j,k,r,s,t,rr,ss,tt);
    return index2cell(i,j,k);
}

void SlownessMesh3d::nodalCellVolume(Array1d<double>& dx, Array1d<double>& dy, Array1d<double>& dz,
				     Array1d<Point3d>& center) const // Modified (used in gravity part)
{
    //
    // determine cell dimensions and center positions
    //
    dx.resize(nnodes);
    dy.resize(nnodes);
    dz.resize(nnodes);
    center.resize(nnodes);

    // first x & y
    for (int k=1; k<=nz; k++){
	center(ser_index(1,1,k)).x(xpos(1)+0.25*dx_vec(1));
	center(ser_index(1,1,k)).y(ypos(1)+0.25*dy_vec(1));
	center(ser_index(nx,1,k)).x(xpos(nx)-0.25*dx_vec(nx-1));
	center(ser_index(nx,1,k)).y(ypos(1)-0.25*dy_vec(1));
	center(ser_index(1,ny,k)).x(xpos(1)+0.25*dx_vec(1));
	center(ser_index(1,ny,k)).y(ypos(ny)-0.25*dy_vec(ny-1));
	center(ser_index(nx,ny,k)).x(xpos(nx)-0.25*dx_vec(nx-1));
	center(ser_index(nx,ny,k)).y(ypos(ny)-0.25*dy_vec(ny-1));
	dx(ser_index(1,1,k)) = 0.5*dx_vec(1);
	dy(ser_index(1,1,k)) = 0.5*dy_vec(1);
	dx(ser_index(nx,1,k)) = 0.5*dx_vec(nx-1);
	dy(ser_index(nx,1,k)) = 0.5*dy_vec(1);
	dx(ser_index(1,ny,k)) = 0.5*dx_vec(1);
	dy(ser_index(1,ny,k)) = 0.5*dy_vec(ny-1);
	dx(ser_index(nx,ny,k)) = 0.5*dx_vec(nx-1);
	dy(ser_index(nx,ny,k)) = 0.5*dy_vec(ny-1);
	for (int i=2; i<nx; i++){
	    for(int j=2; j<ny; j++){
		center(ser_index(i,j,k)).x(xpos(i)+0.5*(dx_vec(i)-dx_vec(i-1)));
		center(ser_index(i,j,k)).y(ypos(j)+0.5*(dy_vec(j)-dy_vec(j-1)));
		dx(ser_index(i,j,k)) = 0.5*(dx_vec(i)+dx_vec(i-1));
		dy(ser_index(i,j,k)) = 0.5*(dy_vec(j)+dy_vec(j-1));
	    }
	}
    }

    // then z
    for (int i=1; i<=nx; i++){
	for(int j=1; j<=ny; j++){
	    center(ser_index(i,j,1)).z(topo(i,j)+zpos(1)+0.25*dz_vec(1)); // top
	    dz(ser_index(i,j,1)) = 0.5*dz_vec(1);
	    center(ser_index(i,j,nz)).z(topo(i,j)+zpos(nz)-0.25*dz_vec(nz-1)); // bottom
	    dz(ser_index(i,j,nz)) = 0.5*dz_vec(nz-1);
	    for (int k=2; k<nz; k++){
		center(ser_index(i,j,k)).z(topo(i,j)+zpos(k)+0.5*(dz_vec(k)-dz_vec(k-1)));
		dz(ser_index(i,j,k)) = 0.5*(dz_vec(k)+dz_vec(k-1));
	    }
	}
    }
}

int SlownessMesh3d::nearest(const Point3d& src) const // Modified
{
    // returns the serial index of the node nearest to
    // the given source location
    int inear=-1, jnear=-1, knear=-1;
    double srcx=src.x(), srcy=src.y(), srcz=src.z();
  
    for (int i=1; i<nx; i++){
	double xnode=xpos(i);
	double xnode2=xpos(i+1);
	double hdx=0.5*(xnode2-xnode);
    
	if (i==1 && (srcx-xnode)<=eps){ // to the left
	    inear = 1;
	}else if (i==(nx-1) && (xnode2-srcx)<=eps){ // to the right
	    inear = nx;
	}else if (fabs(srcx-xnode)<eps){ // on a grid
	    inear = i; 
	}else if (fabs(srcx-xnode2)<eps){ // on a grid
	    inear = i+1; 
	}else if (srcx > xnode && srcx < xnode2){
	    if (fabs(srcx-xnode) < hdx){
		inear = i;
	    }else{
		inear = i+1;
	    }
	}
	if (inear > 0) break;
    }

    for (int j=1; j<ny; j++){
	double ynode=ypos(j);
	double ynode2=ypos(j+1);
	double hdy=0.5*(ynode2-ynode);
    
	if (j==1 && (srcy-ynode)<=eps){ // to the left
	    jnear = 1;
	}else if (j==(ny-1) && (ynode2-srcy)<=eps){ // to the right
	    jnear = ny;
	}else if (fabs(srcy-ynode)<eps){ // on a grid
	    jnear = j; 
	}else if (fabs(srcy-ynode2)<eps){ // on a grid
	    jnear = j+1; 
	}else if (srcy > ynode && srcy < ynode2){
	    if (fabs(srcy-ynode) < hdy){
		jnear = j;
	    }else{
		jnear = j+1;
	    }
	}
	if (jnear > 0) break;
    }
    
    for (int k=1; k<nz; k++){
	double znode=zpos(k)+topo(inear,jnear);
	double znode2=zpos(k+1)+topo(inear,jnear);
	double hdz=0.5*(znode2-znode);
    
	if (k==1 && (srcz-znode)<=eps){ // above the domain
	    knear = 1;
	}else if (k==(nz-1) && (znode2-srcz)<=eps){ // below the domain
	    knear = nz;
	}else if (fabs(srcz-znode)<eps){ // on a grid
	    knear = k; 
	}else if (fabs(srcz-znode2)<eps){ // on a grid
	    knear = k+1; 
	}else if (srcz > znode && srcz < znode2){
	    if (fabs(srcz-znode) < hdz){
		knear = k;
	    }else{
		knear = k+1;
	    }
	}
	if (knear > 0) break;
    }
    return ser_index(inear, jnear, knear);
}
void SlownessMesh3d::nearest(const Interface3d& itf, Array1d<int>& inodes) const // Modified
{
    int n=1;
    int nx = xpos.size();
    int ny = ypos.size();
    inodes.resize(nx*ny);

    for (int i=1; i<=nx; i++){
	for(int j=1; j<=ny; j++){
	    double x=xpos(i);
	    double y=ypos(j);
	    double z=itf.z(x,y);
	    Point3d p(x,y,z);
	    inodes(n) = nearest(p);
	    n++;
	}
    }
}

void SlownessMesh3d::outMesh(ostream& os) const // Modified
{
    os << nx << " " << ny << " " << nz << " "
       << 1.0/p_water << " " << 1.0/p_air << '\n';

    for (int i=1; i<=nx; i++) os << xpos(i) << " ";
    os << '\n';
    for (int j=1; j<=ny; j++) os << ypos(j) << " ";
    os << '\n';
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    os << topo(i,j) << " ";
	}
	os << '\n';
    }
    for (int k=1; k<=nz; k++) os << zpos(k) << " ";
    os << '\n';
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    for (int k=1; k<=nz; k++){
		//os.precision(20);
		os << 1.0/pgrid(i,j,k) << " ";
	    }
	    os << '\n';
	}
//    os << '\n';
    }
}

bool SlownessMesh3d::in_water(const Point3d& pos, const Index3d& guess) const // Modified
{
    if (pos.z()<0.0) return false; // negative z means above horizon
    
    int i=guess.i(), j=guess.j();
    double rx = (pos.x()-xpos(i))/dx_vec(i);
    double ry = (pos.y()-ypos(j))/dy_vec(j);
    double zsurf = topo(i,j)+bx_vec(i,j)*rx+by_vec(i,j)*ry+(bx_vec(i,j+1)-bx_vec(i,j))*rx*ry;
    if (pos.z()<zsurf){
	return true;
    }else{
	return false;
    }
}

bool SlownessMesh3d::in_air(const Point3d& pos, const Index3d& guess) const // Modified
{
    if (pos.z()>=0.0) return false; // positive z means below horizon
    
    int i=guess.i(), j=guess.j();
    double rx = (pos.x()-xpos(i))/dx_vec(i);
    double ry = (pos.y()-ypos(j))/dy_vec(j);
    double zsurf = topo(i,j)+bx_vec(i,j)*rx+by_vec(i,j)*ry+(bx_vec(i,j+1)-bx_vec(i,j))*rx*ry;
    if (pos.z()<zsurf){
	return true;
    }else{
	return false;
    }
}

bool SlownessMesh3d::inWater(const Point3d& pos) const // Modified
{
    Index3d guess(1,1,1);
    upperleft(pos,guess);
    return in_water(pos,guess);
}

bool SlownessMesh3d::inAir(const Point3d& pos) const // Modified
{
    Index3d guess(1,1,1);
    upperleft(pos,guess);
    return in_air(pos,guess);
}

void SlownessMesh3d::printElements(ostream& os) const // Modified
{
    for (int i=1; i<nx; i++){
	for (int j=1; j<ny; j++){
	    for (int k=1; k<nz; k++){
		os << ">\n";
		os << xpos(i) << " " << ypos(j)<< " " << zpos(k)+topo(i,j) << '\n';
		os << xpos(i+1) << " " << ypos(j) << " " << zpos(k)+topo(i+1,j) << '\n';
		os << xpos(i) << " " << ypos(j+1) << " " << zpos(k)+topo(i,j+1) << '\n';
		os << xpos(i) << " " << ypos(j) << " " << zpos(k+1)+topo(i,j) << '\n';
		os << xpos(i+1) << " " << ypos(j+1) << " " << zpos(k)+topo(i+1,j+1) << '\n';
		os << xpos(i+1) << " " << ypos(j) << " " << zpos(k+1)+topo(i+1,j) << '\n';
		os << xpos(i) << " " << ypos(j+1) << " " << zpos(k+1)+topo(i,j+1) << '\n';
		os << xpos(i+1) << " " << ypos(j+1) << " " << zpos(k+1)+topo(i+1,j+1) << '\n';
	    }
	}
    }
}

void SlownessMesh3d::printVGrid(ostream& os, bool printAW) const // Modified
{
    double min_topo=0.0;
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    if (topo(i,j)<min_topo) min_topo=topo(i,j);
	}
    }
    min_topo -= 1.0;		// upper bound for courtesy grid making
				// (extra 1km)
    double dz=dz_vec(1)*0.5;
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    for (int k=1; k<=nz; k++){
		os << xpos(i) << " "
		   << ypos(j) << " "
		   << zpos(k)+topo(i,j) << " "
		   << 1.0/pgrid(i,j,k) << '\n';
	    }
	    if (printAW){
		// some extra grid points for grid file 
		double z = topo(i,j)-dz;
		while(z>=0){
		    os << xpos(i) << " "
		       << ypos(j) << " "
		       << z << " " << 1.0/p_water << '\n';
		    z-=dz;
		}
		while(z>=min_topo){
		    os << xpos(i) << " "
		       << ypos(j) << " "
		       << z << " " << 1.0/p_air << '\n';
		    z-=dz;
		}
	    }
	}
    }
}

void SlownessMesh3d::printVGrid(ostream& os,
				double x0, double x1,
				double y0, double y1,
				double z0, double z1,
				double dx, double dy, double dz) const // Modified
{
    Index3d guess = nodeIndex(nearest(Point3d(x0,y0,z0)));

    int nxx = int((x1-x0)/dx+1);
    int nyy = int((y1-y0)/dy+1);
    int nzz = int((z1-z0)/dz+1);

    for (int ix=1; ix<=nxx; ix++){
	double x = x0+(ix-1)*dx;
	for (int iy=1; iy<=nyy; iy++){
	    double y = y0+(iy-1)*dy;
	    for (int iz=1; iz<=nzz; iz++){
		double z = z0+(iz-1)*dz;
		double p = at(Point3d(x,y,z),guess);
		os << x << " " << y << " " << z << " " << 1.0/p << '\n';
	    }
	}
    }
}

void SlownessMesh3d::printMaskGrid(ostream& os,
				   const Array1d<int>& valid_node) const // Modified, used?
{
    if (valid_node.size() != nnodes)
	error("SlownessMesh3d::printMaskGrid - size mismatch");

    os << nx << " " << ny << " " << nz << " "
       << 1.0/p_water << " " << 1.0/p_air << '\n';

    for (int i=1; i<=nx; i++) os << xpos(i) << " ";
    os << '\n';
    for (int j=1; j<=ny; j++) os << ypos(j) << " ";
    os << '\n';
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    os << topo(i,j) << " ";
	}
    }
    os << '\n';
    for (int k=1; k<=nz; k++) os << zpos(k) << " ";
    os << '\n';
    int inode=1;
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    for (int k=1; k<=nz; k++){
		double val=0.0;
		if (valid_node(inode)>0){
		    val = 1.0/pgrid(i,j,k);
		}
		os << val << " ";
		inode++;
	    }
	    os << '\n';
	}
//    os << '\n';
    }
}

// output DWS
void SlownessMesh3d::printMaskGrid(ostream& os,
				   const Array1d<double>& dws) const // Modified
{
    if (dws.size() != nnodes)
	error("SlownessMesh3d::printMaskGrid - size mismatch");

    int inode=1;
    for (int i=1; i<=nx; i++){
	double x=xpos(i);
	for (int j=1; j<=ny; j++){
	    double y=ypos(j), t=topo(i,j);
	    for (int k=1; k<=nz; k++){
		double z=zpos(k)+t;
		os << x << " " << y << " " << z << " " << dws(inode) << "\n";
		inode++;
	    }
	}
    }
}
