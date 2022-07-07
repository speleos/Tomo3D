/*
 * corrlen3d.cc based on corrlen.cc
 * correlation length functions' implementation
 * 
 * Adria Melendez and Jun Korenaga
 * Fall 2011
 */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include "util3d.h"
#include "corrlen3d.h"

CorrelationLength2d::CorrelationLength2d(const char* fn)
    : eps(1e-6)
{
    ifstream in(fn);
    int nx=0, ny=0;
    in >> nx >> ny;
    xpos.resize(nx); ypos.resize(ny); xval.resize(nx,ny); yval.resize(nx,ny);
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    in >> xpos(i) >> ypos(j) >> xval(i,j) >> yval(i,j);
	}
    }
}

double CorrelationLength2d::atx(double x, double y) const
{
    int ii=0, jj=0;
    double xratio=0., yratio=0.;
    int nx=xpos.size();
    int ny=ypos.size();

    if (x<xpos.front()){
	ii = 1; xratio = 0.;
    }else if (x>xpos.back()){
	ii = nx-1; xratio = 1.;
    }else{
	for (int i=1; i<=nx; i++){
	    if(fabs(x-xpos(i))<eps){
		if(i!=nx){
		    ii = i; xratio = 0.;
		    break;
		}else{
		    ii=nx-1, xratio=1.;
		    break;
		}
	    }
	}
	if(ii==0){
	    for (int i=1; i<nx; i++){
		if(x>xpos(i) && x<xpos(i+1)){
		    ii = i;
		    xratio = (x-xpos(i))/(xpos(i+1)-xpos(i));
		    break;
		}
	    }
	}
    }

    if (y<ypos.front()){
	jj = 1; yratio = 0.;
    }else if (y>ypos.back()){
	jj = ny-1; yratio = 1.;
    }else{
	for (int j=1; j<=ny; j++){
	    if(fabs(y-ypos(j))<eps){
		if(j!=ny){
		    jj=j; yratio=0.;
		    break;
		}else{
		    jj=ny-1; yratio=1.;
		    break;
		}
	    }
	}
	if(jj==0){
	    for (int j=1; j<ny; j++){
		if(y>ypos(j) && y<ypos(j+1)){
		    jj = j;
		    yratio = (y-ypos(j))/(ypos(j+1)-ypos(j));
		    break;
		}
	    }
	}
    }

    return xval(ii,jj)*(1-xratio-yratio+xratio*yratio)+xval(ii+1,jj)*(xratio-xratio*yratio)
	+xval(ii,jj+1)*(yratio-xratio*yratio)+xval(ii+1,jj+1)*xratio*yratio;

    error("CorrelationLengh2d::atx() - impossible!");
}

double CorrelationLength2d::aty(double x, double y) const
{
    int ii=0, jj=0;
    double xratio=0., yratio=0.;
    int nx=xpos.size();
    int ny=ypos.size();

    if (x<xpos.front()){
	ii = 1; xratio = 0.;
    }else if (x>xpos.back()){
	ii = nx-1; xratio = 1.;
    }else{
	for (int i=1; i<=nx; i++){
	    if(fabs(x-xpos(i))<eps){
		if(i!=nx){
		    ii = i; xratio = 0.;
		    break;
		}else{
		    ii=nx-1, xratio=1.;
		    break;
		}
	    }
	}
	if(ii==0){
	    for (int i=1; i<nx; i++){
		if(x>xpos(i) && x<xpos(i+1)){
		    ii = i;
		    xratio = (x-xpos(i))/(xpos(i+1)-xpos(i));
		    break;
		}
	    }
	}
    }

    if (y<ypos.front()){
	jj = 1; yratio = 0.;
    }else if (y>ypos.back()){
	jj = ny-1; yratio = 1.;
    }else{
	for (int j=1; j<=ny; j++){
	    if(fabs(y-ypos(j))<eps){
		if(j!=ny){
		    jj=j; yratio=0.;
		    break;
		}else{
		    jj=ny-1; yratio=1.;
		    break;
		}
	    }
	}
	if(jj==0){
	    for (int j=1; j<ny; j++){
		if(y>ypos(j) && y<ypos(j+1)){
		    jj = j;
		    yratio = (y-ypos(j))/(ypos(j+1)-ypos(j));
		    break;
		}
	    }
	}
    }

    return yval(ii,jj)*(1-xratio-yratio+xratio*yratio)+yval(ii+1,jj)*(xratio-xratio*yratio)
	+yval(ii,jj+1)*(yratio-xratio*yratio)+yval(ii+1,jj+1)*xratio*yratio;

    error("CorrelationLengh2d::aty() - impossible!");
}

CorrelationLength3d::CorrelationLength3d(const char* fn)
{
    ifstream s(fn);
    if (!s){
	cerr << "CorrelationLength3d::cannot open " << fn << "\n";
	exit(1);
    }

    s >> nx >> ny >> nz;
    xpos.resize(nx);
    ypos.resize(ny);
    topo.resize(nx,ny);
    zpos.resize(nz);
    hxgrid.resize(nx,ny,nz);
    hygrid.resize(nx,ny,nz);
    vgrid.resize(nx,ny,nz);

    for (int i=1; i<=nx; i++) s >> xpos(i);
    for (int j=1; j<=ny; j++) s >> ypos(j);
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    s >> topo(i,j);
	}
    }
    for (int k=1; k<=nz; k++) s >> zpos(k);
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    for (int k=1; k<=nz; k++){
		s >> hxgrid(i,j,k);
	    }
	}
    }
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    for (int k=1; k<=nz; k++){
		s >> hygrid(i,j,k);
	    }
	}
    }
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    for (int k=1; k<=nz; k++){
		s >> vgrid(i,j,k);
	    }
	}
    }
    
    // sanity check
    for (int i=2; i<=nx; i++){
	if (xpos(i)<=xpos(i-1)){
	    cerr << "CorrelationLength3d: illegal ordering of x nodes at xpos(" << i << ")="  << xpos(i)
		 << "<=" << "xpos(" << i-1 << ")=" << xpos(i-1) << '\n';
	    exit(1);
	}
    }
    for (int j=2; j<=ny; j++){
	if (ypos(j)<=ypos(j-1)){
	    cerr << "CorrelationLength3d: illegal ordering of y nodes at ypos(" << j << ")="  << ypos(j)
		 << "<=" << "ypos(" << j-1 << ")=" << ypos(j-1) << '\n';
	    exit(1);
	}
    }
    for (int i=2; i<=nz; i++){
	if (zpos(i)<=zpos(i-1)){
	    cerr << "CorrelationLength3d: illegal ordering of z nodes zpos(" << i << ")="  << zpos(i)
		 << "<=" << "zpos(" << i-1 << ")=" << zpos(i-1) << '\n';
	    exit(1);
	}
    }

    rdx_vec.resize(nx-1);
    rdy_vec.resize(ny-1);
    rdz_vec.resize(nz-1);
    bx_vec.resize(nx-1,ny);
    by_vec.resize(nx,ny-1);
    for (int i=1; i<nx; i++){
	rdx_vec(i) = 1.0/(xpos(i+1)-xpos(i));
    }
    for (int j=1; j<ny; j++){
	rdy_vec(j) = 1.0/(ypos(j+1)-ypos(j));
    }
    for (int k=1; k<nz; k++){
	rdz_vec(k) = 1.0/(zpos(k+1)-zpos(k));
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
}

void CorrelationLength3d::upperleft(const Point3d& p, Index3d& guess) const
{
    // note: this code guarantees that the final guess index is bounded by
    //       valid range (1...nx-1)(1...ny-1)(1...nz-1).
    int i_guess=1, j_guess=1, k_guess=1;
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

void CorrelationLength3d::calc_local(const Point3d& pos, int i, int j, int k,
				     double& r, double& s, double& t, double& rr, double& ss, double& tt) const
{
    double rDX = rdx_vec(i);
    double rDY = rdy_vec(j);
    double rDZ = rdz_vec(k);

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

void CorrelationLength3d::at(const Point3d& p, double& hx, double& hy, double& v) const
{
    Index3d index;
    upperleft(p,index);
    int i=index.i(), j=index.j(), k=index.k();
    double r,s,t,rr,ss,tt;
    calc_local(p,i,j,k,r,s,t,rr,ss,tt);

    hx=rr*ss*tt*hxgrid(i,j,k)+rr*ss*t*hxgrid(i,j,k+1)+rr*s*tt*hxgrid(i,j+1,k)+r*ss*tt*hxgrid(i+1,j,k)
	+r*s*tt*hxgrid(i+1,j+1,k)+r*ss*t*hxgrid(i+1,j,k+1)+rr*s*t*hxgrid(i,j+1,k+1)+r*s*t*hxgrid(i+1,j+1,k+1);

    hy=rr*ss*tt*hygrid(i,j,k)+rr*ss*t*hygrid(i,j,k+1)+rr*s*tt*hygrid(i,j+1,k)+r*ss*tt*hygrid(i+1,j,k)
	+r*s*tt*hygrid(i+1,j+1,k)+r*ss*t*hygrid(i+1,j,k+1)+rr*s*t*hygrid(i,j+1,k+1)+r*s*t*hygrid(i+1,j+1,k+1);

    v=rr*ss*tt*vgrid(i,j,k)+rr*ss*t*vgrid(i,j,k+1)+rr*s*tt*vgrid(i,j+1,k)+r*ss*tt*vgrid(i+1,j,k)
	+r*s*tt*vgrid(i+1,j+1,k)+r*ss*t*vgrid(i+1,j,k+1)+rr*s*t*vgrid(i,j+1,k+1)+r*s*t*vgrid(i+1,j+1,k+1);
}

DampingWeight3d::DampingWeight3d(const char* fn)
{
    ifstream s(fn);
    if (!s){
	cerr << "DampingWeight3d::cannot open " << fn << "\n";
	exit(1);
    }

    s >> nx >> ny >> nz;
    xpos.resize(nx);
    ypos.resize(ny);
    topo.resize(nx,ny);
    zpos.resize(nz);
    wgrid.resize(nx,ny,nz);

    for (int i=1; i<=nx; i++) s >> xpos(i);
    for (int j=1; j<=ny; j++) s >> ypos(j);
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++) s >> topo(i,j);
    }
    for (int k=1; k<=nz; k++) s >> zpos(k);
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    for (int k=1; k<=nz; k++){
		s >> wgrid(i,j,k);
	    }
	}
    }

    // sanity check
    for (int i=2; i<=nx; i++){
	if (xpos(i)<=xpos(i-1)){
	    cerr << "DampingWeight3d: illegal ordering of x nodes at xpos(" << i << ")="  << xpos(i)
		 << "<=" << "xpos(" << i-1 << ")=" << xpos(i-1) << '\n';
	    exit(1);
	}
    }
    for (int j=2; j<=ny; j++){
	if (ypos(j)<=ypos(j-1)){
	    cerr << "DampingWeight3d: illegal ordering of y nodes at ypos(" << j << ")="  << ypos(j)
		 << "<=" << "ypos(" << j-1 << ")=" << ypos(j-1) << '\n';
	    exit(1);
	}
    }
    for (int i=2; i<=nz; i++){
	if (zpos(i)<=zpos(i-1)){
	    cerr << "DampingWeight3d: illegal ordering of z nodes zpos(" << i << ")="  << zpos(i)
		 << "<=" << "zpos(" << i-1 << ")=" << zpos(i-1) << '\n';
	    exit(1);
	}
    }

    rdx_vec.resize(nx-1);
    rdy_vec.resize(ny-1);
    rdz_vec.resize(nz-1);
    bx_vec.resize(nx-1,ny);
    by_vec.resize(nx,ny-1);
    for (int i=1; i<nx; i++){
	rdx_vec(i) = 1.0/(xpos(i+1)-xpos(i));
    }
    for (int j=1; j<ny; j++){
	rdy_vec(j) = 1.0/(ypos(j+1)-ypos(j));
    }
    for (int k=1; k<nz; k++){
	rdz_vec(k) = 1.0/(zpos(k+1)-zpos(k));
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
}

void DampingWeight3d::upperleft(const Point3d& p, Index3d& guess) const
{
    // note: this code guarantees that the final guess index is bounded by
    //       valid range (1...nx-1)(1...ny-1)(1...nz-1).
    int i_guess=1, j_guess=1, k_guess=1;
  
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

void DampingWeight3d::calc_local(const Point3d& pos, int i, int j, int k,
				 double& r, double& s, double& t, double& rr, double& ss, double& tt) const
{
    double rDX = rdx_vec(i);
    double rDY = rdy_vec(j);
    double rDZ = rdz_vec(k);

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

void DampingWeight3d::at(const Point3d& p, double& dw) const
{
    Index3d index;
    upperleft(p,index);
    int i=index.i(), j=index.j(), k=index.k();
    double r,s,t,rr,ss,tt;
    calc_local(p,i,j,k,r,s,t,rr,ss,tt);

    dw=rr*ss*tt*wgrid(i,j,k)+rr*ss*t*wgrid(i,j,k+1)+rr*s*tt*wgrid(i,j+1,k)+r*ss*tt*wgrid(i+1,j,k)
	+r*s*tt*wgrid(i+1,j+1,k)+r*ss*t*wgrid(i+1,j,k+1)+rr*s*t*wgrid(i,j+1,k+1)+r*s*t*wgrid(i+1,j+1,k+1);
}
