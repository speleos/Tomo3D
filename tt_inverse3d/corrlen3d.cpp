/*
 * corrlen3d.cc based on corrlen.cc
 * correlation length functions' implementation
 * 
 * Adria Melendez and Jun Korenaga
 * Fall 2011
 */

#include <iostream>
#include <fstream>
#include <limits>
#include <cstdlib>
#include "corrlen3d.hpp"

CorrelationLength2d::CorrelationLength2d(const char* fn)
    : xpos(), ypos(),
      xval(), yval()
{
    ifstream in(fn);
    int nx=0, ny=0;
    in >> nx >> ny;
    xpos.resize(nx); ypos.resize(ny); xval.resize(nx,ny); yval.resize(nx,ny);
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    in >> xpos[i-1] >> ypos[j-1] >> xval(i,j) >> yval(i,j);
	}
    }
}

std::pair<int,double>
locate_coord(double x, std::vector<double> const& coords) {
    if (coords.size() == 1) {
        return std::make_pair(1,0.0);
    }
    if (x <= coords.front()){
	return std::make_pair(1, 0.0);
    }else if (x >= coords.back()){
	return std::make_pair(coords.size()-1, 1.0);
    }else{
        for (int i=1; i<coords.size(); i++){
            if(x >= coords[i-1] && x < coords[i]){
                return std::make_pair(i, (x-coords[i-1])/(coords[i]-coords[i-1]));
            }
	}
    }
    std::abort();
    return std::make_pair(-1, std::numeric_limits<double>::quiet_NaN());
}

double CorrelationLength2d::atx(double x, double y) const
{
    std::pair<int,double> xloc = locate_coord(x,xpos);
    std::pair<int,double> yloc = locate_coord(y,ypos);
    int ii = xloc.first;
    int jj = yloc.first;
    double xratio = xloc.second;
    double yratio = yloc.second;
    
    bool xflat = xpos.size() == 1;
    bool yflat = ypos.size() == 1;
    bool real_2d = !xflat && !yflat;
    
    double v = xval(ii,jj)*(1-xratio-yratio+xratio*yratio);
    
    if (!xflat) {
        v += xval(ii+1,jj)*(xratio-xratio*yratio);
    }
    if (!yflat) {
        v += xval(ii,jj+1)*(yratio-xratio*yratio);
    }
    if (real_2d) {
        v += xval(ii+1,jj+1)*xratio*yratio;
    }
    return v;
}

double CorrelationLength2d::aty(double x, double y) const
{
    std::pair<int,double> xloc = locate_coord(x,xpos);
    std::pair<int,double> yloc = locate_coord(y,ypos);
    int ii = xloc.first;
    int jj = yloc.first;
    double xratio = xloc.second;
    double yratio = yloc.second;
    
    bool xflat = xpos.size() == 1;
    bool yflat = ypos.size() == 1;
    bool real_2d = !xflat && !yflat;

    double v = yval(ii,jj)*(1-xratio-yratio+xratio*yratio);
    if (!xflat) {
        v += yval(ii+1,jj)*(xratio-xratio*yratio);
    }
    if (!yflat) {
	v += yval(ii,jj+1)*(yratio-xratio*yratio);
    }
    if (real_2d) {
        v += yval(ii+1,jj+1)*xratio*yratio;
    }
    return v;
}

CorrelationLength3d::CorrelationLength3d(std::string const& fn)
{
    ifstream s(fn.c_str());
    if (!s){
	cerr << "CorrelationLength3d::cannot open " << fn << "\n";
        std::abort();
    }
    
    read_dimensions(s);
    
    read_xpos(s);
    read_ypos(s);
    read_topo(s);
    read_zpos(s);

    hxgrid.resize(nx(),ny(),nz());
    hygrid.resize(nx(),ny(),nz());
    vgrid.resize(nx(),ny(),nz());
    
    for (int i=1; i<=nx(); i++){
	for (int j=1; j<=ny(); j++){
	    for (int k=1; k<=nz(); k++){
		s >> hxgrid(i,j,k);
	    }
	}
    }
    for (int i=1; i<=nx(); i++){
	for (int j=1; j<=ny(); j++){
	    for (int k=1; k<=nz(); k++){
		s >> hygrid(i,j,k);
	    }
	}
    }
    for (int i=1; i<=nx(); i++){
	for (int j=1; j<=ny(); j++){
	    for (int k=1; k<=nz(); k++){
		s >> vgrid(i,j,k);
	    }
	}
    }
}

void CorrelationLength3d::at(const Point3d& p, double& hx, double& hy, double& v) const
{
    Index3d index = upper_left(p);
    Point3d r = ratios(p,index);
    hx = interpolate(r, index,hxgrid);    
    hy = interpolate(r, index,hygrid);
    v  = interpolate(r, index,vgrid);
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

    for (int i=1; i<=nx; i++) {
        s >> xpos[i-1];
    }
    for (int j=1; j<=ny; j++) {
        s >> ypos[j-1];
    }
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++) s >> topo(i,j);
    }
    for (int k=1; k<=nz; k++) {
        s >> zpos[k-1];
    }
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    for (int k=1; k<=nz; k++){
		s >> wgrid(i,j,k);
	    }
	}
    }

    // sanity check
    for (int i=2; i<=nx; i++){
	if (xpos[i-1]<=xpos[i-2]){
	    cerr << "DampingWeight3d: illegal ordering of x nodes at xpos(" << i << ")="  << xpos[i-1]
		 << "<=" << "xpos(" << i-1 << ")=" << xpos[i-2] << '\n';
	    exit(1);
	}
    }
    for (int j=2; j<=ny; j++){
	if (ypos[j-1]<=ypos[j-2]){
	    cerr << "DampingWeight3d: illegal ordering of y nodes at ypos(" << j << ")="  << ypos[j-1]
		 << "<=" << "ypos(" << j-1 << ")=" << ypos[j-2] << '\n';
	    exit(1);
	}
    }
    for (int i=2; i<=nz; i++){
	if (zpos[i-1]<=zpos[i-2]){
	    cerr << "DampingWeight3d: illegal ordering of z nodes zpos(" << i << ")="  << zpos[i-1]
		 << "<=" << "zpos(" << i-1 << ")=" << zpos[i-2] << '\n';
	    exit(1);
	}
    }

    rdx_vec.resize(nx-1);
    rdy_vec.resize(ny-1);
    rdz_vec.resize(nz-1);
    bx_vec.resize(nx-1,ny);
    by_vec.resize(nx,ny-1);
    for (int i=1; i<nx; i++){
	rdx_vec(i) = 1.0/(xpos[i]-xpos[i-1]);
    }
    for (int j=1; j<ny; j++){
	rdy_vec(j) = 1.0/(ypos[j]-ypos[j-1]);
    }
    for (int k=1; k<nz; k++){
	rdz_vec(k) = 1.0/(zpos[k]-zpos[k-1]);
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
  
    if (xpos[i_guess-1]<=p.x()){
	if (i_guess < nx) i_guess++;
	while (xpos[i_guess-1]<=p.x() && i_guess<nx) i_guess++;
	i_guess--;
    }else{
	while (xpos[i_guess-1]>p.x() && i_guess>1) i_guess--;
    }
  
    if (ypos[j_guess-1]<=p.y()){
	if (j_guess < ny) j_guess++;
	while (ypos[j_guess-1]<=p.y() && j_guess<ny) j_guess++;
	j_guess--;
    }else{
	while (ypos[j_guess-1]>p.y() && j_guess>1) j_guess--;
    }

    double rx = (p.x()-xpos[i_guess-1])*rdx_vec(i_guess);
    double ry = (p.y()-ypos[j_guess-1])*rdy_vec(j_guess);
    double dz = rx*bx_vec(i_guess,j_guess)+ry*by_vec(i_guess,j_guess)+topo(i_guess,j_guess);
    if (zpos[k_guess-1]+dz<=p.z()){
	if (k_guess < nz) k_guess++;
	while (zpos[k_guess-1]+dz<=p.z() && k_guess<nz) k_guess++;
	k_guess--;
    }else{
	while (zpos[k_guess-1]+dz>p.z() && k_guess>1) k_guess--;
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
    double x = pos.x()-xpos[i-1];
    double y = pos.y()-ypos[j-1];
    double z = pos.z()-zpos[k-1];
	
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

    dw = (rr*ss*tt*wgrid(i,j,k)
          +rr*ss*t*wgrid(i,j,k+1)
          +rr*s*tt*wgrid(i,j+1,k)
          +r*ss*tt*wgrid(i+1,j,k)
          +r*s*tt*wgrid(i+1,j+1,k)
          +r*ss*t*wgrid(i+1,j,k+1)
          +rr*s*t*wgrid(i,j+1,k+1)
          +r*s*t*wgrid(i+1,j+1,k+1));
}
