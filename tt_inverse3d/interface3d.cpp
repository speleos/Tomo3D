/*
 * interface3d.cc based on interface.cc
 *
 * Adria Melendez and Jun Korenaga
 * Spring 2011
 */

#include <iostream>
#include <fstream>
#include <limits>

#include "boost/thread.hpp"

#include "interface3d.hpp"

#include "error.hpp"

const double Interface3d::eps = 1e-4;

Interface3d::Interface3d(const SlownessMesh3d& smesh) //modified
{
    nxr = smesh.xpos().size();
    nyr = smesh.ypos().size();
    xpos.resize(nxr); ypos.resize(nyr);
    zpos.resize(nxr,nyr);
    xpos = smesh.xpos();
    ypos = smesh.ypos();
    zpos = smesh.topo();
    calc_slope();
}

Interface3d::Interface3d(std::string const& fn) //modified
{
    ifstream in(fn.c_str());
    if(!in){
	cerr<< "Interface3d::cannot open "<<fn<<"\n";
	exit(1);
    }

    // First option: Input interface file as an xyz file
    // Row 1: number of x and y nodes

    in >> nxr >> nyr;
    xpos.resize(nxr);
    ypos.resize(nyr);
    zpos.resize(nxr,nyr);
    nnodes=nxr*nyr;
    ncells=(nxr-1)*(nyr-1);
    ser_index_r.resize(nxr,nyr);
    node_index_r.resize(nnodes);
    cell_index_r.resize(ncells);
    index2cell_r.resize(nxr-1,nyr-1);
    int N=1;
    for(int i=1; i<=nxr; i++){
	for(int j=1; j<=nyr; j++){
	    in >> xpos[i-1] >> ypos[j-1] >> zpos(i,j);
	    ser_index_r(i,j)=N;
	    node_index_r(N).set(i,j);
	    N++;
	}
    }

    int icell=1;
    for (int i=1; i<nxr; i++){
	for(int j=1; j<nyr; j++){
	    index2cell_r(i,j) = icell;
	    cell_index_r[icell++-1] = ser_index_r(i,j);
	}
    }

    // sanity check
    for (int i=2; i<=nxr; i++){
	if (xpos[i-1]<=xpos[i-2]){
	    cerr << "Interface3d: illegal ordering of x nodes at i=" << i 
                 << ", positions: " << xpos[i-2] << ", "<< xpos[i-1] << ".\n";
	    exit(1);
	}
    }
    for (int j=2; j<=nyr; j++){
	if (ypos[j-1]<=ypos[j-1-1]){
	    cerr << "Interface3d: illegal ordering of y nodes at j=" << j
                 << ", positions: " << ypos[j-2] << ", "<< ypos[j-1] << ".\n";
	    exit(1);
	}
    }
    calc_slope();
}


const Index2d& Interface3d::nodeIndex(int i) const { return node_index_r(i); }

int Interface3d::nodeIndex(int i, int j) const { return ser_index_r(i,j); } // Added

void Interface3d::calc_slope() //modified
{
    nxr=xpos.size();
    nyr=ypos.size();  

    slope_x.resize(nxr-1,nyr);
    slope_y.resize(nxr,nyr-1);
    for (int j=1; j<=nyr; j++){
	for (int i=1; i<nxr; i++) {
            slope_x(i,j) = (zpos(i+1,j)-zpos(i,j))/(xpos[i+2]-xpos[i-1]);
        }
    }
    for (int i=1; i<=nxr; i++){
	for (int j=1; j<nyr; j++) {
            slope_y(i,j) = (zpos(i,j+1)-zpos(i,j))/(ypos[j+2]-ypos[j-1]);
        }
    }
}

double Interface3d::xmin() const { return xpos.front(); }
double Interface3d::xmax() const { return xpos.back(); }
double Interface3d::ymin() const { return ypos.front(); }
double Interface3d::ymax() const { return ypos.back(); }

namespace {
    typedef std::pair<int,double> dpos;
}
std::ostream&
operator<<(std::ostream& out, dpos const& p) {
    out << "{idx:" << p.first << ",rat:" << p.second << '}';
    return out;
}

/// \brief Compute the position of a coordinate in a sequence of discrete coordinates.
/// \param x the coordinate
/// \param discrete an increasing sequence of coordinates
/// \returns a (i,r) pair such that x \in [i, i+1[ and r the position nto the interval (0 if x on i, 1 if x on 1+i).  
dpos
position(double x, std::vector<double> const& discrete) {
    assert(discrete.size() > 0);
    bool flat = discrete.size() == 1;
    if (flat) {
        return dpos(0,0);
    } else {
        if (x <= discrete.front()) {
            return dpos(0,0);
        }
        for(int i = 0; i < (discrete.size() - 1); ++i) {
            double p0 = discrete[i];
            double p1 = discrete[i+1];
            if (x >= p0 && x < p1) {
                return dpos(i,(x-p0)/(p1-p0));
            }
        }
        return dpos(discrete.size()-2, 1);
    }
}

namespace {
    void
    check_delta(char dir, double dz_old, double dz_new, double x, double y, dpos xpos, dpos ypos) {
        static boost::mutex mutex;
        static double delta_max = 0;
        double delta = std::abs(dz_new -dz_old);
        {
            boost::lock_guard<boost::mutex> lock(mutex);
            if (delta > delta_max) {
                delta_max = delta;
            }
        }
        if (delta > 0.0001) {
            boost::lock_guard<boost::mutex> lock(mutex);
            std::cout << "\nNew dzd" << dir << " error reccord: "<< delta_max
                      << " dzd" << dir << " went from " << dz_old << " to " << dz_new 
                      << " at (" << x << ',' << y << "), that is (" << xpos << ',' << ypos << ')'
                      << std::endl;
        }
    }
}

boost::tuple<double, double>
Interface3d::dzdxy(double x, double y) const
{
    double dzdy;
    double dzdx;
    // Yes, I know, it's replicated
    if (xflat()) {
        dpos yloc = position(y, ypos);
        double dy = ypos[yloc.first + 1] - ypos[yloc.first];
        frozenx<double> zline(zpos, xpos.size(), ypos.size(), 0);
        double zj0 = zline[yloc.first];
        double zj1 = zline[yloc.first+1];
        dzdy = (zj1 - zj0)/dy;
        dzdx = 0;
    } else if (yflat()) {
        dpos xloc = position(x, xpos);
        double dx = xpos[xloc.first + 1] - xpos[xloc.first];
        frozeny<double> zline(zpos, xpos.size(), ypos.size(), 0);
        double zi0 = zline[xloc.first];
        double zi1 = zline[xloc.first+1];
        dzdx = (zi1 - zi0)/dx;
        dzdy = 0;
    } else {
        dpos xloc = position(x, xpos);
        dpos yloc = position(y, ypos);
        double dx = xpos[xloc.first + 1] - xpos[xloc.first];
        double dy = ypos[yloc.first + 1] - ypos[yloc.first];
        
        plane<double> zplane(zpos, xpos.size(), ypos.size());
        
        double zi0j0 = zplane(xloc.first,   yloc.first);
        double zi0j1 = zplane(xloc.first,   yloc.first+1, 0);
        double zi1j0 = zplane(xloc.first+1, yloc.first,   0);
        double zi1j1 = zplane(xloc.first+1, yloc.first+1, 0);
        double xratio = xloc.second;
        double yratio = yloc.second;
        
        double dxdj0 = (zi1j0 - zi0j0)/dx;
        double dxdj1 = (zi1j1 - zi0j1)/dx;
        double dydi0 = (zi0j1 - zi0j0)/dy;
        double dydi1 = (zi1j1 - zi1j0)/dy;
        double dzdx = dxdj0 + (dxdj1 - dxdj0) * yratio;
        double dzdy = dydi0 + (dydi1 - dydi0) * xratio;
        if (false) {
            check_delta('x', this->dzdx(x, y), dzdx, x, y, xloc, yloc);
            check_delta('y', this->dzdy(x, y), dzdy, x, y, xloc, yloc);
        }
    }
    return boost::make_tuple(dzdx,dzdy);
}

double 
Interface3d::z(double x, double y) const
{
    dpos xloc = position(x, xpos);
    dpos yloc = position(y, ypos);
    plane<double> zplane(zpos, xpos.size(), ypos.size());
    
    double zi0j0 = zplane(xloc.first,   yloc.first);
    double zi0j1 = zplane(xloc.first,   yloc.first+1, 0);
    double zi1j0 = zplane(xloc.first+1, yloc.first,   0);
    double zi1j1 = zplane(xloc.first+1, yloc.first+1, 0);
    double xratio = xloc.second;
    double yratio = yloc.second;

    return (zi0j0   * ( 1 - xratio - yratio + xratio*yratio)
            + zi1j0 * ( xratio - xratio*yratio)
            + zi0j1 * ( yratio - xratio*yratio)
            + zi1j1 * xratio*yratio);
}

namespace details {
    struct zfct : public std::binary_function<double,double,double> {
        zfct(Interface3d const& i) : interface(i) {}
        double operator()(double x, double y) const {
            return interface.z(x,y);
        }
        Interface3d const& interface;
    };
    
    template <class ZFCT>
    double slope(ZFCT const& fct, double x1, double x2) {
        return (fct(x2) - fct(x1)) / (x2-x1);
    }
    
    template <class ZFCT>
    double
    dzdv(double x, ZFCT const&zfct, std::vector<double> const& pos, double epsilon) {
        bool out_of_bounds = x < pos.front() || x > pos.back();
        if (out_of_bounds) {
            return 0;
        }
        int nx = pos.size();
        bool flat = nx == 1;
        if (flat) {
            return 0;
        }
        bool almost_on_last = (pos.back() - x) < epsilon;
        if (almost_on_last) {
            return slope(zfct, pos[nx-2], x);
        }
        bool almost_on_first = (x-pos.front()) < epsilon;
        if (almost_on_first) {
            return slope(zfct, x, pos[1]);
        }
        for(int i=2; i<nx; ++i){
            bool almost_on_node = fabs(x-pos[i-1]) < epsilon;
            if(almost_on_node) {
                return (slope(zfct, pos[i-2], x) 
                        + slope(zfct, x, pos[i])) / 2;
            }
        }
        for(int i=1; i<nx; ++i){
            if(x>pos[i-1] && x<pos[i]){
                return slope(zfct, pos[i-1], pos[i]);
	    }
        }
        std::abort();
    }
    
}

double Interface3d::dzdx(double x, double y) const
{
    return details::dzdv(x, std::bind2nd(details::zfct(*this), y), xpos, eps);
}

double Interface3d::dzdy(double x, double y) const
{
    return details::dzdv(y, std::bind1st(details::zfct(*this), x), ypos, eps);
}

namespace details {
    Interface3d::opt_idx 
    position(double x, std::vector<double> const& pos, double epsilon) {
        int const nx = pos.size();
        assert(nx);
        if (nx == 1) {
            return 1;
        }
        for (int i=1; i<=nx; ++i){
            if(fabs(x-pos[i-1])<epsilon){
                return std::min(i, nx-1);
            }
        }
        
        for (int i=1; i<nx; ++i){
            if (x > pos[i-1] && x < pos[i]){
                return i;
            }
        }
        
        if (x < pos.front() || x > pos.back()){ // out of x-bounds
            return Interface3d::opt_idx(boost::none);
        }
        std::abort();
    }
}

void Interface3d::locate_in_plane(double x, double y,
                                  opt_idx& i0j0, opt_idx& i1j0, opt_idx& i0j1, opt_idx& i1j1) const
{
    opt_idx i = details::position(x, xpos, eps);
    opt_idx j = details::position(y, ypos, eps);
    
    if ( i && j) {
        i0j0 = ser_index_r(*i,*j);
        if (!xflat()) {
            i1j0 = ser_index_r(*i+1,*j);
        }
        if (!yflat()) {
            i0j1 = ser_index_r(*i,*j+1);
        }
        if (!(xflat() || yflat())) {
            i1j1 = ser_index_r(*i+1,*j+1);
        }
    }
}

double Interface3d::x(int i) const
{
    Index2d node=node_index_r(i);
    int inode=node.i();
    assert(inode>0 && inode<=xpos.size());
    return xpos[inode-1];
}

double Interface3d::y(int i) const
{
    Index2d node=node_index_r(i);
    int jnode=node.j();
    assert(jnode>0 && jnode<=ypos.size());
    return ypos[jnode-1];
}

void Interface3d::cellNodes(int icell, int& j1, int& j2, int& j3, int& j4) const // Added
{
    j1 = cell_index_r[icell-1];
    Index2d index=node_index_r(j1);
    int i=index.i(), j=index.j();
    j2 = ser_index_r(i+1,j);
    j3 = ser_index_r(i,j+1);
    j4 = ser_index_r(i+1,j+1);
}

void Interface3d::cellNormKernel(int icell, double& fac) const // Added
{
    Index2d index=node_index_r(cell_index_r[icell-1]);
    int i=index.i(), j=index.j();
    double dx=xpos[i]-xpos[i-1];
    double dy=ypos[j]-ypos[j-1];
    fac = sqrt(dx*dy);
}

void Interface3d::set(const Array1d<double>& a)
{
    if (a.size() != zpos.nRow()*zpos.nCol())
	error("Interface3d::set - size mismatch");
    int Nr=1;
    for (int i=1; i<=zpos.nRow(); i++){
	for (int j=1; j<=zpos.nCol(); j++){
	    zpos(i,j) = a(Nr);
	    Nr++;
	}
    }
    calc_slope();
}

void Interface3d::get(Array1d<double>& a) const
{
    if (a.size() != zpos.nRow()*zpos.nCol())
	error("Interface3d::get - size mismatch");
    int Nr=1;
    for (int i=1; i<=zpos.nRow(); i++){
	for (int j=1; j<=zpos.nCol(); j++){
	    a(Nr) = zpos(i,j);
	    Nr++;
	}
    }
}

ostream&
operator<<(ostream& os, Interface3d const& itf)
{
    std::vector<double> const& xpos = itf.xpos;
    std::vector<double> const& ypos = itf.ypos;
    std::vector<double> const& zpos = itf.zpos;
    os << xpos.size() <<" "<< ypos.size() <<'\n';
    int const ncol = itf.zpos.nCol();
    for (int i=0; i < xpos.size(); ++i){
	for (int j=0; j < ypos.size(); ++j){
	    os << xpos[i] << " " << ypos[j] << " " << zpos[i*ncol + j] << '\n';
	}
    }
    return os;
}
    
// output DWS for reflector
void Interface3d::printMaskGridRefl(ostream& os,
                                   const Array1d<double>& dwsr) const // Modified
{
    assert(dwsr.size() == nnodes);

    int inode=1;
    for (int i=1; i<=xpos.size(); i++){
        double x=xpos[i-1];
        for (int j=1; j<=ypos.size(); j++){
            double y=ypos[j-1];
	    os << x << " " << y << " " << dwsr(inode) << "\n";
	    inode++;
	}
    }
}
