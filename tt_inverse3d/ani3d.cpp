/*
 * ani3d.cpp - anisotropy mesh implementation
 * based on smesh3d.cpp
 *
 * Adria Melendez, Alain Miniussi and Jun Korenaga
 * Winter 2016
 */
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include "ani3d.hpp"
#include "d_nr.hpp"

namespace detail {
    int
    find_interval(double x, std::vector<double> const& coords) {
        assert(coords.size() != 0);
        typedef std::vector<double>::const_iterator iter;
        iter b = coords.begin();
        iter e = coords.end();
        // All position below the first interval are considered in 
        // the first interval. As if the first position was -inf
        ++b;
        // All position beyond the last interval are considered in 
        // the last interval. As if the last position was +inf
        --e;
        int pos;
        // From the preceeding comments, if there is only one interval,
        // it's considered covering -inf,+inf
        if (coords.size()==1) {
            pos = 0;
        } else {
            iter lb = std::upper_bound(b, e, x);
            pos = lb - b;
        }
        return pos;        
    }
}

AnisotropyMesh3d::AnisotropyMesh3d(AnisotropyMesh3d const&  o)
:   Mesh3d(o),
    a_water(o.a_water),
    a_air(o.a_air), 
    agrid(o.agrid),
    ser_index(o.ser_index),
    index2cell(o.index2cell),
    node_index(o.node_index),
    cell_index(o.cell_index),
    Sm_H1(o.Sm_H1),
    Sm_H2(o.Sm_H2),
    Sm_V(o.Sm_V),
    T_common(o.T_common) {}

AnisotropyMesh3d::AnisotropyMesh3d(std::string const& fname) // Modified
:   Mesh3d()
{
    ifstream s(fname.c_str());
    if (!s){
        cerr << "AnisotropyMesh3d::cannot open " << fname << "\n";
        exit(1);
    }

    int szx, szy, szz;
    s >> szx >> szy >> szz >> a_water >> a_air;
    Mesh3d::resize(szx,szy,szz);
    Mesh3d::read_xpos(s);
    Mesh3d::read_ypos(s);
    Mesh3d::read_topo(s);
    Mesh3d::read_zpos(s);
    
    if (a_water != 0.0)
        cerr << "AnisotropyMesh3d::invalid water anisotropy. Must be set to 0.\n";
    if (a_air != 0.0)
        cerr << "AnisotropyMesh3d::invalid air anisotropy. Must be set to 0.\n";
    
    agrid.resize(nx(),ny(),nz());
    ser_index.resize(nx(),ny(),nz());
    node_index.resize(nb_nodes());
    //even with flat geometry, we need one indirection fake cell
    cell_index.resize(std::max(nx()-1,1)*std::max(ny()-1,1)*(nz()-1));
    index2cell.resize(std::max(nx()-1,1),std::max(ny()-1,1),nz()-1);

    int N=1;
    agrid = 0.0;
    for (int i=1; i<=nx(); i++){
        for(int j=1; j<=ny(); j++){
            for (int k=1; k<=nz(); k++){
                double ani;
                s >> ani;
                agrid(i,j,k) = ani;
                ser_index(i,j,k) = N;
                node_index(N).set(i,j,k);
                N++;
            }
        }
    }
    int icell=1;
    for (int i=1; i<std::max(nx(),2); i++){
        for(int j=1; j<std::max(ny(),2); j++){
            for (int k=1; k<nz(); k++){
                index2cell(i,j,k) = icell;
                cell_index[icell++-1] = ser_index(i,j,k);
            }
        }
    }

    // sanity check
    for (int i=2; i<=nx(); i++){
        if (xpos(i)<=xpos(i-1)){
            cerr << "AnisotropyMesh3d: illegal ordering of x nodes at xpos(" << i << ")="  << xpos(i)
                 << "<=" << "xpos(" << i-1 << ")=" << xpos(i-1) << '\n';
            exit(1);
        }
    }
    for (int j=2; j<=ny(); j++){
        if (ypos(j)<=ypos(j-1)){
            cerr << "AnisotropyMesh3d: illegal ordering of y nodes at ypos(" << j << ")="  << ypos(j)
                 << "<=" << "ypos(" << j-1 << ")=" << ypos(j-1) << '\n';
            exit(1);
        }
    }
    for (int i=2; i<=nz(); i++){
        if (zpos(i)<=zpos(i-1)){
            cerr << "AnisotropyMesh3d: illegal ordering of z nodes zpos(" << i << ")="  << zpos(i)
                 << "<=" << "zpos(" << i-1 << ")=" << zpos(i-1) << '\n';
            exit(1);
        }
    }
    
    commonNormKernel();
}

boost::optional<double>
AnisotropyMesh3d::xexp(Point3d const& p, int i, int j, int k, double lx) const {
    if (xflat()) {
        return boost::optional<double>(1.0);
    } else {
        int bnode = nodeIndex(i,j,k);
        double dx = nodePos(bnode).x()-p.x();
        if (abs(dx) <= lx){
            return boost::optional<double>(std::exp(-(dx*dx)/(lx*lx)));
        } else {
            return boost::optional<double>();
        }
    }
}

boost::optional<double>
AnisotropyMesh3d::yexp(Point3d const& p, int i, int j, int k, double ly) const {
    if (yflat()) {
        return boost::optional<double>(1.0);
    } else {
        int bnode = nodeIndex(i,j,k);
        double dy = nodePos(bnode).y()-p.y();
        if (abs(dy) <= ly){
            return boost::optional<double>(std::exp(-(dy*dy)/(ly*ly)));
        } else {
            return boost::optional<double>();
        }
    }
}

boost::optional<double>
AnisotropyMesh3d::zexp(Point3d const& p, int nodeIndex, double lz) const {
    double dz = nodePos(nodeIndex).z()-p.z();
    if (abs(dz) <= lz){
        return boost::optional<double>(std::exp(-(dz*dz)/(lz*lz)));
    } else {
        return boost::optional<double>();
    }
}
        
void AnisotropyMesh3d::set(const Array1d<double>& u) // Modified
{
    assert(u.size() == nb_nodes());
    //assert(a.size() == nb_nodes());

    int N=1;
    for (int i=1; i<=nx(); i++){
        for (int j=1; j<=ny(); j++){
            for (int k=1; k<=nz(); k++){
                agrid(i,j,k) = u(N);
                N++;
            }
        }
    }
}

void AnisotropyMesh3d::get(Array1d<double>& u) const // Modified
{
    std::vector<double> anisotropy_model;
    u.reserve(nb_nodes());
    for (int i=1; i<=nx(); ++i){
        for (int j=1; j<=ny(); ++j){
            for (int k=1; k<=nz(); ++k){
                anisotropy_model.push_back(agrid(i,j,k));
            }
        }
    }
    u.stl().swap(anisotropy_model);
}

void AnisotropyMesh3d::vget(Array1d<double>& v) const // Modified
{
    assert(v.size() == nb_nodes());
    //assert(a.size() == nb_nodes());

    int N=1;
    for (int i=1; i<=nx(); i++){
        for (int j=1; j<=ny(); j++){
            for (int k=1; k<=nz(); k++){
                v(N) = agrid(i,j,k);
                N++;
            }
        }
    }
}

const Index3d&
AnisotropyMesh3d::nodeIndex(int i) const { return node_index(i); }

int AnisotropyMesh3d::nodeIndex(int i, int j, int k) const { return ser_index(i,j,k); } // Modified

Point3d AnisotropyMesh3d::nodePos(int i) const // Modified
{
    int node_i = node_index(i).i();
    int node_j = node_index(i).j();
    int node_k = node_index(i).k();
    Point3d p(xpos(node_i),ypos(node_j),zpos(node_k)+topo(node_i,node_j));

    return p;
}

double AnisotropyMesh3d::at(const Point3d& pos, Index3d& guess) const
{
    upper_left(pos,guess);
    if (in_water(pos,guess)) {
        return a_water;
    }
    if (in_air(pos,guess)) {
        return a_air;
    }
    return Mesh3d::at(pos, guess, agrid);
}

double AnisotropyMesh3d::at(const Point3d& pos) const
{
    Index3d guess = nodeIndex(nearest(pos));
    return at(pos,guess);
}

double AnisotropyMesh3d::at(const Point3d& pos, Index3d& guess, Point3d& da) const
{
    upper_left(pos,guess);
    if (in_water(pos,guess)){
        da = Point3d(0,0,0);
        return a_water;
    }
    if (in_air(pos,guess)){
        da = Point3d(0,0,0);
        return a_air;
    }
    
    Point3d ratio = ratios(pos,guess);
    double a = interpolate(ratio, guess, agrid);
    double r = ratio.x();
    double s = ratio.y();
    double t = ratio.z();
    double rr = 1-r;
    double ss = 1-s;
    double tt = 1-t;
    int i = guess.i();
    int j = guess.j();
    int k = guess.k();
    
    double dadt, dadtt;
    double dadx, dady, dadz;
    
    if (xflat()) {
        dadt  = ss*agrid(i,j,k+1) + s*agrid(i,j+1,k+1);
        dadtt = ss*agrid(i,j,k)   + s*agrid(i,j+1,k);
    } else if (yflat()) {
        dadt  = rr*agrid(i,j,k+1) + r*agrid(i+1,j,k+1);
        dadtt = rr*agrid(i,j,k)   + r*agrid(i+1,j,k);
    } else {
        dadt  = (rr*ss*agrid(i,j,k+1) + rr*s*agrid(i,j+1,k+1)
                 + r*ss*agrid(i+1,j,k+1) + r*s*agrid(i+1,j+1,k+1));
        dadtt = (rr*ss*agrid(i,j,k)   + r*ss*agrid(i+1,j,k)
                 + rr*s*agrid(i,j+1,k)   + r*s*agrid(i+1,j+1,k));
    }
    
    dadz = (dadt-dadtt)/dz(k);

    if (xflat()) {
        dadx = 0;
        double dads  = tt*agrid(i,j+1,k) + t*agrid(i,j+1,k+1);
        double dadss = tt*agrid(i,j,k)   + t*agrid(i,j,k+1);
        dady = ((dads-dadss)-by(i,j)*dadz)/dy(j);
    } else if (yflat()) {
        dady = 0;
        double dadr  = tt*agrid(i+1,j,k) + t*agrid(i+1,j,k+1);
        double dadrr = tt*agrid(i,j,k)   + t*agrid(i,j,k+1);
        dadx = ((dadr-dadrr)-bx(i,j)*dadz)/dx(i);
    } else {
        double dadr  = (ss*tt*agrid(i+1,j,k) + ss*t*agrid(i+1,j,k+1) 
                        + s*tt*agrid(i+1,j+1,k) + s*t*agrid(i+1,j+1,k+1));
        double dads  = (rr*tt*agrid(i,j+1,k) + rr*t*agrid(i,j+1,k+1)
                        + r*tt*agrid(i+1,j+1,k) + r*t*agrid(i+1,j+1,k+1));
        double dadrr = (ss*tt*agrid(i,j,k)   + s*tt*agrid(i,j+1,k)
                        + ss*t*agrid(i,j,k+1)   + s*t*agrid(i,j+1,k+1));
        double dadss = (rr*tt*agrid(i,j,k)   + r*tt*agrid(i+1,j,k)
                        + rr*t*agrid(i,j,k+1)   + r*t*agrid(i+1,j,k+1));
        dadx = ((dadr-dadrr)-(s*bx(i,j)+ss*bx(i,j+1))*dadz)/dx(i);
        dady = ((dads-dadss)-(r*by(i,j)+rr*by(i+1,j))*dadz)/dy(j);
    }
    da = Point3d(dadx,dady,dadz);
    return a;
}

void
AnisotropyMesh3d::cellNodes(int icell, std::vector<int>& n) const
{
    n[0] = cell_index[icell-1];
    Index3d index=node_index(n[0]);
    int i=index.i(), j=index.j(), k=index.k();
    if (xflat()) {
        n[1] = ser_index(i,j+1,k);
        n[2] = ser_index(i,j,k+1);
        n[3] = ser_index(i,j+1,k+1);
    } else if (yflat()) {
        n[1] = ser_index(i+1,j,k);
        n[2] = ser_index(i,j,k+1);
        n[3] = ser_index(i+1,j,k+1);
    } else {
        n[1] = ser_index(i+1,j,k);
        n[2] = ser_index(i,j+1,k);
        n[3] = ser_index(i,j,k+1);
        n[4] = ser_index(i,j+1,k+1);
        n[5] = ser_index(i+1,j,k+1);
        n[6] = ser_index(i+1,j+1,k);
        n[7] = ser_index(i+1,j+1,k+1);
    }
}

void AnisotropyMesh3d::cellNormKernel(int icell, Array2d<double>& T) const // Modified
{
    Index3d index=node_index(cell_index[icell-1]);
    double x = xflat() ? 1 : dx(index.i());
    double y = yflat() ? 1 : dy(index.j());
    double z = dz(index.k());
    T = T_common*sqrt(x*y*z);
}

void AnisotropyMesh3d::commonNormKernel() // Modified
{
    Array2d<double>& T = T_common;
    if (flat()) {
        T.resize(4,4);
        T(1,1) = 4; T(1,2) = 2; T(1,3) = 1; T(1,4) = 2;
        T(2,1) = 2; T(2,2) = 4; T(2,3) = 2; T(2,4) = 1;
        T(3,1) = 1; T(3,2) = 2; T(3,3) = 4; T(3,4) = 2;
        T(4,1) = 2; T(4,2) = 1; T(4,3) = 1; T(4,4) = 4;
        T *= (1.0/36.0);
    } else {
        T.resize(8,8);
        T(1,1) = 24; T(1,2) = 4;  T(1,3) = 4;  T(1,4) = 2;  T(1,5) = 4;  T(1,6) = 2;  T(1,7) = 2;  T(1,8) = 1;
        T(2,1) = 4;  T(2,2) = 24; T(2,3) = 2;  T(2,4) = 4;  T(2,5) = 2;  T(2,6) = 4;  T(2,7) = 1;  T(2,8) = 2;
        T(3,1) = 4;  T(3,2) = 2;  T(3,3) = 24; T(3,4) = 4;  T(3,5) = 2;  T(3,6) = 1;  T(3,7) = 4;  T(3,8) = 2;
        T(4,1) = 2;  T(4,2) = 4;  T(4,3) = 4;  T(4,4) = 24; T(4,5) = 1;  T(4,6) = 2;  T(4,7) = 2;  T(4,8) = 4;
        T(5,1) = 4;  T(5,2) = 2;  T(5,3) = 2;  T(5,4) = 1;  T(5,5) = 24; T(5,6) = 4;  T(5,7) = 4;  T(5,8) = 2;
        T(6,1) = 2;  T(6,2) = 4;  T(6,3) = 1;  T(6,4) = 2;  T(6,5) = 4;  T(6,6) = 24; T(6,7) = 2;  T(6,8) = 4;
        T(7,1) = 2;  T(7,2) = 1;  T(7,3) = 4;  T(7,4) = 2;  T(7,5) = 4;  T(7,6) = 2;  T(7,7) = 24; T(7,8) = 4;
        T(8,1) = 1;  T(8,2) = 2;  T(8,3) = 2;  T(8,4) = 4;  T(8,5) = 2;  T(8,6) = 4;  T(8,7) = 4;  T(8,8) = 24;
        
        T *= (1.0/216.0);
    }

    Array1d<double> p(T.nCol());
    d_choldc(T.toRecipe(), p.size(), p.toRecipe());
    for (int i=1; i <= p.size(); i++){
        T(i,i) = p(i);
        for (int j=i+1; j<=p.size(); j++){
            T(i,j) = T(j,i);
        }
    }
}

int AnisotropyMesh3d::locate_in_cell(const Point3d& p, Index3d& guess,
                                 opt_idx& jLFU, opt_idx& jRFU, opt_idx& jLBU,
                                 opt_idx& jLFL, opt_idx& jLBL, opt_idx& jRFL,
                                 opt_idx& jRBU, opt_idx& jRBL,
                                 Point3d& r) const
{
    upper_left(p,guess);
    if (in_water(p,guess)) {
        return -1;
    }
    if (in_air(p,guess)) {
        return -2;
    }
    
    int i=guess.i(), j=guess.j(), k=guess.k();

    jLFU = ser_index(i,j,k);
    jLFL = ser_index(i,j,k+1);
    if (!xflat()) {
        jRFU = ser_index(i+1,j,k);
        jRFL = ser_index(i+1,j,k+1);
    }
    if (!yflat()) {
        jLBU = ser_index(i,j+1,k);
        jLBL = ser_index(i,j+1,k+1);
    }
    if (!(xflat() || yflat())) {
        jRBU = ser_index(i+1,j+1,k);
        jRBL = ser_index(i+1,j+1,k+1);
    }
    
    r = ratios(p, guess);
    return index2cell(i,j,k);
}

namespace detail {
    int
    near(double x, std::vector<double> const& pos) {
        if (pos.size() == 1) {
            return 1;
        } else {
            for (std::size_t i=1; i<pos.size(); i++){
                double halfway = pos[i-1] + 0.5*(pos[i] - pos[i-1]);
                if ( x < halfway ) {
                    return i;
                }
            }
            return pos.size();
        }
    }
}

int AnisotropyMesh3d::nearest(const Point3d& src) const // Modified
{
    // returns the serial index of the node nearest to
    // the given source location
    int inear = detail::near(src.x(), xpos());
    int jnear = detail::near(src.y(), ypos());
    int knear = detail::near(src.z()-topo(inear,jnear), zpos());
    
    return ser_index(inear, jnear, knear);
}

void
AnisotropyMesh3d::outMesh(ostream& os) const {
    os << nx() << " " << ny() << " " << nz() << " "
       << a_water << " " << a_air << '\n';
    
    for (int i=1; i<=nx(); i++) {
        os << xpos(i) << " ";
    }
    os << '\n';
    for (int j=1; j<=ny(); j++) {
        os << ypos(j) << " ";
    }
    os << '\n';
    for (int i=1; i<=nx(); i++){
        for (int j=1; j<=ny(); j++){
            os << topo(i,j) << " ";
        }
        os << '\n';
    }
    for (int k=1; k<=nz(); k++) {
        os << zpos(k) << " ";
    }
    os << '\n';
    for (int i=1; i<=nx(); i++){
        for (int j=1; j<=ny(); j++){
            for (int k=1; k<=nz(); k++){
                os << agrid(i,j,k) << " ";
            }
            os << '\n';
        }
    }
}

bool AnisotropyMesh3d::in_water(const Point3d& pos, const Index3d& guess) const
{
    if (pos.z()<0.0) {
        return false; // negative z means above horizon
    } else {
        int i = guess.i();
        int j = guess.j();
        double zsurf = topo(i,j);

        double rx = 0; 
        double ry = 0;
        if (!xflat()) {
            zsurf += bx(i,j)*xratio(pos,i);
        }
        if (!yflat()) {
            zsurf += by(i,j)*yratio(pos,j);
        }
        if (!(xflat() || yflat())) {
            zsurf += (bx(i,j+1)-bx(i,j))*rx*ry;
        }
        return pos.z()<zsurf;
    }
}

bool AnisotropyMesh3d::in_air(const Point3d& pos, const Index3d& guess) const // Modified
{
    if (pos.z()>=0.0) {
        return false; // positive z means below horizon
    }
    
    int i=guess.i(), j=guess.j();
    double rx = xratio(pos, i);
    double ry = yratio(pos, j);
    double zsurf = topo(i,j);
    if (!xflat()) {
        zsurf += bx(i,j)*rx;
    }
    if (!yflat()) {
        zsurf += by(i,j)*ry;
    }
    if (!(xflat() || yflat())) {
        zsurf += (bx(i,j+1)-bx(i,j))*rx*ry;
    }
    return pos.z()<zsurf;
}

bool AnisotropyMesh3d::inWater(const Point3d& pos) const // Modified
{
    return in_water(pos,upper_left(pos));
}

bool AnisotropyMesh3d::inAir(const Point3d& pos) const // Modified
{
    return in_air(pos,upper_left(pos));
}

void AnisotropyMesh3d::printElements(ostream& os) const
{
    for (int i=1; i<nx(); i++){
        for (int j=1; j<ny(); j++){
            for (int k=1; k<nz(); k++){
                os << ">\n";
                os << xpos(i) << " " << ypos(j)<< " " << zpos(k)+topo(i,j) << '\n';
                if (!xflat()) {
                    os << xpos(i+1) << " " << ypos(j) << " " << zpos(k)+topo(i+1,j) << '\n';
                }
                if (!yflat()) {
                    os << xpos(i) << " " << ypos(j+1) << " " << zpos(k)+topo(i,j+1) << '\n';
                }
                os << xpos(i) << " " << ypos(j) << " " << zpos(k+1)+topo(i,j) << '\n';
                if (!(xflat() || yflat())) {
                    os << xpos(i+1) << " " << ypos(j+1) << " " << zpos(k)+topo(i+1,j+1) << '\n';
                }
                if (!xflat()) {
                    os << xpos(i+1) << " " << ypos(j) << " " << zpos(k+1)+topo(i+1,j) << '\n';
                }
                if (!yflat()) {
                    os << xpos(i) << " " << ypos(j+1) << " " << zpos(k+1)+topo(i,j+1) << '\n';
                }
                if (!(xflat() || yflat())) {
                    os << xpos(i+1) << " " << ypos(j+1) << " " << zpos(k+1)+topo(i+1,j+1) << '\n';
                }
            }
        }
    }
}

void AnisotropyMesh3d::printAGrid(ostream& os, bool printAW) const // Modified
{
    double min_topo=0.0;
    for (int i=1; i<=nx(); i++){
        for (int j=1; j<=ny(); j++){
            if (topo(i,j)<min_topo) min_topo=topo(i,j);
        }
    }
    min_topo -= 1.0;        // upper bound for courtesy grid making
    // (extra 1km)
    double dz2 = dz(1)/2;
    for (int i=1; i<=nx(); i++){
        for (int j=1; j<=ny(); j++){
            for (int k=1; k<=nz(); k++){
                os << xpos(i) << " "
                   << ypos(j) << " "
                   << zpos(k)+topo(i,j) << " "
                   << agrid(i,j,k) << '\n';
            }
            if (printAW){
                // some extra grid points for grid file 
                double z = topo(i,j)-dz2;
                while(z>=0){
                    os << xpos(i) << " "
                       << ypos(j) << " "
                       << z << " " << a_water << '\n';
                    z -= dz2;
                }
                while(z>=min_topo){
                    os << xpos(i) << " "
                       << ypos(j) << " "
                       << z << " " << a_air << '\n';
                    z -= dz2;
                }
            }
        }
    }
}

void AnisotropyMesh3d::printAGrid(ostream& os,
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

void AnisotropyMesh3d::printMaskGrid(ostream& os,
                                   const std::vector<int>& valid_node) const // Modified, used?
{
    assert(valid_node.size() == nb_nodes());

    os << nx() << " " << ny() << " " << nz() << " "
       << a_water << " " << a_air << '\n';

    for (int i=1; i<=nx(); i++) os << xpos(i) << " ";
    os << '\n';
    for (int j=1; j<=ny(); j++) os << ypos(j) << " ";
    os << '\n';
    for (int i=1; i<=nx(); i++){
        for (int j=1; j<=ny(); j++){
            os << topo(i,j) << " ";
        }
    }
    os << '\n';
    for (int k=1; k<=nz(); k++) os << zpos(k) << " ";
    os << '\n';
    int inode=1;
    for (int i=1; i<=nx(); i++){
        for (int j=1; j<=ny(); j++){
            for (int k=1; k<=nz(); k++){
                double val=0.0;
                if (valid_node[inode-1]>0){
                    val = agrid(i,j,k);
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
void AnisotropyMesh3d::printMaskGrid(ostream& os,
                                   const Array1d<double>& dws) const // Modified
{
    assert(dws.size() == nb_nodes());

    int inode=1;
    for (int i=1; i<=nx(); i++){
        double x=xpos(i);
        for (int j=1; j<=ny(); j++){
            double y=ypos(j), t=topo(i,j);
            for (int k=1; k<=nz(); k++){
                double z=zpos(k)+t;
                os << x << " " << y << " " << z << " " << dws(inode) << "\n";
                inode++;
            }
        }
    }
}
