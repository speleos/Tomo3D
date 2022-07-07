#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <limits>

#include "mesh.hpp"

Mesh3d::Mesh3d() {}

Mesh3d::Mesh3d(Mesh3d const& m)
:   my_xpos(m.my_xpos),
    my_ypos(m.my_ypos),
    my_zpos(m.my_zpos),
    my_topo(m.my_topo) {}

Mesh3d::Mesh3d( Array1d<double>& x, Array1d<double>& y, Array1d<double>& z,
                Array2d<double>& topo) 
:   my_xpos(x),
    my_ypos(y),
    my_zpos(z),
    my_topo(topo) {}

bool
Mesh3d::read_dimensions(std::istream& in) {
    if (!in.good()) {
        std::abort();
    }
    std::string line;
    std::getline(in,line);
    std::vector<int> dims;
    std::istringstream is(line);
    std::copy(std::istream_iterator<double>(is), std::istream_iterator<double>(),
              std::back_inserter(dims));
    if (dims.size() == 2) {
        dims.insert(dims.begin(), 1);
    }
    if (dims.size() != 3) {
        std::abort();
    }
    resize(dims[0], dims[1], dims[2]);
    return in.good();
}

void
Mesh3d::resize(int nx, int ny, int nz) {
    my_xpos.resize(nx);
    my_ypos.resize(ny);
    my_zpos.resize(nz);
    my_topo.resize(nx,ny);
}

bool
Mesh3d::read(std::istream& is, Array1d<double>& data) {
    for(std::size_t i = 0; i < data.stl().size(); ++i) {
        is >> (data.stl()[i]);
        assert(is.good());
    }
    return is.good();
}

bool
Mesh3d::read_xpos(std::istream& is) {
    assert(!my_xpos.empty());
    return read(is, my_xpos);
}

bool
Mesh3d::read_ypos(std::istream& is) {
    assert(!my_ypos.empty());
    return read(is, my_ypos);
}

bool
Mesh3d::read_zpos(std::istream& is) {
    assert(!my_zpos.empty());
    return read(is, my_zpos);
}

bool
Mesh3d::read_topo(std::istream& is) {
    std::vector<double>& v = my_topo;
    assert(!v.empty());
    assert(v.size() == my_xpos.size() * my_ypos.size());
    for (std::size_t i = 0; i < v.size(); ++i) {
        is >> v[i];
        assert(is.good());
    }
    return is.good();
}

namespace {
    inline
    double
    delta_coords( std::vector<double> const& coords, int i) {
        assert(i >=0 && i < (coords.size()-1));
        return coords[i+1] - coords[i];
    }
    
    static bool const normalize_ratios = true;

    inline
    double
    coord_ratio( double x, std::vector<double> const& coords, int i) {
        double r = ( coords.size() == 1
                     ? 0
                     : (x - coords[i])/(delta_coords(coords, i)));
        if (normalize_ratios) {
            if (r < 0) { r = 0; }
            if (r > 1) { r = 1; }
        }
        return r;
    }
}

double
Mesh3d::bx(int i, int j) const {
    return my_topo(i+1,j) - my_topo(i,j);
}

double
Mesh3d::by(int i, int j) const {
    return my_topo(i,j+1) - my_topo(i,j);
}

double
Mesh3d::dx(int i) const {
    assert(!xflat());
    return delta_coords(my_xpos, i-1);
}

double
Mesh3d::dy(int i) const {
    assert(!yflat());
    return delta_coords(my_ypos, i-1);
}

double
Mesh3d::dz(int i) const {
    return delta_coords(my_zpos, i-1);
}

double
Mesh3d::xratio(double x, int i) const {
    return coord_ratio(x, my_xpos.stl(), i-1);
}
        
double
Mesh3d::yratio(double y, int j) const {
    return coord_ratio(y, my_ypos.stl(), j-1);
}

double
Mesh3d::zratio(double z, int k) const {
    return coord_ratio(z, my_zpos.stl(), k-1);
}

double 
Mesh3d::xmin() const { return my_xpos.front(); }

double
Mesh3d::xmax() const { return my_xpos.back(); }

double
Mesh3d::ymin() const { return my_ypos.front(); }

double
Mesh3d::ymax() const { return my_ypos.back(); }

double
Mesh3d::zmin() const
{
    double tmin = std::numeric_limits<double>::max();
    int nx = my_xpos.size();
    int ny = my_ypos.size();
    for (int i=1; i<= nx; i++){
	for (int j=1; j<=ny; j++){
	    if (my_topo(i,j) < tmin) {
                tmin = my_topo(i,j);
            }
        }
    }
    return tmin+my_zpos.front();
}

double 
Mesh3d::zmax() const
{
    double tmax = std::numeric_limits<double>::min();
    int nx = my_xpos.size();
    int ny = my_ypos.size();
    for (int i=1; i<=nx; i++){
	for (int j=1; j<=ny; j++){
	    if (my_topo(i,j) > tmax) {
                tmax = my_topo(i,j);
            }
        }
    }
    return tmax+my_zpos.back();
}


namespace {
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

void
Mesh3d::upper_left(Point3d const& p, Index3d& guess) const 
{
    // note: this code guarantees that the final guess index is bounded by
    //       valid range (1...nx-1)(1...ny-1)(1...nz-1).
    int i = find_interval(p.x(), my_xpos.stl())+1;
    int j = find_interval(p.y(), my_ypos.stl())+1;
    
    double dz = my_topo(i,j);
    if (!xflat()) {
        dz += xratio(p, i) * bx(i,j);
    }
    if (!yflat()) {
        dz += yratio(p, j) * by(i,j);
    }
    double z = p.z() - dz;
    int k = find_interval(z, my_zpos.stl())+1;
    guess.set(i,j,k);
}

Index3d
Mesh3d::upper_left(Point3d const& p) const
{
    // note: this code guarantees that the final guess index is bounded by
    //       valid range (1...nx-1)(1...ny-1)(1...nz-1).
    Index3d idx;
    upper_left(p,idx);
    return idx;
}

Point3d
Mesh3d::ratios(const Point3d& pos, Index3d const& idx)  const {
    assert(my_zpos.size() > 1);
    // no interpolation on flat dimensions;
    double r = xratio(pos, idx.i());
    double s = yratio(pos, idx.j());
    
    double rr = 1-r;
    double ss = 1-s;
    
    double b = my_topo(idx.i(), idx.j())*rr*ss;
    if (!xflat()) {
        b += my_topo(idx.i()+1, idx.j())*r*ss;
    }
    if (!yflat()) {
        b += my_topo(idx.i(), idx.j()+1)*rr*s;
    }
    if (!(xflat() || yflat())) {
        b += my_topo(idx.i()+1, idx.j()+1)*r*s;
    }
    double z = pos.z()-my_zpos(idx.k());
    double t = (z-b)/dz(idx.k());
//    BOOST_ASSERT_MSG(std::abs(t) < 1e5,t); //ALOUREIRO assert with debug msg
//    cout<< "Point3d (z): "<< z << "Point3d (t): "<< t <<"\n";
    assert(std::abs(t) < 1e5);
    return Point3d(r,s,t);
}

double
Mesh3d::interpolate(Point3d const& r, Index3d const& ul, Array3d<double> const& nodes) const {
    Point3d rr = Point3d(1,1,1) - r;
    int i = ul.i();
    int j = ul.j();
    int k = ul.k();
    
    double u = (rr.x()*rr.y()*rr.z()*nodes(i,j,k) 
                + rr.x()*rr.y()*r.z()*nodes(i,j,k+1));
    if (!xflat()) {
        u += (r.x()*rr.y()*rr.z()*nodes(i+1,j,k) 
              + r.x()*rr.y()*r.z()*nodes(i+1,j,k+1));
    }
    if (!yflat()) {
        u += (rr.x()*r.y()*rr.z()*nodes(i,j+1,k)
              + rr.x()*r.y()*r.z()*nodes(i,j+1,k+1));
    }
    if (!(xflat() || yflat())) {
        u += (r.x()*r.y()*rr.z()*nodes(i+1,j+1,k)
              + r.x()*r.y()*r.z()*nodes(i+1,j+1,k+1));
    }
    assert(std::isnormal(u) || u == 0.0);// need to move to c++11 fpclassify
    return u;
}

double
Mesh3d::at(Point3d const& pos, Index3d const& ul, Array3d<double> const& nodes) const
{
    Point3d r = ratios(pos, ul);
    return interpolate(r, ul, nodes);
}
