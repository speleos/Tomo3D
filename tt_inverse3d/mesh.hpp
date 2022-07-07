#ifndef _TOMO_MESH3D_H_
#define _TOMO_MESH3D_H_

#include "array.hpp" // This was already there.
#include "geom3d.hpp"
#include "index3d.hpp"

class Mesh3d {
public:
    Mesh3d();
    Mesh3d( Mesh3d const& m);
    Mesh3d( Array1d<double>& x, Array1d<double>& y, Array1d<double>& z,
            Array2d<double>& topo);
    
    double dx(int i) const;
    double dy(int j) const;
    double dz(int k) const;

    double bx(int i, int j) const;
    double by(int i, int j) const;
    
    std::vector<double> const& xpos() const { return my_xpos.stl(); }
    std::vector<double> const& ypos() const { return my_ypos.stl(); }
    std::vector<double> const& zpos() const { return my_zpos.stl(); }
    Array2d<double> const& topo() const { return my_topo; }
    
    int nx() const { return my_xpos.size(); }
    int ny() const { return my_ypos.size(); }
    int nz() const { return my_zpos.size(); }

    int nb_nodes() const { return nx()*ny()*nz(); }
    

    double xpos(int i) const { return my_xpos(i); }
    double ypos(int j) const { return my_ypos(j); }
    double zpos(int k) const { return my_zpos(k); }
    double topo(int i, int j) const { return my_topo(i,j); }

    double xratio(double x, int i) const;
    double yratio(double x, int j) const;
    double zratio(double x, int j) const;
    double xratio(Point3d const& p, int i) const { return xratio(p.x(),i); }
    double yratio(Point3d const& p, int j) const { return yratio(p.y(),j); }
    double zratio(Point3d const& p, int k) const { return zratio(p.z(),k); }

    double xmin() const;
    double xmax() const;
    double ymin() const;
    double ymax() const;
    double zmin() const;
    double zmax() const;

    bool xflat() const { return my_xpos.size() == 1; }
    bool yflat() const { return my_ypos.size() == 1; }
    bool flat() const { return xflat() || yflat(); }
    
public:
    
    double at(Point3d const& pos, Index3d const& ul, Array3d<double> const& nodes) const;
    void resize(int nx, int ny, int nz);
    
    bool read(std::istream& i, Array1d<double>& data);
    bool read_xpos(std::istream& i);
    bool read_ypos(std::istream& i);
    bool read_zpos(std::istream& i);
    bool read_topo(std::istream& i);
    bool read(std::istream& i, Array2d<double>& data);

    bool read_dimensions(std::istream& in);

    Point3d ratios(Point3d const& p, Index3d const& pos) const;
    double  interpolate(Point3d const& r, Index3d const& ul, Array3d<double> const& nodes) const;
    void    upper_left(Point3d const& p, Index3d& guess) const;
    Index3d upper_left(Point3d const& p) const;

protected:
    Array1d<double> my_xpos, my_ypos, my_zpos;
    Array2d<double> my_topo;
};

#endif /* _TOMO_MESH3D_H_ */
