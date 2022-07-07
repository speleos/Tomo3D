/*
 * interface3d.h based on interface.h
 *
 * Adria Melendez and Jun Korenaga
 * Spring 2011
 */

#ifndef _TOMO_INTERFACE3D_H_
#define _TOMO_INTERFACE3D_H_

#include "boost/tuple/tuple.hpp"

#include "array.hpp"
#include "smesh3d.hpp"

//class SlownessMesh3d;

class Interface3d {
public:
    Interface3d(){}
//    Interface3d(SlownessMesh3d const&);
    Interface3d(const SlownessMesh3d&);
//    Interface3d(const SlownessMesh3d const&);
    Interface3d(std::string const& fn);

    const Index2d& nodeIndex(int i) const;
    int nodeIndex(int i, int j) const;
    double z(double x, double y) const;
    boost::tuple<double,double> dzdxy(double x, double y) const;

    typedef boost::optional<int> opt_idx;
    void locate_in_plane(double, double, opt_idx& i0j0, opt_idx& i1j0, opt_idx& i0i1, opt_idx& i1j1) const;
    double x(int) const;
    double y(int) const;
    double xmin() const;
    double xmax() const;
    double ymin() const;
    double ymax() const;
    int numNodes() const { return nnodes; }
    int numCells() const { return ncells; }
    void cellNodes(int, int&, int&, int&, int&) const;
    void cellNormKernel(int, double&) const;
    void set(const Array1d<double>&);
    void get(Array1d<double>&) const;

    friend ostream&
    operator<<(ostream&, const Interface3d&);
    void printMaskGridRefl(ostream&, const Array1d<double>&) const;
    
    bool xflat() const { return xpos.size() == 1; }
    bool yflat() const { return ypos.size() == 1; }
private:
    double dzdx(double x, double y) const;
    double dzdy(double x, double y) const;

private:
    int nxr, nyr, nnodes, ncells;
    void calc_slope();
    static const double eps;
    std::vector<double> xpos;
    std::vector<double> ypos;
    Array2d<double> zpos;
    Array2d<double> slope_x, slope_y;
    Array2d<int> ser_index_r, index2cell_r;
    Array1d<Index2d> node_index_r;
    std::vector<int> cell_index_r;
};

#endif /* _TOMO_INTERFACE3D_H_ */

