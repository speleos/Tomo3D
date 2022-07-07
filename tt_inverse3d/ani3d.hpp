/*
 * ani3d.hpp anisotropy mesh interface
 * based on smesh3d.hpp
 *
 * Adria Melendez, Alain Miniussi and Jun Korenaga
 * Winter 2016
 */

#ifndef _TOMO_ANI3D_H_
#define _TOMO_ANI3D_H_

#include <list>
#include "boost/optional.hpp"
#include "array.hpp"
#include "mesh.hpp"
#include "geom3d.hpp"
#include "heap_deque.hpp"
#include "index3d.hpp"

class AnisotropyMesh3d : public Mesh3d {
 public:
    AnisotropyMesh3d(std::string const& fname);
    AnisotropyMesh3d(AnisotropyMesh3d const& other);

    void set(const Array1d<double>&);
    void get(Array1d<double>&) const;
    void vget(Array1d<double>&) const;
    
    const Index3d& nodeIndex(int i) const;
    int nodeIndex(int i, int j, int k) const;
    Point3d nodePos(int i) const;

    // anisotropy interpolation
    double at(const Point3d& pos) const;
    double at(const Point3d& pos, Index3d& guess) const;
    double at(const Point3d& pos, Index3d& guess, Point3d& da) const;
    double atWater() const { return a_water; }
									   
    // cell-oriented functions
    int numCells() const { return cell_index.size(); }
    void cellNodes(int cell, std::vector<int>& neighb) const;
    int  cellIndex(int i, int j, int k) const { return index2cell(i,j,k); }
    void cellNormKernel(int, Array2d<double>&) const;

    typedef boost::optional<int> opt_idx;    
    int locate_in_cell( const Point3d&, Index3d&,
                        opt_idx&, opt_idx&, opt_idx&, opt_idx&, opt_idx&, opt_idx&, opt_idx&, opt_idx&, 
                        Point3d& r) const;

    // miscellaneous
    int nearest(const Point3d& src) const;
    bool inWater(const Point3d& pos) const;
    bool inAir(const Point3d& pos) const;

    // output functions
    void outMesh(ostream&) const;
    void printElements(ostream&) const;
    void printAGrid(ostream&,bool) const;
    void printAGrid(ostream&,
                    double,double,double,double,double,double,double,double,double) const;
    void printMaskGrid(ostream&, const std::vector<int>&) const;
    void printMaskGrid(ostream&, const Array1d<double>&) const;

    // smoothing-related functions
    boost::optional<double> xexp(Point3d const& p, int i, int j, int k, double lx) const;
    boost::optional<double> yexp(Point3d const& p, int i, int j, int k, double ly) const;
    boost::optional<double> zexp(Point3d const& p, int nodeIndex, double lz) const;
    
private:
    void commonNormKernel();
    bool in_water(const Point3d& pos, const Index3d& guess) const;
    bool in_air(const Point3d& pos, const Index3d& guess) const;
    //    double almost_exact_ttime(double v1, double v2,
    //			      double dpath) const;

    int ncells;
    double a_water; /* anisotropy of water column */
    double a_air; /* anisotropy of air column */
    Array3d<double> agrid;
    Array3d<int> ser_index, index2cell;
    Array1d<Index3d> node_index;
    std::vector<int> cell_index;
    Array2d<double> Sm_H1, Sm_H2, Sm_V, T_common;
};

#endif
