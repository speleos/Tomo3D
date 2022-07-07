/*
 * corrlen3d.h based on corrlen.h
 * correlation length functions' interface
 *
 * Adria Melendez and Jun Korenaga
 * Fall 2011
 */

#ifndef _TOMO_CORRLEN3D_H_
#define _TOMO_CORRLEN3D_H_

#include <string>

#include "array.hpp"
#include "mesh.hpp"
#include "geom3d.hpp"
#include "index3d.hpp"

class CorrelationLength2d {
public:
    CorrelationLength2d(const char *fn);
    double atx(double, double) const;
    double aty(double, double) const;

private:
    std::vector<double> xpos, ypos;
    Array2d<double> xval, yval;
};

class CorrelationLength3d : public Mesh3d {
public:
    CorrelationLength3d(std::string const& fn);
    void at(const Point3d&, double&, double&, double&) const;

private:
    void upperleft(const Point3d& pos, Index3d& guess) const;
    void calc_local(const Point3d& pos, int, int, int,
		    double&, double&, double&, double&, double&, double&) const;

    Array3d<double> hxgrid, hygrid, vgrid;
};

class DampingWeight3d {
public:
    DampingWeight3d(const char *fn);
    void at(const Point3d&, double&) const;

private:
    void upperleft(const Point3d& pos, Index3d& guess) const;
    void calc_local(const Point3d& pos, int, int, int,
		    double&, double&, double&, double&, double&, double&) const;

    int nx, ny, nz;
    Array3d<double> wgrid;
    Array2d<double> topo, bx_vec, by_vec;
    std::vector<double> xpos, ypos, zpos;
    Array1d<double> rdx_vec, rdy_vec, rdz_vec;
};

#endif /* _TOMO_CORRLEN3D_H_ */
