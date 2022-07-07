/*
 * bend3d.h - ray-bending solver interface
 * based on bend.h
 *
 * Adria Melendez and Jun Korenaga
 * Spring 2011
 */

#ifndef _TOMO_BEND3D_H_
#define _TOMO_BEND3D_H_

#include "array.hpp" // from mconv
#include "geom3d.hpp"
#include "smesh3d.hpp"
#include "betaspline3d.hpp"
#include "interface3d.hpp"

class BendingSolver3d {
public:
    BendingSolver3d(const SlownessMesh3d& s, const BetaSpline3d& bs,
		    double tol1=1e-4, double tol2=1e-7);
  int refine(Array1d<Point3d>& path, double& orig_time, double& new_time, bool ani);
    int refine(Array1d<Point3d>& path, double& orig_time, double& new_time,
	       const std::vector<int>& start_i, const std::vector<int>& end_i,
	       const Array1d<const Interface3d*>& interf, bool ani);
    
private:
    typedef double (BendingSolver3d::*PF1DIM)(double);

    int check_path_size(Array1d<Point3d>& path) const;
    int check_path_size(Array1d<Point3d>& path,
			const std::vector<int>& start_i,
			const std::vector<int>& end_i,
			const Array1d<const Interface3d*>& interf) const;
    void adjust_dTdV(Array1d<Point3d>& dTdV,
		     const Array1d<Point3d>& path,
		     const std::vector<int>& start_i,
		     const std::vector<int>& end_i,
		     const Array1d<const Interface3d*>& interf);
    double line_min(Array1d<Point3d>& path, const Array1d<Point3d>& direc);
    double line_min(Array1d<Point3d>& path, const Array1d<Point3d>& direc,
		    const std::vector<int>& start_i, const std::vector<int>& end_i,
		    const Array1d<const Interface3d*>& interf);
    void mnbrak(double *ax, double *bx, double *cx,
		double *fa, double *fb, double *fc, PF1DIM);
    double brent(double ax, double bx, double cx, double *xmin, PF1DIM);
    double f1dim(double x);
    double f1dim_interf(double x);
		    
    const SlownessMesh3d& smesh;
    const BetaSpline3d& bs;
    const int nintp;
    const double cg_tol, brent_tol;

  bool is_Ani;

    Array1d<Point3d> Q, dQdu;
    Array1d<const Point3d*> pp, new_pp;
    std::vector<Point3d> const *point_p, *direc_p;
    std::vector<int> const     *start_i_p, *end_i_p;
    std::vector<const Interface3d*> const *interf_p;
    Array1d<Point3d> new_point;
    Array1d<Point3d> dTdV, new_dTdV, direc;
};

#endif /* _TOMO_BEND3D_H_ */
