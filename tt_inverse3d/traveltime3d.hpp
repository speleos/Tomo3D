/*
 * traveltime3d.h - traveltime related helper functions
 * based on traveltime.h
 *
 * Adria Melendez and Jun Korenaga
 * Spring 2011
 */

#ifndef _TOMO_TRAVELTIME3D_H_
#define _TOMO_TRAVELTIME3D_H_

#include "array.hpp" // from mconv
#include "geom3d.hpp"
#include "smesh3d.hpp"
#include "betaspline3d.hpp"

double calcTravelTime(const SlownessMesh3d& m, const Array1d<Point3d>& path,
		      const BetaSpline3d& bs,
		      const Array1d<const Point3d*>& pp,
		      Array1d<Point3d>& Q, bool ani);

void calc_dTdV(const SlownessMesh3d& m, const Array1d<Point3d>& path,
	       const BetaSpline3d& bs, Array1d<Point3d>& dTdV,
	       const Array1d<const Point3d*>& pp,
	       Array1d<Point3d>& Q, Array1d<Point3d>& dQdu, bool ani);

#endif /* _TOMO_TRAVELTIME3D_H_ */
