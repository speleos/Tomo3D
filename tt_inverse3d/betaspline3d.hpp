/*
 * betaspline3d.h - beta spline interface
 * based on betaspline.h
 *
 * Adria Melendez and Jun Korenaga
 * Spring 2011
 */

#ifndef _TOMO_BETASPLINE3D_H_
#define _TOMO_BETASPLINE3D_H_

#include <vector>
#include <ostream>
#include "geom3d.hpp"

class BetaSpline3d {
    // note: quite ineficient, generate a lot of useless temps
public:
    BetaSpline3d(double beta1, double beta2, int nintp);

    int numIntp() const { return nintp; }
    void resetNIntp(int);
    void interpolate( Point3d const& p0, Point3d const& p1,
                      Point3d const& p2, Point3d const& p3,
                      std::vector<Point3d>& Q) const;
    void interpolate(const Point3d&, const Point3d&,
		     const Point3d&, const Point3d&,
		     std::vector<Point3d>&, std::vector<Point3d>&) const;

    // valid i range = -2,-1,0,1
    double coeff_b(int i, int iintp) const { return b[i+2][iintp-1]; }
    double coeff_dbdu(int i, int iintp) const { return dbdu[i+2][iintp-1]; }
    
private:
    void calc_c(double, double);
    void calc_u();
    void calc_B();
    void calc_dBdu();
    
    int nintp;
    std::vector<double> u, b[4], dbdu[4];
    double c[4][4];
};

// helper functions
void makeBSpoints(const std::vector<Point3d>& orig, std::vector<const Point3d*>& pp);
void printCurve( std::ostream& os,
                 const std::vector<Point3d>& orig, const BetaSpline3d& bs);
    
#endif /* _TOMO_BETASPLINE3D_H_ */
