/*
 * bend_brent.cc - Brent's line minimization for bending solver
 *                 (taken from Numerical Recipe's brent.c)
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#include "bend3d.hpp"
#include <cmath>
#include <limits>

namespace {
    int const ITMAX = 100;
    double const CGOLD = 0.3819660;
    double const ZEPS = 1.0e-10;
    
    template <typename N> 
    void shft(N& a, N& b, N& c, N const& d) {
        a = b;
        b = c;
        c = d;
    }
    
    template <typename F, typename S>
    F sign(F a, S b) {
        return ((b >= 0) == (a >= 0)) ? a : -a;
    }
}

double BendingSolver3d::brent(double ax, double bx, double cx,
			      double *xmin, PF1DIM pfunc)
{
    double e = 0.0;
    double d = 0;
    double a = ax < cx ? ax : cx;
    double b = ax > cx ? ax : cx;
    double x = bx;
    double w = bx;
    double v = bx;
    double fw = (this->*pfunc)(bx);
    double fv = fw;
    double fx = fw;
    double u  = std::numeric_limits<double>::quiet_NaN();
    double fu = std::numeric_limits<double>::quiet_NaN();

    for (int iter=1;iter<=ITMAX;iter++) {
        double xm = 0.5*(a+b);
        double tol1 = brent_tol*fabs(x)+ZEPS;
        double tol2 = 2.0*tol1;
        if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
            *xmin=x;
            return fx;
        }
        if (fabs(e) > tol1) {
            double r = (x-w)*(fx-fv);
            double q=(x-v)*(fx-fw);
            double p=(x-v)*q-(x-w)*r;
            q = 2.0*(q-r);
            if (q > 0.0) {
                p = -p;
            } else {
                q = -q;
            }
            double etemp = e;
            e=d;
            if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) {
                e = x >= xm ? a-x : b-x;
                d = CGOLD*e;
            } else {
                d=p/q;
                u=x+d;
                if (u-a < tol2 || b-u < tol2)
                    d = sign(tol1,xm-x);
            }
        } else {
            d=CGOLD*(e=(x >= xm ? a-x : b-x));
        }
        u = fabs(d) >= tol1 ? x+d : x+sign(tol1,d);
        fu=(this->*pfunc)(u);
        if (fu <= fx) {
            if (u >= x)
                a=x; 
            else
                b=x;
            shft(v,w,x,u);
            shft(fv,fw,fx,fu);
        } else {
            if (u < x) a=u;
            else b=u;
            if (fu <= fw || w == x) {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            } else if (fu <= fv || v == x || v == w) {
                v=u;
                fv=fu;
            }
        }
    }
    error("BendingSolver3d::Too many iterations in brent");
    *xmin=x;
    return fx;
}
