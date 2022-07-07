/*
 * bend_mnbrak.cc - initial bracket search for bending solver
 *                  (taken from Numerical Recipe's mnbrak.cc)
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#include "bend3d.hpp"
#include <math.h>

namespace {
    double const GOLD   = 1.618034;
    double const GLIMIT = 100.0;
    double const TINY   = 1.0e-20;
    template <typename V, typename S>
    V sign(V const& v, S const& s) {
        return (v >= 0) == (s >= 0) ? v : -v;
    }
    
    inline
    void
    shift( double& a, double& b, double& c, double const& d) {
        a = b;
        b = c;
        c = d;
    }
}

void BendingSolver3d::mnbrak(double *ax, double *bx, double *cx,
			     double *fa, double *fb, double *fc,
			     PF1DIM pfunc)
{
    *fa=(this->*pfunc)(*ax);
    *fb=(this->*pfunc)(*bx);
    if (*fb > *fa) {
        std::swap(*ax, *bx);
        std::swap(*fb, *fa);
    }
    *cx=(*bx)+GOLD*(*bx-*ax);
    *fc=(this->*pfunc)(*cx);
    while (*fb > *fc) {
        double r=(*bx-*ax)*(*fb-*fc);
        double q=(*bx-*cx)*(*fb-*fa);
        double u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
            (2.0*sign(max(fabs(q-r),TINY),q-r));
        double ulim=(*bx)+GLIMIT*(*cx-*bx);
        double fu;
        if ((*bx-u)*(u-*cx) > 0.0) {
            double lfu = (this->*pfunc)(u);
            if (lfu < *fc) {
                *ax = *bx;
                *bx = u;
                *fa = *fb;
                *fb = lfu;
                return;
            } else if (lfu > *fb) {
                *cx = u;
                *fc = lfu;
                return;
            } else {
                u=(*cx)+GOLD*(*cx-*bx);
                fu=(this->*pfunc)(u);
            }
        } else if ((*cx-u)*(u-ulim) > 0.0) {
            fu=(this->*pfunc)(u);
            if (fu < *fc) {
                *bx = *cx;
                *cx = u;
                u = *cx+GOLD*(*cx-*bx);
                *fb = *fc;
                *fc = fu;
                fu = (this->*pfunc)(u);
            }
        } else if ((u-ulim)*(ulim-*cx) >= 0.0) {
            u=ulim;
            fu=(this->*pfunc)(u);
        } else {
            u=(*cx)+GOLD*(*cx-*bx);
            fu=(this->*pfunc)(u);
        }
        shift(*ax,*bx,*cx,u);
        shift(*fa,*fb,*fc,fu);
    }
}
