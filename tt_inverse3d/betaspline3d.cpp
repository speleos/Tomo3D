/*
 * betaspline3d.cc - beta spline implementation
 * based on betaspline.cc
 *
 * Adria Melendez and Jun Korenaga
 * Spring 2011
 */

#include <cassert>

#include "betaspline3d.hpp"
#include "error.hpp" // from mconv

BetaSpline3d::BetaSpline3d(double beta1, double beta2, int n)
  : nintp(n)
{
  if (nintp<=2) error("BetaSpline3d::invalid nintp");
  calc_c(beta1,beta2);
  calc_u();
  calc_B();
  calc_dBdu();
}

void BetaSpline3d::interpolate( Point3d const& p0, Point3d const& p1,
                                Point3d const& p2, Point3d const& p3,
                                vector<Point3d>& Q) const
{
    assert(Q.size() == nintp);
    for (int i=0; i<nintp; ++i) {
        Q[i] = b[0][i]*p0+b[1][i]*p1+b[2][i]*p2+b[3][i]*p3;
    }
}

void BetaSpline3d::interpolate(const Point3d& p0, const Point3d& p1,
			       const Point3d& p2, const Point3d& p3,
			       vector<Point3d>& Q, vector<Point3d>& dQdu) const
{
    assert(Q.size() == nintp);
    interpolate(p0, p1, p2, p3, Q);
    for (int i=0; i<nintp; ++i){
        dQdu[i] = dbdu[0][i]*p0 +dbdu[1][i]*p1 +dbdu[2][i]*p2 +dbdu[3][i]*p3;
    }
}

void BetaSpline3d::resetNIntp(int n)
{
    if (n<=2) error("BetaSpline3d::invalid nintp");
    nintp = n;
    calc_u();
    calc_B();
    calc_dBdu();
}

void BetaSpline3d::calc_c(double beta1, double beta2)
{
    double beta1_2 = beta1*beta1;
    double beta1_3 = beta1_2*beta1;
    double delta = 2.0*beta1_3+4.0*(beta1_2+beta1)+beta2+2.0;
    double rdelta = 1.0/delta;
    
    c[0][0] = 2.0*beta1_3*rdelta;
    c[0][1] = (4.0*(beta1_2+beta1)+beta2)*rdelta;
    c[0][2] = 2.0*rdelta;
    c[0][3] = 0.0;
    
    c[1][0] = -6.0*beta1_3*rdelta;
    c[1][1] = 6.0*beta1*(beta1_2-1.0)*rdelta;
    c[1][2] = 6.0*beta1*rdelta;
    c[1][3] = 0.0;
    
    c[2][0] = -1.0*c[1][0];
    c[2][1] = 3.0*(-2.0*beta1_3-2.0*beta1_2-beta2)*rdelta;
    c[2][2] = 3.0*(2.0*beta1_2+beta2)*rdelta;
    c[2][3] = 0.0;
    
    c[3][0] = c[1][0]/3.0;
    c[3][1] = 2.0*(beta1_3+beta1_2+beta1+beta2)*rdelta;
    c[3][2] = -2.0*(beta1_2+beta1+beta2+1.0)*rdelta;
    c[3][3] = c[0][2];
}

void BetaSpline3d::calc_u()
{
    u.resize(nintp);
    
    // set local coordinate
    double du = 1.0/(nintp-1);
    for (int i=0; i<nintp; i++){
        u[i] = i*du;
    }
}

void BetaSpline3d::calc_B()
{
  for (int k=0; k<4; k++)
      b[k].resize(nintp,0.0);

  for (int i=0; i<nintp; i++){
      double ug=1.0;
      for (int j=0; j<4; j++){
          for (int k=0; k<4; k++){
              b[k][i] += c[j][k]*ug;
          }
          ug *= u[i];
      }
  }
}    

void BetaSpline3d::calc_dBdu()
{
    for (int k=0; k<4; k++)
        dbdu[k].resize(nintp,0.0);
    
    for (int i=0; i<nintp; i++){
        double ug=1.0;
        for (int j=1; j<=3; j++){
            for (int k=0; k<4; k++){
                dbdu[k][i] += c[j][k]*j*ug;
            }
            ug *= u[i];
        }
    }
} 

// helper functions
void makeBSpoints(const vector<Point3d>& orig, vector<const Point3d*>& pp)
{
    pp.clear();
    int np=orig.size();
    pp.reserve(orig.size()+4);
    pp.push_back(&orig.front());
    pp.push_back(&orig.front());
    
    for (vector<Point3d>::const_iterator it = orig.begin();
         it != orig.end(); ++it) {
        pp.push_back(&(*it));
    }
    pp.push_back(&orig.back());
    pp.push_back(&orig.back());
}

void printCurve(ostream& os,
                const vector<Point3d>& orig, const BetaSpline3d& bs)
{
    vector<const Point3d*> pp;
    makeBSpoints(orig,pp);
    
    int np=orig.size();
    int nintp=bs.numIntp();
    vector<Point3d> Q(nintp);
    for (int i=0; i<= np; ++i){
        bs.interpolate(*pp[i],*pp[i+1],*pp[i+2],*pp[i+3],Q);
        for (int j=0; j<nintp; j++){
            os << Q[j].x() << " " << Q[j].y() << " " << Q[j].z() << '\n';
        }
    }
}
