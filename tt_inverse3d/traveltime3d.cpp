 /*
 * traveltime3d.cc - travel time integration along a ray path
 * based on traveltime.cc
 *
 * Adria Melendez and Jun Korenaga
 * Spring 2011
 */

#include "traveltime3d.hpp"
#include "interface3d.hpp"
#include <iostream>
#include <fstream>

double calcTravelTime(const SlownessMesh3d& smesh,
		      const Array1d<Point3d>& path,
		      const BetaSpline3d& bs,
		      const Array1d<const Point3d*>& pp,
		      Array1d<Point3d>& Q, bool ani)
{
  int np = path.size();
  int nintp = bs.numIntp();
  Index3d guess_index = smesh.nodeIndex(smesh.nearest(*pp(1)));

  double ttime=0.0;
  for (int i=1; i<=np+1; i++){
    int j1=i;
    int j2=i+1;
    int j3=i+2;
    int j4=i+3;
    bs.interpolate(*pp(j1),*pp(j2),*pp(j3),*pp(j4),Q);

    if (smesh.inWater(Q(nintp/2))){
      // rectangular integration
      for (int j=2; j<=nintp; j++){
	Point3d midp = 0.5*(Q(j-1)+Q(j));
	double u0 = smesh.at(midp,guess_index);
	double dist = Q(j).distance(Q(j-1));

	if (ani) {
	  double phi = smesh.Phi;
	  double delta;
	  delta = smesh.anid->at(midp);
	  double epsilon;
	  if (smesh.useVh){
	      double vh = smesh.anie->at(midp);
	      epsilon = (vh*u0)-1;
	  }else{
	      epsilon = smesh.anie->at(midp);
	  }
	  Point3d h(0.0,0.0,1.0);
	  Point3d ray = Q(j)-Q(j-1);
	  double cos_theta = ray.inner_product(h)/ray.norm();
	  double theta = acos(cos_theta);
	  double gamma, au0;
	  if (theta > asin(1.0)){
	    theta = acos(-1.0) - theta;
	  }
	  if(theta >= phi){
	    gamma = theta - phi;
	  }else if(theta < phi){
	    gamma = phi - theta;
	  }
	  au0 = u0/(1+delta*sin(gamma)*sin(gamma)*cos(gamma)*cos(gamma)
		    +epsilon*sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma));
	  ttime += au0*dist;
	}else{
	  ttime += u0*dist;
	}
      }		
    }else{
      // trapezoidal integration along one segment
      double u0 = smesh.at(Q(1),guess_index);
      for (int j=2; j<=nintp; j++){
	double u1 = smesh.at(Q(j),guess_index);
	double dist= Q(j).distance(Q(j-1));

	if (ani) {
	  double phi = smesh.Phi;
	  double delta0, delta1;
	  delta0 = smesh.anid->at(Q(j-1));
	  delta1 = smesh.anid->at(Q(j));
	  double epsilon0, epsilon1;
	  if(smesh.useVh){
	      double vh0 = smesh.anie->at(Q(j-1));
	      double vh1 = smesh.anie->at(Q(j));
	      epsilon0 = (vh0*u0)-1;
	      epsilon1 = (vh1*u1)-1;
	  }else{
	      epsilon0 = smesh.anie->at(Q(j-1));
	      epsilon1 = smesh.anie->at(Q(j));
	  }
	  Point3d h(0.0,0.0,1.0);
	  Point3d ray = Q(j)-Q(j-1);
	  double cos_theta = ray.inner_product(h)/ray.norm();
	  double theta = acos(cos_theta);
	  double gamma, au;
	  if (theta > asin(1.0)){
	    theta = acos(-1.0) - theta;
	  }
	  if(theta >= phi){
	    gamma = theta - phi;
	  }else if(theta < phi){
	    gamma = phi - theta;
	  }
	  double au0 = u0/(1+delta0*sin(gamma)*sin(gamma)*cos(gamma)*cos(gamma)
			   +epsilon0*sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma));
	  double au1 = u1/(1+delta1*sin(gamma)*sin(gamma)*cos(gamma)*cos(gamma)
			   +epsilon1*sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma));
	  au = 0.5*(au0+au1);
	  ttime += au*dist;
	}else{
	  ttime += 0.5*(u0+u1)*dist;
	}
	u0 = u1;
      }
    }
  }
  return ttime;
}

void calc_dTdV(const SlownessMesh3d& smesh, const Array1d<Point3d>& path,
	       const BetaSpline3d& bs, Array1d<Point3d>& dTdV,
	       const Array1d<const Point3d*>& pp,
	       Array1d<Point3d>& Q, Array1d<Point3d>& dQdu, bool ani)
{

  int np = path.size();
  int nintp = bs.numIntp();
  double du = 1.0/(nintp-1);

  dTdV(1).set(0.0,0.0,0.0); // end points are fixed
  dTdV(np).set(0.0,0.0,0.0);

  Index3d guess_index = smesh.nodeIndex(smesh.nearest(*pp(1)));
  for (int i=2; i<np; i++){
    double dTdVx=0.0, dTdVy=0.0, dTdVz=0.0;

    // loop over four segments affected by the point Vi
    int iseg_start = i-1;
    int iseg_end = i+2; // min(i+2,np-1);
    for (int j=iseg_start; j<=iseg_end; j++){
      int j1=j;
      int j2=j+1;
      int j3=j+2;
      int j4=j+3;
      bs.interpolate(*pp(j1),*pp(j2),*pp(j3),*pp(j4),Q,dQdu);
      int i_j = i-j;

      // trapezoidal integration along one segment
      Point3d startp = 0.9*Q(1)+0.1*Q(2); // in order to avoid singularity
      Point3d du0;
      double u0, au0;
      Point3d dad0, dae0, dad1, dae1;

      Point3d da0, da1;
      if (ani) {
	double phi = smesh.Phi;
	double delta0;
	delta0 = smesh.anid->at(startp,guess_index,dad0);
	double epsilon0;
	if(smesh.useVh){
	    double vh0 = smesh.anie->at(startp,guess_index,dae0);
	    u0 = smesh.at(startp,guess_index,du0);
//	    double v0 = 1/u0;
	    epsilon0 = (vh0*u0)-1;
//	    dae0 = (v0*dae0-vh0*du0)/(v0*v0);
	    dae0 = u0*(dae0-vh0*du0*u0);
	}else{
	    u0 = smesh.at(startp,guess_index,du0);
	    epsilon0 = smesh.anie->at(startp,guess_index,dae0);
	}
	Point3d h(0.0,0.0,1.0);
	Point3d ray = Q(2)-startp;
	double cos_theta = ray.inner_product(h)/ray.norm();
	double theta = acos(cos_theta);
	double gamma;
	if (theta > asin(1.0)){
	  theta = acos(-1.0) - theta;
	}
	if(theta >= phi){
	  gamma = theta - phi;
	}else if(theta < phi){
	  gamma = phi - theta;
	}
	double ps0 = 1+delta0*sin(gamma)*sin(gamma)*cos(gamma)*cos(gamma)
	    +epsilon0*sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma);
	au0 = u0/ps0;
	da0 = (sin(gamma)*sin(gamma)*cos(gamma)*cos(gamma)*dad0
	       +sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma)*dae0);
	du0 = (du0/ps0)-(u0*da0)/(ps0*ps0);
	u0 = au0;
      }else{
	u0 = smesh.at(startp,guess_index,du0);
      }
      double dQnorm0 = dQdu(1).norm();
      double integx0, integy0, integz0;

      if (dQnorm0 > 0){
	double tmp1 = u0*bs.coeff_dbdu(i_j,1)/dQnorm0;
	double tmp2 = dQnorm0*bs.coeff_b(i_j,1);
	integx0 = tmp1*dQdu(1).x()+tmp2*du0.x();
	integy0 = tmp1*dQdu(1).y()+tmp2*du0.y();
	integz0 = tmp1*dQdu(1).z()+tmp2*du0.z();
      } else {
	integx0 = integy0 = integz0 = 0.0;
      }
      
      for (int m=2; m<=nintp; m++){
	Point3d du1;
	double u1, au1;
	double ps1 = -1.0;
	if (m==nintp){
	  Point3d endp = 0.1*Q(m-1)+0.9*Q(m); // in order to avoid singularity
	  if (ani) {
	    double phi = smesh.Phi;
	    double delta1;
	    delta1 = smesh.anid->at(endp,guess_index,dad1);
	    double epsilon1;
	    if(smesh.useVh){
		double vh1 = smesh.anie->at(endp,guess_index,dae1);
		double u1 = smesh.at(endp,guess_index,du1);
//		double v1 = 1/u1;
		epsilon1 = (vh1*u1)-1;
//		dae1 = (v1*dae1-vh1*du1)/(v1*v1);
		dae1 = u1*(dae1-vh1*du1*u1);
	    }else{
		u1 = smesh.at(endp,guess_index,du1);
		epsilon1 = smesh.anie->at(endp,guess_index,dae1);
	    }
	    Point3d h(0.0,0.0,1.0);
	    Point3d ray = endp-Q(m-1); // this ray is parallel to the previous one by definition: could not think of any alternative.
	    double cos_theta = ray.inner_product(h)/ray.norm();
	    double theta = acos(cos_theta);
	    double gamma;
	    if (theta > asin(1.0)){
	      theta = acos(-1.0) - theta;
	    }
	    if(theta >= phi){
	      gamma = theta - phi;
	    }else if(theta < phi){
	      gamma = phi - theta;
	    }
	    ps1 = 1+delta1*sin(gamma)*sin(gamma)*cos(gamma)*cos(gamma)
		   +epsilon1*sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma);
	    au1 = u1/ps1;
	    da1 = (sin(gamma)*sin(gamma)*cos(gamma)*cos(gamma)*dad1
		   +sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma)*dae1);
	    du1 = (du1/ps1)-(u1*da1)/(ps1*ps1);
	    u1 = au1;
	  }else{
	    u1 = smesh.at(endp,guess_index,du1);
	  }
	}else{
	  if (ani) {
	    double phi = smesh.Phi;
	    double delta1;
	    delta1 = smesh.anid->at(Q(m),guess_index,dad1);
	    double epsilon1;
	    if(smesh.useVh){
		double vh1 = smesh.anie->at(Q(m),guess_index,dae1);
		u1 = smesh.at(Q(m),guess_index,du1);
//		double v1 = 1/u1;
		epsilon1 = (vh1*u1)-1;
//		dae1 = (v1*dae1-vh1*du1)/(v1*v1);
		dae1 = u1*(dae1-vh1*du1*u1);
	    }else{
		u1 = smesh.at(Q(m),guess_index,du1);
		epsilon1 = smesh.anie->at(Q(m),guess_index,dae1);
	    }
	    Point3d h(0.0,0.0,1.0);
	    Point3d ray = Q(m+1)-Q(m);
	    double cos_theta = ray.inner_product(h)/ray.norm();
	    double theta = acos(cos_theta);
	    double gamma;
	    if (theta > asin(1.0)){
	      theta = acos(-1.0) - theta;
	    }
	    if(theta >= phi){
	      gamma = theta - phi;
	    }else if(theta < phi){
	      gamma = phi - theta;
	    }
	    ps1 = 1+delta1*sin(gamma)*sin(gamma)*cos(gamma)*cos(gamma)
		   +epsilon1*sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma);
	    au1 = u1/ps1;
	    da1 = (sin(gamma)*sin(gamma)*cos(gamma)*cos(gamma)*dad1
		   +sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma)*dae1);
	    du1 = (du1/ps1)-(u1*da1)/(ps1*ps1);
	    u1 = au1;
	  }else{
	    u1 = smesh.at(Q(m),guess_index,du1);
	  }
	}

	double dQnorm1 = dQdu(m).norm();
	double integx1, integy1, integz1;
	if (dQnorm1 > 0){
	  double tmp1 = u1*bs.coeff_dbdu(i_j,m)/dQnorm1;
	  double tmp2 = dQnorm1*bs.coeff_b(i_j,m);
	  integx1 = tmp1*dQdu(m).x()+tmp2*du1.x();
	  integy1 = tmp1*dQdu(m).y()+tmp2*du1.y();
	  integz1 = tmp1*dQdu(m).z()+tmp2*du1.z();
	} else {
	  integx1 = integy1 = integz1 = 0.0;
	}
          
	dTdVx += 0.5*(integx0+integx1)*du;
	dTdVy += 0.5*(integy0+integy1)*du;
	dTdVz += 0.5*(integz0+integz1)*du;
          
	integx0 = integx1;
	integy0 = integy1;
	integz0 = integz1;
      }
    }
    dTdV(i).set(dTdVx,dTdVy,dTdVz);
  }
}
