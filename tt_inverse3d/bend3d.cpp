/*
 * bend3d.cc - ray-bending solver implementation
 *           direct minimization of travel times by conjugate gradient
 *           (rays are represented as beta-splines)
 * based on bend.cc
 *
 * Adria Melendez and Jun Korenaga
 * Spring 2011
 */

#include <iostream> // remove this 2
#include <fstream>
#include <cmath>
#include "bend3d.hpp"
#include "traveltime3d.hpp"

BendingSolver3d::BendingSolver3d(const SlownessMesh3d& s,
				 const BetaSpline3d& b,
				 double tol1, double tol2) // Modified
    : smesh(s), bs(b), nintp(bs.numIntp()),
      cg_tol(tol1), brent_tol(tol2),
      point_p(0), direc_p(0),
      start_i_p(0), end_i_p(0),
      interf_p(0)
{
    Q.resize(nintp);
    dQdu.resize(nintp);

    int max_np = int(3*pow(double(smesh.nb_nodes()),1.0/3.0));
    dTdV.reserve(max_np);
    new_dTdV.reserve(max_np);
    direc.reserve(max_np);

    pp.reserve(max_np+4);
    new_pp.reserve(max_np+4);
}

int BendingSolver3d::refine(Array1d<Point3d>& path,
			    double& orig_time, double& new_time, bool ani)
{
    int np = check_path_size(path);
    makeBSpoints(path,pp);
    new_point.resize(np); 
    makeBSpoints(new_point,new_pp);
    is_Ani = ani;
    
    dTdV.resize(np);
    new_dTdV.resize(np);
    direc.resize(np);
    
    orig_time=calcTravelTime(smesh, path, bs, pp, Q, ani);
    double old_time = orig_time;
    
    calc_dTdV(smesh, path, bs, dTdV, pp, Q, dQdu, ani);
    direc = -dTdV; // p0
    // conjugate gradient search
    int iter_max = np*10;
    for (int iter=1; iter<=iter_max; iter++){
        new_time=line_min(path, direc);
        
        if ((old_time-new_time)<=cg_tol){
            return iter;
        } else {
	  calc_dTdV(smesh, path, bs, new_dTdV, pp, Q, dQdu, ani);
            double gg=0.0, dgg=0.0;
            for (int i=1; i<=np; i++){
                double x0=dTdV(i).x(),y0=dTdV(i).y(),z0=dTdV(i).z();
                double x1=new_dTdV(i).x(),y1=new_dTdV(i).y(),z1=new_dTdV(i).z();
                gg += x0*x0+y0*y0+z0*z0;
                dgg += x1*x1+y1*y1+z1*z1;
            }
            if (abs(gg)<1e-10){
                return iter;
            }else{
                dgg /= gg;
            }
            
            for (int i=1; i<=np; i++){
                direc(i) = -new_dTdV(i)+dgg*direc(i);
            }
            
            dTdV = new_dTdV;
            old_time = new_time;
        }
    }
    return -iter_max; // too many iterations
}

int BendingSolver3d::refine(Array1d<Point3d>& path,
			    double& orig_time, double& new_time,
			    const std::vector<int>& start_i,
			    const std::vector<int>& end_i,
			    const Array1d<const Interface3d*>& interf, bool ani)
{
    int np=check_path_size(path,start_i,end_i,interf);
    makeBSpoints(path,pp);
    new_point.resize(np); 
    makeBSpoints(new_point,new_pp);
    is_Ani = ani;
    
    dTdV.resize(np);
    new_dTdV.resize(np);
    direc.resize(np);
    orig_time=calcTravelTime(smesh, path, bs, pp, Q, ani);
    double old_time = orig_time;

    calc_dTdV(smesh, path, bs, dTdV, pp, Q, dQdu, ani);
    adjust_dTdV(dTdV, path, start_i, end_i, interf);
    direc = -dTdV; // p0
    
    // conjugate gradient search
    int iter_max = np*10;
    for (int iter=1; iter<=iter_max; iter++){
        new_time=line_min(path, direc, start_i, end_i, interf);
        if ((old_time-new_time)<=cg_tol){
            return iter;
        }

        calc_dTdV(smesh, path, bs, new_dTdV, pp, Q, dQdu, ani);
        adjust_dTdV(new_dTdV, path, start_i, end_i, interf);

        double gg=0.0, dgg=0.0;
        for (int i=1; i<=np; i++){
            Point3d const& di    = dTdV(i);
            Point3d const& newdi = new_dTdV(i);
            double x0 = di.x();
            double y0 = di.y();
            double z0 = di.z();
            double x1 = newdi.x();
            double y1 = newdi.y();
            double z1 = newdi.z();

            gg  += x0*x0+y0*y0+z0*z0;
            dgg += x1*x1+y1*y1+z1*z1;
        }
        //	if (gg==0.0){ // 2000.5.8 
        // this may be the cause for "-NaN" reported by
        // Allegra. (due to division by extremely small gg)
        if (abs(gg)<1e-10){
            return iter;
        }else{
            dgg /= gg;
        }
        for (int i=1; i<=np; i++){
            direc(i) = -new_dTdV(i)+dgg*direc(i);
        }
        dTdV = new_dTdV;
        old_time = new_time;
    }
    return -iter_max; // too many iterations
}

void BendingSolver3d::adjust_dTdV(Array1d<Point3d>& dTdV,
				  const Array1d<Point3d>& path,
				  const std::vector<int>& start_i,
				  const std::vector<int>& end_i,
				  const Array1d<const Interface3d*>& interf) // Is it necessary for dTdV(i).z() to be zero? YES
{
    for (int a=1; a<=interf.size(); a++){
        boost::tuple<double,double> slopes = interf(a)->dzdxy(path(start_i[a-1]).x(),path(start_i[a-1]).y());
        double slope_x = slopes.get<0>();
        double slope_y = slopes.get<1>();
        double total_x = 0.0;
        double total_y = 0.0;

        for (int i=start_i[a-1]; i<=end_i[a-1]; i++){
            Point3d& di = dTdV(i);
            double x2 = std::pow(di.x(),2);
            double y2 = std::pow(di.y(),2);
            if((x2+y2) > 0) {
                double r = x2/(x2+y2);
                total_x += di.x() + slope_x*di.z()*r;
                total_y += di.y() + slope_y*di.z()*(1-r);
            }else{
                total_x += slope_x*di.z()/2;
                total_y += slope_y*di.z()/2;
            }
            di.z() = 0.0;
        }
        total_x /= (end_i[a-1]-start_i[a-1]+1);
        total_y /= (end_i[a-1]-start_i[a-1]+1);

        for (int i=start_i[a-1]; i<=end_i[a-1]; i++){
            dTdV(i).x() = total_x; // all points have the same direction now.
            dTdV(i).y() = total_y;
        }
    }
}

int BendingSolver3d::check_path_size(Array1d<Point3d>& path) const
{
    assert(path.size() >= 2);
    return path.size();
}

int BendingSolver3d::check_path_size(Array1d<Point3d>& path,
				     const std::vector<int>& start_i,
				     const std::vector<int>& end_i,
				     const Array1d<const Interface3d*>& interf) const
{
    assert(path.size() >= 2);
    int np=path.size();

    int nintf = interf.size();
    if (start_i.size() != nintf || end_i.size() != nintf)
        error("BendingSolver3d::refine - size mismatch in interface spec");
    for (int i=1; i<=nintf; i++){
        if (start_i[i-1] > end_i[i-1] ||
            start_i[i-1] < 1 || start_i[i-1] > np ||
            end_i[i-1] < 1 || end_i[i-1] > np)
            error("BendingSolver3d::refine - invalid interface points");
    }
    return np;
}

double BendingSolver3d::line_min(Array1d<Point3d>& path, const Array1d<Point3d>& direc) // No need to be modified
{
    point_p = &path; // for f1dim()
    direc_p = &direc;

    double ax = 0.0;
    double xx = 1.0;
    double bx, fa, fx, fb, xmin;
    mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, &BendingSolver3d::f1dim);
    double new_val = brent(ax,xx,bx,&xmin, &BendingSolver3d::f1dim);
    point_p = 0;
    direc_p = 0;
    
    for (int i=1; i<=path.size(); i++){
        path(i) += xmin*direc(i);
    }
    return new_val;
}

double BendingSolver3d::f1dim(double x) // No need to be modified
{
    for (int i=0; i < point_p->size(); ++i){
        new_point.stl()[i] = (*point_p)[i]+x*(*direc_p)[i];
    }
    double d = calcTravelTime(smesh, new_point, bs, new_pp, Q, is_Ani);
    return d;
}

double BendingSolver3d::line_min(Array1d<Point3d>& path, const Array1d<Point3d>& direc,
				 const std::vector<int>& start_i, const std::vector<int>& end_i,
				 const Array1d<const Interface3d*>& interf) // Modified
{
    point_p   = &path; // for f1dim_interf()
    direc_p   = &direc;
    start_i_p = &start_i;
    end_i_p   = &end_i;
    interf_p  = &interf;

    //    double ax = 0.0;
    //    double xx = 1.0;
    double ax = -0.1;
    double xx = 0.1;
    double bx, fa, fx, fb, xmin;
    mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, &BendingSolver3d::f1dim_interf);
    double new_val = brent(ax,xx,bx,&xmin, &BendingSolver3d::f1dim_interf);
    point_p = 0;
    direc_p = 0;
    interf_p = 0;
    start_i_p = 0;
    end_i_p   = 0;

    for (int i=1; i<=path.size(); i++){
        path(i) += xmin*direc(i);
    }
    for (int a=1; a<=interf.size(); a++){
        double new_z = interf(a)->z(path(start_i[a-1]).x(),path(start_i[a-1]).y());
        for (int i=start_i[a-1]; i<=end_i[a-1]; i++){
            path(i).z() = new_z;
        }
    }
    return new_val;
}

double BendingSolver3d::f1dim_interf(double x) // Modified
{    
    std::vector<Point3d>&  new_points = new_point.stl();
    
    for (int i=0; i < point_p->size(); ++i) {
        new_points[i] = (*point_p)[i] + x*(*direc_p)[i];
    }

    for (int a=0; a < interf_p->size(); ++a){
        int start = (*start_i_p)[a];
        int end   = (*end_i_p)[a];
        Point3d new_point = new_points[start-1];
        double new_z = (*interf_p)[a]->z(new_point.x(),new_point.y());
        for (int i = start; i < end; ++i){
            new_points[i-1].z() = new_z;
        }
    }
    double d = calcTravelTime(smesh, new_point, bs, new_pp, Q, is_Ani);
    return d;
}

