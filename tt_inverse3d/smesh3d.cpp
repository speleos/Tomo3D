/*
 * smesh3d.cpp - slowness mesh implementation
 * based on smesh.cc
 *
 * Adria Melendez and Jun Korenaga
 * Spring 2011
 */
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include "smesh3d.hpp"
#include "interface3d.hpp"
#include "d_nr.hpp"

namespace details {
    int
    find_interval(double x, std::vector<double> const& coords) {
        assert(coords.size() != 0);
        typedef std::vector<double>::const_iterator iter;
        iter b = coords.begin();
        iter e = coords.end();
        // All position below the first interval are considered in 
        // the first interval. As if the first position was -inf
        ++b;
        // All position beyond the last interval are considered in 
        // the last interval. As if the last position was +inf
        --e;
        int pos;
        // From the preceeding comments, if there is only one interval,
        // it's considered covering -inf,+inf
        if (coords.size()==1) {
            pos = 0;
        } else {
            iter lb = std::upper_bound(b, e, x);
            pos = lb - b;
        }
        return pos;        
    }
}

SlownessMesh3d::SlownessMesh3d(SlownessMesh3d const&  o)
:   Mesh3d(o),
    p_water(o.p_water),
    p_air(o.p_air), 
    pgrid(o.pgrid),
    vgrid(o.vgrid),
    ser_index(o.ser_index),
    index2cell(o.index2cell),
    node_index(o.node_index),
    cell_index(o.cell_index),
    Sm_H1(o.Sm_H1),
    Sm_H2(o.Sm_H2),
    Sm_V(o.Sm_V),
    T_common(o.T_common),
    Phi(o.Phi),
    anid(o.anid),
    anie(o.anie),
    useVh(o.useVh) {}

SlownessMesh3d::SlownessMesh3d(std::string const& fname) // Modified
:   Mesh3d()
{
    ifstream s(fname.c_str());
    if (!s){
        cerr << "SlownessMesh3d::cannot open " << fname << "\n";
        exit(1);
    }

    int szx, szy, szz;
    s >> szx >> szy >> szz >> p_water >> p_air;
    Mesh3d::resize(szx,szy,szz);
    Mesh3d::read_xpos(s);
    Mesh3d::read_ypos(s);
    Mesh3d::read_topo(s);
    Mesh3d::read_zpos(s);
    
    if (p_water <= 0.0)
        cerr << "SlownessMesh3d::invalid water velocity\n";
    if (p_air <= 0.0)
        cerr << "SlownessMesh3d::invalid air velocity\n";
    p_water = 1.0/p_water;
    p_air = 1.0/p_air;
    
    pgrid.resize(nx(),ny(),nz());
    vgrid.resize(nx(),ny(),nz());
    ser_index.resize(nx(),ny(),nz());
    node_index.resize(nb_nodes());
    //even whit flat geometry, we need one inderection fake cell
    cell_index.resize(std::max(nx()-1,1)*std::max(ny()-1,1)*(nz()-1));
    index2cell.resize(std::max(nx()-1,1),std::max(ny()-1,1),nz()-1);

    int N=1;
    vgrid = 0.0;
    for (int i=1; i<=nx(); i++){
        for(int j=1; j<=ny(); j++){
            for (int k=1; k<=nz(); k++){
                double vel;
                s >> vel;
                vgrid(i,j,k) = vel;
                pgrid(i,j,k) = 1.0/vel;
                ser_index(i,j,k) = N;
                node_index(N).set(i,j,k);
                N++;
            }
        }
    }
    int icell=1;
    for (int i=1; i<std::max(nx(),2); i++){
        for(int j=1; j<std::max(ny(),2); j++){
            for (int k=1; k<nz(); k++){
                index2cell(i,j,k) = icell;
                cell_index[icell++-1] = ser_index(i,j,k);
            }
        }
    }

    // sanity check
    for (int i=2; i<=nx(); i++){
        if (xpos(i)<=xpos(i-1)){
            cerr << "SlownessMesh3d: illegal ordering of x nodes at xpos(" << i << ")="  << xpos(i)
                 << "<=" << "xpos(" << i-1 << ")=" << xpos(i-1) << '\n';
            exit(1);
        }
    }
    for (int j=2; j<=ny(); j++){
        if (ypos(j)<=ypos(j-1)){
            cerr << "SlownessMesh3d: illegal ordering of y nodes at ypos(" << j << ")="  << ypos(j)
                 << "<=" << "ypos(" << j-1 << ")=" << ypos(j-1) << '\n';
            exit(1);
        }
    }
    for (int i=2; i<=nz(); i++){
        if (zpos(i)<=zpos(i-1)){
            cerr << "SlownessMesh3d: illegal ordering of z nodes zpos(" << i << ")="  << zpos(i)
                 << "<=" << "zpos(" << i-1 << ")=" << zpos(i-1) << '\n';
            exit(1);
        }
    }
    for (int i=1; i<=nx(); i++){
        for(int j=1; j<=ny(); j++){
            for (int k=1; k<=nz(); k++){
                if (vgrid(i,j,k)<=0){
                    cerr << "SlownessMesh3d: non-positive velocity v("
                         << i << "," << j << "," << k << ")=" << vgrid(i,j,k) << '\n';
                    exit(1);
                }
            }
        }
    }
    
    commonNormKernel();
}

boost::optional<double>
SlownessMesh3d::xexp(Point3d const& p, int i, int j, int k, double lx) const {
    if (xflat()) {
        return boost::optional<double>(1.0);
    } else {
        int bnode = nodeIndex(i,j,k);
        double dx = nodePos(bnode).x()-p.x();
        if (abs(dx) <= lx){
            return boost::optional<double>(std::exp(-(dx*dx)/(lx*lx)));
        } else {
            return boost::optional<double>();
        }
    }
}

boost::optional<double>
SlownessMesh3d::yexp(Point3d const& p, int i, int j, int k, double ly) const {
    if (yflat()) {
        return boost::optional<double>(1.0);
    } else {
        int bnode = nodeIndex(i,j,k);
        double dy = nodePos(bnode).y()-p.y();
        if (abs(dy) <= ly){
            return boost::optional<double>(std::exp(-(dy*dy)/(ly*ly)));
        } else {
            return boost::optional<double>();
        }
    }
}

boost::optional<double>
SlownessMesh3d::zexp(Point3d const& p, int nodeIndex, double lz) const {
    double dz = nodePos(nodeIndex).z()-p.z();
    if (abs(dz) <= lz){
        return boost::optional<double>(std::exp(-(dz*dz)/(lz*lz)));
    } else {
        return boost::optional<double>();
    }
}
        
void SlownessMesh3d::set(const Array1d<double>& u) // Modified
{
    assert(u.size() == nb_nodes());

    int N=1;
    for (int i=1; i<=nx(); i++){
        for (int j=1; j<=ny(); j++){
            for (int k=1; k<=nz(); k++){
	      if(pgrid(i,j,k) <= 0.0){
		cerr << "SlownessMesh3d::set: Zero or negative velocity encountered:\n";
		cerr << "Node (i,j,k)=(" << i << "," << j << "," << k << ")\n";
		cerr << "Terminating.\n";
		exit(1);
	      }else{
                pgrid(i,j,k) = u(N);
                vgrid(i,j,k) = 1.0/pgrid(i,j,k);
                N++;
	      }
            }
        }
    }
}

void SlownessMesh3d::get(Array1d<double>& u) const // Modified
{
    std::vector<double> velocity_model;
    u.reserve(nb_nodes());
    for (int i=1; i<=nx(); ++i){
        for (int j=1; j<=ny(); ++j){
            for (int k=1; k<=nz(); ++k){
	      if(pgrid(i,j,k) <= 0.0){
		cerr << "SlownessMesh3d::get: Zero or negative velocity encountered:\n";
		cerr << "Node (i,j,k)=(" << i << "," << j << "," << k << ")\n";
		cerr << "Terminating.\n";
		exit(1);
	      }else{
                velocity_model.push_back(pgrid(i,j,k));
	      }
            }
        }
    }
    u.stl().swap(velocity_model);
}

void SlownessMesh3d::vget(Array1d<double>& v) const // Modified
{
    assert(v.size() == nb_nodes());

    int N=1;
    for (int i=1; i<=nx(); i++){
        for (int j=1; j<=ny(); j++){
            for (int k=1; k<=nz(); k++){
	      if(vgrid(i,j,k) <= 0.0){
		cerr << "SlownessMesh3d::vget: Zero or negative velocity encountered:\n";
		cerr << "Node (i,j,k)=(" << i << "," << j << "," << k << ")\n";
		cerr << "Terminating.\n";
		exit(1);
	      }else{
                v(N) = vgrid(i,j,k);
                N++;
	      }
            }
        }
    }
}

const Index3d&
SlownessMesh3d::nodeIndex(int i) const { return node_index(i); }

int SlownessMesh3d::nodeIndex(int i, int j, int k) const { return ser_index(i,j,k); } // Modified

Point3d SlownessMesh3d::nodePos(int i) const // Modified
{
    int node_i = node_index(i).i();
    int node_j = node_index(i).j();
    int node_k = node_index(i).k();
    Point3d p(xpos(node_i),ypos(node_j),zpos(node_k)+topo(node_i,node_j));

    return p;
}

void SlownessMesh3d::tiltAngle(double x)
{
  cerr << "x = " << x << '\n';
  //  cerr << "Phi = " << Phi << '\n';
  if(x>=0 && x<acos(-1.0)){
    Phi = x;
    //    cerr << "Phi = " << Phi << '\n';
  }else{
    cerr << "TomographicInversion3d::tiltAngle - tilt must be 0 or positive and smaller than Pi. Check tilt definition in the user guide.\n";
  }
  //  cerr << "Phi = " << Phi << '\n';
  //  cerr << "this = " << this->Phi << '\n';
}

void SlownessMesh3d::add_anisotropy(boost::shared_ptr<AnisotropyMesh3d> const& anidp,
				    boost::shared_ptr<AnisotropyMesh3d> const& aniep)
{
    anid = anidp;
    anie = aniep;
}

void SlownessMesh3d::usingVh()
{
    useVh = true;
}

double SlownessMesh3d::calc_ttime(int nsrc, int nrcv, bool ani) const // Modified
{

    if (nsrc<1 || nsrc>nb_nodes() || nrcv<1 || nrcv>nb_nodes())    error("SlownessMesh3d::calc_ttime - invalid input");

    Index3d src=node_index(nsrc);
    Index3d rcv=node_index(nrcv);
    int isrc=src.i(), jsrc=src.j(), ksrc=src.k();
    int ircv=rcv.i(), jrcv=rcv.j(), krcv=rcv.k();
    double xs=xpos(isrc), xr=xpos(ircv), ys=ypos(jsrc), yr=ypos(jrcv);
    double zs=topo(isrc,jsrc)+zpos(ksrc), zr=topo(ircv,jrcv)+zpos(krcv);
    Point3d sloc(xs,ys,zs);
    Point3d rloc(xr,yr,zr);
    int di = ircv-isrc;
    int dj = jrcv-jsrc;
    int dk = krcv-ksrc;

    /*
     *The following if-condition selects the plane in which we will interpolate velocity.
     *In this version 1.9, we maximize the number of interpolations.
     *We interpolate velocity at the planes perpendicular to the direction with max(abs(di),abs(dj)).
     *We may prefer to minimize interpolations by selecting min(abs(di),abs(dj)).
     *
     *We discard the z-direction because it is not possible to proceed as for both horizontal
     *directions (using proportionality) due to relief.
     */

    double ttime=0.0;
    //    cerr << "ani = " << ani << '\n';
    //    cerr << "Phi = " << Phi << '\n';
   
    if(abs(di)>abs(dj)){
        
        if(di<0){ //exchange src and rcv to get the first for-loop to work properly.
            swap(xs,xr);
            swap(ys,yr);
            swap(zs,zr);
            swap(isrc,ircv);
            swap(jsrc,jrcv);
            swap(ksrc,krcv);
            swap(sloc,rloc);
            di = -di;
            dj = -dj;
            dk = -dk;
        }

        double distx = xr-xs;
        double disty = yr-ys;
        double distz = zr-zs;
        double ratioyx = disty/distx; //Calculate proportions.
        double ratiozx = distz/distx;
        Point3d precut = sloc;
        double pave1 = pgrid(isrc,jsrc,ksrc);

        for(int i=isrc+1; i<=ircv; i++){
            double incx = (xpos(i)-xs);
            double ycut = ys+ratioyx*incx; //Find intersections.
            double zcut = zs+ratiozx*incx;
            Point3d cut(xpos(i),ycut,zcut);
            double lng = cut.distance(precut); //Distance between previous and present intersection points.

            int jstore = ( yflat()
                           ? 1 
                           : std::min( details::find_interval(ycut, my_ypos.stl())+1, ny()-1) );

            //Interpolate bathymetry at ycut to remove it from zcut and be able to locate it between vertical nodes.
	    double znb = zcut;
	    if(yflat()) {
	      znb -= topo(i,jstore);
	    } else {
	      double ratiob = (ycut-ypos(jstore))/(ypos(jstore+1)-ypos(jstore));
	      double bat = (1-ratiob)*topo(i,jstore)+ratiob*topo(i,jstore+1);
	      znb -= bat;
	    }
    
            //Velocity interpolation

            double pave2;
            if(znb<zpos(1)){
	      pave2 = zcut>=0 ? p_water : p_air;
            }else if(znb>zpos(nz())){
                double wy = (ycut-ypos(jstore))/(ypos(jstore+1)-ypos(jstore)); 
                pave2 = wy*pgrid(i,jstore+1,nz())+(1-wy)*pgrid(i,jstore,nz());
            }else{
                int kstore = std::min( details::find_interval(znb, my_zpos.stl())+1, nz()-1);

                double wy = yratio(ycut,jstore); 
                double wz = zratio(znb,kstore);
		// if(wz < 0){
		//   cerr << "wy=" << wy << ", wz=" << wz <<'\n';
		// }
                pave2 = (1-wy-wz+wy*wz)*pgrid(i,jstore,kstore)+(wz-wy*wz)*pgrid(i,jstore,kstore+1);
		if(!yflat()) {
		  pave2 += (wy-wy*wz)*pgrid(i,jstore+1,kstore)+wy*wz*pgrid(i,jstore+1,kstore+1) ;
		}
            }
	    if(lng<0){
	      cerr << "pave1=" << pave1 << ", pave2=" << pave2 << '\n';
	      cerr << "lng=" << lng << '\n';
	    }

	    if (ani) {
	      double delta1, delta2;
	      delta1 = anid->at(precut);
	      delta2 = anid->at(cut);
	      double epsilon1, epsilon2;
	      if(useVh){
		  double vh1 = anie->at(precut);
		  double vh2 = anie->at(cut);
		  epsilon1 = (vh1*pave1)-1;
		  epsilon2 = (vh2*pave2)-1;
	      }else{
		  epsilon1 = anie->at(precut);
		  epsilon2 = anie->at(cut);
	      }
	      Point3d h(0.0,0.0,1.0);
	      Point3d ray = cut-precut;
	      double cos_theta = ray.inner_product(h)/ray.norm();
	      double theta = acos(cos_theta);
	      double apave1, apave2;
	      double gamma, apave;
	      if (theta > asin(1.0)){
		theta = acos(-1.0) - theta;
	      }
	      if(theta >= Phi){
		gamma = theta - Phi;
	      }else if(theta < Phi){
		gamma = Phi - theta;
	      }
	      apave1 = pave1/(1+delta1*sin(gamma)*sin(gamma)*cos(gamma)*cos(gamma)
			    +epsilon1*sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma));
	      apave2 = pave2/(1+delta2*sin(gamma)*sin(gamma)*cos(gamma)*cos(gamma)
			    +epsilon2*sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma));	      
	      apave = 0.5*(apave1+apave2);
	      ttime += lng*apave;
	    }else{
	      double pave = (pave1+pave2)*0.5;
	      ttime += lng*pave;
	    }
            pave1=pave2;
            precut=cut;
        }

    }else if(abs(di)<abs(dj)){

        if(dj<0){ //exchange src and rcv to get the first for loop to work properly.
            swap(xs,xr);
            swap(ys,yr);
            swap(zs,zr);
            swap(isrc,ircv);
            swap(jsrc,jrcv);
            swap(ksrc,krcv);
            swap(sloc,rloc);
            di=-di;
            dj=-dj;
            dk=-dk;
        }

        double distx = xr-xs;
        double disty = yr-ys;
        double distz = zr-zs;
        double ratioxy = distx/disty;
        double ratiozy = distz/disty;
        Point3d precut = sloc;
        double pave1 = pgrid(isrc,jsrc,ksrc);

        for(int j=jsrc+1; j<=jrcv; j++){
            double incy = (ypos(j)-ys);
            double xcut = xs+ratioxy*incy;
            double zcut = zs+ratiozy*incy;
            Point3d cut(xcut,ypos(j),zcut);
            double lng = cut.distance(precut);
     
            int istore = ( xflat()
                           ? 1 
                           : std::min( details::find_interval(xcut, my_xpos.stl())+1, nx()-1) );
            //Interpolate bathymetry at xcut to remove it from zcut and be able to locate it between vertical nodes.

            double znb = zcut;
            if (xflat()) {
                znb -= topo(istore,j);
            } else {
                double ratiob = (xcut-xpos(istore))/(xpos(istore+1)-xpos(istore));
                double bat = (1-ratiob)*topo(istore,j)+ratiob*topo(istore+1,j);
                znb -= bat;
            }
            //Velocity interpolation

            double pave2;
            if(znb<zpos(1)){
                pave2 = zcut>=0 ? p_water : p_air;
            }else if(znb>zpos(nz())){
                double wx = (xcut-xpos(istore))/(xpos(istore+1)-xpos(istore)); 
                pave2 = wx*pgrid(istore+1,j,nz())+(1-wx)*pgrid(istore,j,nz());
            }else{
                int kstore = std::min( details::find_interval(znb, my_zpos.stl())+1, nz()-1);
                double wz = zratio(znb, kstore);
                double wx = xratio(xcut, istore);
		// if(wz <0){
		//   cerr << "wx=" << wx << ", wz=" << wz <<'\n';
		// }
                pave2 = ( (1-wx-wz+wx*wz)*pgrid(istore,j,kstore)
                          +(wz-wx*wz)*pgrid(istore,j,kstore+1) );
                if (!xflat()) {
                    pave2 += ((wx-wx*wz)*pgrid(istore+1,j,kstore)
                              +wx*wz*pgrid(istore+1,j,kstore+1));
                }
            }
	    if(lng<0){
	      cerr << "pave1=" << pave1 << ", pave2=" << pave2 << '\n';
	      cerr << "lng=" << lng << '\n';
	    }

	    if (ani) {
	      double delta1, delta2;
	      delta1 = anid->at(precut);
	      delta2 = anid->at(cut);
	      double epsilon1, epsilon2;
	      if(useVh){
		  double vh1 = anie->at(precut);
		  double vh2 = anie->at(cut);
		  epsilon1 = (vh1*pave1)-1;
		  epsilon2 = (vh2*pave2)-1;
	      }else{
		  epsilon1 = anie->at(precut);
		  epsilon2 = anie->at(cut);
	      }
	      Point3d h(0.0,0.0,1.0);
	      Point3d ray = cut-precut;
	      double cos_theta = ray.inner_product(h)/ray.norm();
	      double theta = acos(cos_theta);
	      double apave1, apave2;
	      double gamma, apave;
	      if (theta > asin(1.0)){
		theta = acos(-1.0) - theta;
	      }
	      if(theta >= Phi){
		gamma = theta - Phi;
	      }else if(theta < Phi){
		gamma = Phi - theta;
	      }
	      apave1 = pave1/(1+delta1*sin(gamma)*sin(gamma)*cos(gamma)*cos(gamma)
			    +epsilon1*sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma));
	      apave2 = pave2/(1+delta2*sin(gamma)*sin(gamma)*cos(gamma)*cos(gamma)
			    +epsilon2*sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma));	      
	      apave = 0.5*(apave1+apave2);
	      ttime += lng*apave;
	    }else{
	      double pave = (pave1+pave2)*0.5;
	      ttime += lng*pave;
	    }
            pave1=pave2;
            precut=cut;
        }  
    
    }else if(abs(di)==abs(dj)){

        if(di==0){ //vertical path particular case
            if(dk<0){ //exchange src and rcv for the first loop to work.
                swap(isrc,ircv);
                swap(jsrc,jrcv);
                swap(ksrc,krcv);
            }
      
            double zprev = zpos(ksrc);
            double pave1 = pgrid(isrc,jsrc,ksrc);
            for(int k=ksrc+1; k<=krcv; k++){
                double lng = fabs(zpos(k)-zprev); //Distance between previous and present nodes.

		if (ani) {
		  Point3d Zprev(0.0,0.0,zprev);
		  Point3d Zpos(0.0,0.0,zpos(k));
		  double delta1, delta2;
		  delta1 = anid->at(Zprev);
		  delta2 = anid->at(Zpos);
		  double epsilon1, epsilon2;
		  if(useVh){
		      double vh1 = anie->at(Zprev);
		      double vh2 = anie->at(Zpos);
		      epsilon1 = (vh1*pave1)-1;
		      epsilon2 = (vh2*pgrid(isrc,jsrc,k))-1;
		  }else{
		      epsilon1 = anie->at(Zprev);
		      epsilon2 = anie->at(Zpos);
		  }
		  double apave1, apave2;
		  double apave;
    
		  double gamma = Phi;
		  apave1 = pave1/(1+delta1*sin(gamma)*sin(gamma)*cos(gamma)*cos(gamma)
				  +epsilon1*sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma));
		  apave2 = pgrid(isrc,jsrc,k)/(1+delta2*sin(gamma)*sin(gamma)*cos(gamma)*cos(gamma)
					       +epsilon2*sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma));	      
		  apave = 0.5*(apave1+apave2);
		  ttime += lng*apave;
		}else{
		  double pave = (pave1+pgrid(isrc,jsrc,k))*0.5;
		  ttime += lng*pave;
		}
                zprev=zpos(k);
                pave1=pgrid(isrc,jsrc,k);
            }
        }else{
            if(di<0){ //exchange src and rcv to get the first for loop to work properly.
                swap(xs,xr);
                swap(ys,yr);
                swap(zs,zr);
                swap(isrc,ircv);
                swap(jsrc,jrcv);
                swap(ksrc,krcv);
                swap(sloc,rloc);
                di=-di;
                dj=-dj;
                dk=-dk;
            }
      
            double distx = xr-xs;
            double disty = yr-ys;
            double distz = zr-zs;
            double ratioyx = disty/distx; //Calculate proportions.
            double ratiozx = distz/distx;
            Point3d precut = sloc;
            double pave1 = pgrid(isrc,jsrc,ksrc);
      
            for(int i=isrc+1; i<=ircv; i++){
                double incx = (xpos(i)-xs);
                double ycut = ys+ratioyx*incx; //Find intersections.
                double zcut = zs+ratiozx*incx;
                Point3d cut(xpos(i),ycut,zcut);
                double lng = cut.distance(precut); //Distance between previous and present intersection points.
    
                int jstore = ( yflat()
                               ? 1 
                               : std::min( details::find_interval(ycut, my_ypos.stl())+1, ny()-1) );

                //Interpolate bathymetry at ycut to remove it from zcut and be able to locate it between vertical nodes.
                double znb = zcut;
                if (yflat()) {
                    znb -= topo(i,jstore);
                } else {
                    double ratiob = (ycut-ypos(jstore))/(ypos(jstore+1)-ypos(jstore));
                    double bat = (1-ratiob)*topo(i,jstore)+ratiob*topo(i,jstore+1);
                    znb -= bat;
                }
                //Velocity interpolation

                double pave2;
                if(znb<zpos(1)){
                    pave2 = zcut >= 0 ? p_water : p_air;
                }else if(znb>zpos(nz())){
                    double wy = (ycut-ypos(jstore))/(ypos(jstore+1)-ypos(jstore)); 
                    pave2 = wy*pgrid(i,jstore+1,nz())+(1-wy)*pgrid(i,jstore,nz());
                }else{
                    int kstore = std::min( details::find_interval(znb, my_zpos.stl())+1, nz()-1);
                    double wy = yratio(ycut,jstore); 
                    double wz = zratio(znb,kstore);
		    // if (wz<0){
		    //   cerr << "wy=" << wy << ", wz=" << wz <<'\n';
		    // }
                    pave2 = (1-wy-wz+wy*wz)*pgrid(i,jstore,kstore)+(wz-wy*wz)*pgrid(i,jstore,kstore+1);
		    if (!yflat()) {
		      pave2 += (wy-wy*wz)*pgrid(i,jstore+1,kstore)+wy*wz*pgrid(i,jstore+1,kstore+1);
		    }
                }
		if(lng<0){
		  cerr << "pave1=" << pave1 << ", pave2=" << pave2 << '\n';
		  cerr << "lng=" << lng << '\n';
		}

		if (ani) {
		  double delta1, delta2;
		  delta1 = anid->at(precut);
		  delta2 = anid->at(cut);
		  double epsilon1, epsilon2;
		  if(useVh){
		      double vh1 = anie->at(precut);
		      double vh2 = anie->at(cut);
		      epsilon1 = (vh1*pave1)-1;
		      epsilon2 = (vh2*pave2)-1;
		  }else{
		      epsilon1 = anie->at(precut);
		      epsilon2 = anie->at(cut);
		  }
		  Point3d h(0.0,0.0,1.0);
		  Point3d ray = cut-precut;
		  double cos_theta = ray.inner_product(h)/ray.norm();
		  double theta = acos(cos_theta);
		  double apave1, apave2;
		  double gamma, apave;
		  if (theta > asin(1.0)){
		    theta = acos(-1.0) - theta;
		  }
		  if(theta >= Phi){
		    gamma = theta - Phi;
		  }else if(theta < Phi){
		    gamma = Phi - theta;
		  }
		  apave1 = pave1/(1+delta1*sin(gamma)*sin(gamma)*cos(gamma)*cos(gamma)
				  +epsilon1*sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma));
		  apave2 = pave2/(1+delta2*sin(gamma)*sin(gamma)*cos(gamma)*cos(gamma)
				  +epsilon2*sin(gamma)*sin(gamma)*sin(gamma)*sin(gamma));	      
		  apave = 0.5*(apave1+apave2);
		  ttime += lng*apave;
		}else{
		  double pave = (pave1+pave2)*0.5;
		  ttime += lng*pave;
		}
                pave1=pave2;
                precut=cut;
		
            }
        }
    }

    if (ttime<0){
        cerr << "SlownessMesh3d::calc_ttime - negative traveltime encountered for ("
             << nsrc << ", " << nrcv << ")\n";
	cerr << "ttime=" << ttime << '\n';
        cerr << "pgrid(" << isrc << "," << jsrc << "," << ksrc << ")=" << pgrid(isrc,jsrc,ksrc)
             << ", pgrid(" << ircv << "," << jrcv << "," << krcv << ")=" << pgrid(ircv,jrcv,krcv)
             << '\n';
	cerr << "xsrc=" << xs << ", ysrc=" << ys << ", zsrc=" << zs << ", zsrc in grid=" << zs-topo(isrc,jsrc) << "=" << zpos(ksrc) << '\n';
	cerr << "xrcv=" << xr << ", yrcv=" << yr << ", zrcv=" << zr << ", zrcv in grid=" << zr-topo(ircv,jrcv) << "=" << zpos(krcv) << '\n';
        exit(1);
    }

    return ttime;
}

double SlownessMesh3d::at(const Point3d& pos, Index3d& guess) const
{
    upper_left(pos,guess);
    if (in_water(pos,guess)) {
        return p_water;
    }
    if (in_air(pos,guess)) {
        return p_air;
    }
    return 1/Mesh3d::at(pos, guess, vgrid);
}

double SlownessMesh3d::at(const Point3d& pos) const
{
    Index3d guess = nodeIndex(nearest(pos));
    return at(pos,guess);
}

double SlownessMesh3d::at(const Point3d& pos, Index3d& guess, Point3d& du) const
{
    upper_left(pos,guess);
    if (in_water(pos,guess)){
        du = Point3d(0,0,0);
        return p_water;
    }
    if (in_air(pos,guess)){
        du = Point3d(0,0,0);
        return p_air;
    }
    
    Point3d ratio = ratios(pos,guess);
    double  u     = interpolate(ratio, guess, vgrid);
    double r = ratio.x();
    double s = ratio.y();
    double t = ratio.z();
    double rr = 1-r;
    double ss = 1-s;
    double tt = 1-t;
    int i = guess.i();
    int j = guess.j();
    int k = guess.k();
    
    double dudt, dudtt;
    double dudx, dudy, dudz;
    
    if (xflat()) {
        dudt  = ss*vgrid(i,j,k+1) + s*vgrid(i,j+1,k+1);
        dudtt = ss*vgrid(i,j,k)   + s*vgrid(i,j+1,k);
    } else if (yflat()) {
        dudt  = rr*vgrid(i,j,k+1) + r*vgrid(i+1,j,k+1);
        dudtt = rr*vgrid(i,j,k)   + r*vgrid(i+1,j,k);
    } else {
        dudt  = (rr*ss*vgrid(i,j,k+1) + rr*s*vgrid(i,j+1,k+1)
                 + r*ss*vgrid(i+1,j,k+1) + r*s*vgrid(i+1,j+1,k+1));
        dudtt = (rr*ss*vgrid(i,j,k)   + r*ss*vgrid(i+1,j,k)
                 + rr*s*vgrid(i,j+1,k)   + r*s*vgrid(i+1,j+1,k));
    }
    
    // before here we've worked with velocity but we need the spatial partial derivatives of slowness: factor -1/(u*u)
    dudz = ((dudtt-dudt)/(u*u))/dz(k);
    if (xflat()) {
        dudx = 0;
        double duds  = tt*vgrid(i,j+1,k) + t*vgrid(i,j+1,k+1);
        double dudss = tt*vgrid(i,j,k)   + t*vgrid(i,j,k+1);
        dudy = -((duds-dudss)/(u*u)-by(i,j)*dudz)/dy(j);
    } else if (yflat()) {
        dudy = 0;
        double dudr  = tt*vgrid(i+1,j,k) + t*vgrid(i+1,j,k+1);
        double dudrr = tt*vgrid(i,j,k)   + t*vgrid(i,j,k+1);
        dudx = -((dudr-dudrr)/(u*u)-bx(i,j)*dudz)/dx(i);
    } else {
        double dudr  = (ss*tt*vgrid(i+1,j,k) + ss*t*vgrid(i+1,j,k+1) 
                        + s*tt*vgrid(i+1,j+1,k) + s*t*vgrid(i+1,j+1,k+1));
        double duds  = (rr*tt*vgrid(i,j+1,k) + rr*t*vgrid(i,j+1,k+1)
                        + r*tt*vgrid(i+1,j+1,k) + r*t*vgrid(i+1,j+1,k+1));
        double dudrr = (ss*tt*vgrid(i,j,k)   + s*tt*vgrid(i,j+1,k)
                        + ss*t*vgrid(i,j,k+1)   + s*t*vgrid(i,j+1,k+1));
        double dudss = (rr*tt*vgrid(i,j,k)   + r*tt*vgrid(i+1,j,k)
                        + rr*t*vgrid(i,j,k+1)   + r*t*vgrid(i+1,j,k+1));
        dudx = -((dudr-dudrr)/(u*u)-(s*bx(i,j)+ss*bx(i,j+1))*dudz)/dx(i);
        dudy = -((duds-dudss)/(u*u)-(r*by(i,j)+rr*by(i+1,j))*dudz)/dy(j);
    }
    du = Point3d(dudx,dudy,dudz);
    return 1/u;
}

void
SlownessMesh3d::cellNodes(int icell, std::vector<int>& n) const
{
    n[0] = cell_index[icell-1];
    Index3d index=node_index(n[0]);
    int i=index.i(), j=index.j(), k=index.k();
    if (xflat()) {
        n[1] = ser_index(i,j+1,k);
        n[2] = ser_index(i,j,k+1);
        n[3] = ser_index(i,j+1,k+1);
    } else if (yflat()) {
        n[1] = ser_index(i+1,j,k);
        n[2] = ser_index(i,j,k+1);
        n[3] = ser_index(i+1,j,k+1);
    } else {
        n[1] = ser_index(i+1,j,k);
        n[2] = ser_index(i,j+1,k);
        n[3] = ser_index(i,j,k+1);
        n[4] = ser_index(i,j+1,k+1);
        n[5] = ser_index(i+1,j,k+1);
        n[6] = ser_index(i+1,j+1,k);
        n[7] = ser_index(i+1,j+1,k+1);
    }
}

void SlownessMesh3d::cellNormKernel(int icell, Array2d<double>& T) const // Modified
{
    Index3d index=node_index(cell_index[icell-1]);
    double x = xflat() ? 1 : dx(index.i());
    double y = yflat() ? 1 : dy(index.j());
    double z = dz(index.k());
    T = T_common*sqrt(x*y*z);
}

void SlownessMesh3d::commonNormKernel() // Modified
{
    Array2d<double>& T = T_common;
    if (flat()) {
        T.resize(4,4);
        T(1,1) = 4; T(1,2) = 2; T(1,3) = 1; T(1,4) = 2;
        T(2,1) = 2; T(2,2) = 4; T(2,3) = 2; T(2,4) = 1;
        T(3,1) = 1; T(3,2) = 2; T(3,3) = 4; T(3,4) = 2;
        T(4,1) = 2; T(4,2) = 1; T(4,3) = 1; T(4,4) = 4;
        T *= (1.0/36.0);
    } else {
        T.resize(8,8);
        T(1,1) = 24; T(1,2) = 4;  T(1,3) = 4;  T(1,4) = 2;  T(1,5) = 4;  T(1,6) = 2;  T(1,7) = 2;  T(1,8) = 1;
        T(2,1) = 4;  T(2,2) = 24; T(2,3) = 2;  T(2,4) = 4;  T(2,5) = 2;  T(2,6) = 4;  T(2,7) = 1;  T(2,8) = 2;
        T(3,1) = 4;  T(3,2) = 2;  T(3,3) = 24; T(3,4) = 4;  T(3,5) = 2;  T(3,6) = 1;  T(3,7) = 4;  T(3,8) = 2;
        T(4,1) = 2;  T(4,2) = 4;  T(4,3) = 4;  T(4,4) = 24; T(4,5) = 1;  T(4,6) = 2;  T(4,7) = 2;  T(4,8) = 4;
        T(5,1) = 4;  T(5,2) = 2;  T(5,3) = 2;  T(5,4) = 1;  T(5,5) = 24; T(5,6) = 4;  T(5,7) = 4;  T(5,8) = 2;
        T(6,1) = 2;  T(6,2) = 4;  T(6,3) = 1;  T(6,4) = 2;  T(6,5) = 4;  T(6,6) = 24; T(6,7) = 2;  T(6,8) = 4;
        T(7,1) = 2;  T(7,2) = 1;  T(7,3) = 4;  T(7,4) = 2;  T(7,5) = 4;  T(7,6) = 2;  T(7,7) = 24; T(7,8) = 4;
        T(8,1) = 1;  T(8,2) = 2;  T(8,3) = 2;  T(8,4) = 4;  T(8,5) = 2;  T(8,6) = 4;  T(8,7) = 4;  T(8,8) = 24;
        
        T *= (1.0/216.0);
    }

    Array1d<double> p(T.nCol());
    d_choldc(T.toRecipe(), p.size(), p.toRecipe());
    for (int i=1; i <= p.size(); i++){
        T(i,i) = p(i);
        for (int j=i+1; j<=p.size(); j++){
            T(i,j) = T(j,i);
        }
    }
}

int SlownessMesh3d::locate_in_cell(const Point3d& p, Index3d& guess,
                                 opt_idx& jLFU, opt_idx& jRFU, opt_idx& jLBU,
                                 opt_idx& jLFL, opt_idx& jLBL, opt_idx& jRFL,
                                 opt_idx& jRBU, opt_idx& jRBL,
                                 Point3d& r) const
{
    upper_left(p,guess);
    if (in_water(p,guess)) {
        return -1;
    }
    if (in_air(p,guess)) {
        return -2;
    }
    
    int i=guess.i(), j=guess.j(), k=guess.k();

    jLFU = ser_index(i,j,k);
    jLFL = ser_index(i,j,k+1);
    if (!xflat()) {
        jRFU = ser_index(i+1,j,k);
        jRFL = ser_index(i+1,j,k+1);
    }
    if (!yflat()) {
        jLBU = ser_index(i,j+1,k);
        jLBL = ser_index(i,j+1,k+1);
    }
    if (!(xflat() || yflat())) {
        jRBU = ser_index(i+1,j+1,k);
        jRBL = ser_index(i+1,j+1,k+1);
    }
    
    r = ratios(p, guess);
    return index2cell(i,j,k);
}

namespace details {
    int
    near(double x, std::vector<double> const& pos) {
        if (pos.size() == 1) {
            return 1;
        } else {
            for (std::size_t i=1; i<pos.size(); i++){
                double halfway = pos[i-1] + 0.5*(pos[i] - pos[i-1]);
                if ( x < halfway ) {
                    return i;
                }
            }
            return pos.size();
        }
    }
}

int SlownessMesh3d::nearest(const Point3d& src) const // Modified
{
    // returns the serial index of the node nearest to
    // the given source location
    int inear = details::near(src.x(), xpos());
    int jnear = details::near(src.y(), ypos());
    int knear = details::near(src.z()-topo(inear,jnear), zpos());
    
    return ser_index(inear, jnear, knear);
}
void SlownessMesh3d::nearest(const Interface3d& itf, std::vector<int>& inodes) const // Modified
{
    int n=1;
    inodes.resize(nx()*ny());

    for (int i=1; i<=nx(); i++){
        for(int j=1; j<=ny(); j++){
            double x=xpos(i);
            double y=ypos(j);
            double z=itf.z(x,y);
            Point3d p(x,y,z);
            inodes[n-1] = nearest(p);
            n++;
        }
    }
}

void
SlownessMesh3d::outMesh(ostream& os) const {
    os << nx() << " " << ny() << " " << nz() << " "
       << 1.0/p_water << " " << 1.0/p_air << '\n';
    
    for (int i=1; i<=nx(); i++) {
        os << xpos(i) << " ";
    }
    os << '\n';
    for (int j=1; j<=ny(); j++) {
        os << ypos(j) << " ";
    }
    os << '\n';
    for (int i=1; i<=nx(); i++){
        for (int j=1; j<=ny(); j++){
            os << topo(i,j) << " ";
        }
        os << '\n';
    }
    for (int k=1; k<=nz(); k++) {
        os << zpos(k) << " ";
    }
    os << '\n';
    for (int i=1; i<=nx(); i++){
        for (int j=1; j<=ny(); j++){
            for (int k=1; k<=nz(); k++){
                os << 1.0/pgrid(i,j,k) << " ";
            }
            os << '\n';
        }
    }
}

bool SlownessMesh3d::in_water(const Point3d& pos, const Index3d& guess) const
{
    if (pos.z()<0.0) {
        return false; // negative z means above horizon
    } else {
        int i = guess.i();
        int j = guess.j();
        double zsurf = topo(i,j);

        double rx = 0; 
        double ry = 0;
        if (!xflat()) {
            zsurf += bx(i,j)*xratio(pos,i);
        }
        if (!yflat()) {
            zsurf += by(i,j)*yratio(pos,j);
        }
        if (!(xflat() || yflat())) {
            zsurf += (bx(i,j+1)-bx(i,j))*rx*ry;
        }
        return pos.z()<zsurf;
    }
}

bool SlownessMesh3d::in_air(const Point3d& pos, const Index3d& guess) const // Modified
{
    if (pos.z()>=0.0) {
        return false; // positive z means below horizon
    }
    
    int i=guess.i(), j=guess.j();
    double rx = xratio(pos, i);
    double ry = yratio(pos, j);
    double zsurf = topo(i,j);
    if (!xflat()) {
        zsurf += bx(i,j)*rx;
    }
    if (!yflat()) {
        zsurf += by(i,j)*ry;
    }
    if (!(xflat() || yflat())) {
        zsurf += (bx(i,j+1)-bx(i,j))*rx*ry;
    }
    return pos.z()<zsurf;
}

bool SlownessMesh3d::inWater(const Point3d& pos) const // Modified
{
    return in_water(pos,upper_left(pos));
}

bool SlownessMesh3d::inAir(const Point3d& pos) const // Modified
{
    return in_air(pos,upper_left(pos));
}

void SlownessMesh3d::printElements(ostream& os) const
{
    for (int i=1; i<nx(); i++){
        for (int j=1; j<ny(); j++){
            for (int k=1; k<nz(); k++){
                os << ">\n";
                os << xpos(i) << " " << ypos(j)<< " " << zpos(k)+topo(i,j) << '\n';
                if (!xflat()) {
                    os << xpos(i+1) << " " << ypos(j) << " " << zpos(k)+topo(i+1,j) << '\n';
                }
                if (!yflat()) {
                    os << xpos(i) << " " << ypos(j+1) << " " << zpos(k)+topo(i,j+1) << '\n';
                }
                os << xpos(i) << " " << ypos(j) << " " << zpos(k+1)+topo(i,j) << '\n';
                if (!(xflat() || yflat())) {
                    os << xpos(i+1) << " " << ypos(j+1) << " " << zpos(k)+topo(i+1,j+1) << '\n';
                }
                if (!xflat()) {
                    os << xpos(i+1) << " " << ypos(j) << " " << zpos(k+1)+topo(i+1,j) << '\n';
                }
                if (!yflat()) {
                    os << xpos(i) << " " << ypos(j+1) << " " << zpos(k+1)+topo(i,j+1) << '\n';
                }
                if (!(xflat() || yflat())) {
                    os << xpos(i+1) << " " << ypos(j+1) << " " << zpos(k+1)+topo(i+1,j+1) << '\n';
                }
            }
        }
    }
}

void SlownessMesh3d::printVGrid(ostream& os, bool printAW) const // Modified
{
    double min_topo=0.0;
    for (int i=1; i<=nx(); i++){
        for (int j=1; j<=ny(); j++){
            if (topo(i,j)<min_topo) min_topo=topo(i,j);
        }
    }
    min_topo -= 1.0;        // upper bound for courtesy grid making
    // (extra 1km)
    double dz2 = dz(1)/2;
    for (int i=1; i<=nx(); i++){
        for (int j=1; j<=ny(); j++){
            for (int k=1; k<=nz(); k++){
                os << xpos(i) << " "
                   << ypos(j) << " "
                   << zpos(k)+topo(i,j) << " "
                   << 1.0/pgrid(i,j,k) << '\n';
            }
            if (printAW){
                // some extra grid points for grid file 
                double z = topo(i,j)-dz2;
                while(z>=0){
                    os << xpos(i) << " "
                       << ypos(j) << " "
                       << z << " " << 1.0/p_water << '\n';
                    z -= dz2;
                }
                while(z>=min_topo){
                    os << xpos(i) << " "
                       << ypos(j) << " "
                       << z << " " << 1.0/p_air << '\n';
                    z -= dz2;
                }
            }
        }
    }
}

void SlownessMesh3d::printVGrid(ostream& os,
                                double x0, double x1,
                                double y0, double y1,
                                double z0, double z1,
                                double dx, double dy, double dz) const // Modified
{
    Index3d guess = nodeIndex(nearest(Point3d(x0,y0,z0)));

    int nxx = int((x1-x0)/dx+1);
    int nyy = int((y1-y0)/dy+1);
    int nzz = int((z1-z0)/dz+1);

    for (int ix=1; ix<=nxx; ix++){
        double x = x0+(ix-1)*dx;
        for (int iy=1; iy<=nyy; iy++){
            double y = y0+(iy-1)*dy;
            for (int iz=1; iz<=nzz; iz++){
                double z = z0+(iz-1)*dz;
                double p = at(Point3d(x,y,z),guess);
                os << x << " " << y << " " << z << " " << 1.0/p << '\n';
            }
        }
    }
}

void SlownessMesh3d::printMaskGrid(ostream& os,
                                   const std::vector<int>& valid_node) const // Modified, used?
{
    assert(valid_node.size() == nb_nodes());

    os << nx() << " " << ny() << " " << nz() << " "
       << 1.0/p_water << " " << 1.0/p_air << '\n';

    for (int i=1; i<=nx(); i++) os << xpos(i) << " ";
    os << '\n';
    for (int j=1; j<=ny(); j++) os << ypos(j) << " ";
    os << '\n';
    for (int i=1; i<=nx(); i++){
        for (int j=1; j<=ny(); j++){
            os << topo(i,j) << " ";
        }
    }
    os << '\n';
    for (int k=1; k<=nz(); k++) os << zpos(k) << " ";
    os << '\n';
    int inode=1;
    for (int i=1; i<=nx(); i++){
        for (int j=1; j<=ny(); j++){
            for (int k=1; k<=nz(); k++){
                double val=0.0;
                if (valid_node[inode-1]>0){
                    val = 1.0/pgrid(i,j,k);
                }
                os << val << " ";
                inode++;
            }
            os << '\n';
        }
        //    os << '\n';
    }
}

// output DWS
void SlownessMesh3d::printMaskGrid(ostream& os,
                                   const Array1d<double>& dws) const // Modified
{
    assert(dws.size() == nb_nodes());

    int inode=1;
    for (int i=1; i<=nx(); i++){
        double x=xpos(i);
        for (int j=1; j<=ny(); j++){
            double y=ypos(j), t=topo(i,j);
            for (int k=1; k<=nz(); k++){
                double z=zpos(k)+t;
                os << x << " " << y << " " << z << " " << dws(inode) << "\n";
                inode++;
            }
        }
    }
}
