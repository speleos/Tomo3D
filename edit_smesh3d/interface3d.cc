/*
 * interface3d.cc based on interface.cc
 *
 * Adria Melendez and Jun Korenaga
 * Spring 2011
 */

#include "interface3d.h"
#include <iostream>
#include <fstream>
#include "error.h"
#include "util3d.h"

const double Interface3d::eps = 1e-4;

Interface3d::Interface3d(const SlownessMesh3d& smesh) //modified
{
    nxr=smesh.xpos.size();
    nyr=smesh.ypos.size();
    xpos.resize(nxr); ypos.resize(nyr);
    zpos.resize(nxr,nyr);
    xpos = smesh.xpos;
    ypos = smesh.ypos;
    zpos = smesh.topo;
    calc_slope();
}

Interface3d::Interface3d(const char* fn) //modified
{
//    int nxr, nyr;
    ifstream in(fn);
    if(!fn){
	cerr<< "Interface3d::cannot open "<<fn<<"\n";
	exit(1);
    }

    // First option: Input interface file as an xyz file
    // Row 1: number of x and y nodes

    in >> nxr >> nyr;
    xpos.resize(nxr);
    ypos.resize(nyr);
    zpos.resize(nxr,nyr);
    nnodes=nxr*nyr;
    ncells=(nxr-1)*(nyr-1);
    ser_index_r.resize(nxr,nyr);
    node_index_r.resize(nnodes);
    cell_index_r.resize(ncells);
    index2cell_r.resize(nxr-1,nyr-1);
    int N=1;
    for(int i=1; i<=nxr; i++){
	for(int j=1; j<=nyr; j++){
	    in >> xpos(i) >> ypos(j) >> zpos(i,j);
	    ser_index_r(i,j)=N;
	    node_index_r(N).set(i,j);
	    N++;
	}
    }

    int icell=1;
    for (int i=1; i<nxr; i++){
	for(int j=1; j<nyr; j++){
	    index2cell_r(i,j) = icell;
	    cell_index_r(icell++) = ser_index_r(i,j);
	}
    }


    // Second option: Input interface file based on velocity file format
    // Row 1: number of x and y nodes, Row 2: x, Row 3: y, followed by nxr rows with nyr elements for z

    //    in >> nxr >> nyr;
    //    xpos.resize(nxr);
    //    ypos.resize(nyr);
    //    zpos.resize(nxr,nyr);
    //    for(int i=1; i<=nxr; i++) in >> xpos(i);
    //    for(int j=1; j<=nyr; j++) in >> ypos(j);
    //    for(int i=1; i<=nxr; i++){
    //      for(int j=1; j<=nyr; j++) in >> zpos(i,j);
    //    }

    // sanity check
    for (int i=2; i<=nxr; i++){
	if (xpos(i)<=xpos(i-1)){
	    cerr << "Interface3d: illegal ordering of x nodes at i=" << i << ", positions: " << xpos(i-1) << ", "<< xpos(i) << ".\n";
	    exit(1);
	}
    }
    for (int j=2; j<=nyr; j++){
	if (ypos(j)<=ypos(j-1)){
	    cerr << "Interface3d: illegal ordering of y nodes at j=" << j << ", positions: " << ypos(j-1) << ", "<< ypos(j) << ".\n";
	    exit(1);
	}
    }
    calc_slope();
}


const Index2d& Interface3d::nodeIndex(int i) const { return node_index_r(i); }

int Interface3d::nodeIndex(int i, int j) const { return ser_index_r(i,j); } // Added

void Interface3d::calc_slope() //modified
{
    nxr=xpos.size();
    nyr=ypos.size();  

    slope_x.resize(nxr-1,nyr);
    slope_y.resize(nxr,nyr-1);
    for (int j=1; j<=nyr; j++){
	for (int i=1; i<nxr; i++) slope_x(i,j) = (zpos(i+1,j)-zpos(i,j))/(xpos(i+1)-xpos(i));
    }
    for (int i=1; i<=nxr; i++){
	for (int j=1; j<nyr; j++) slope_y(i,j) = (zpos(i,j+1)-zpos(i,j))/(ypos(j+1)-ypos(j));
    }
}

double Interface3d::xmin() const { return xpos.front(); }
double Interface3d::xmax() const { return xpos.back(); }
double Interface3d::ymin() const { return ypos.front(); }
double Interface3d::ymax() const { return ypos.back(); }

double Interface3d::z(double x, double y) const //modified.
{
    int ii=0, jj=0;
    double xratio=0., yratio=0.;
    int nxr=xpos.size();
    int nyr=ypos.size();

    if (x<xpos.front()){
	ii = 1; xratio = 0.;
    }else if (x>xpos.back()){
	ii = nxr-1; xratio = 1.;
    }else{
	for (int i=1; i<=nxr; i++){
	    if(fabs(x-xpos(i))<eps){
		if(i!=nxr){
		    ii = i; xratio = 0.;
		    break;
		}else{
		    ii=nxr-1; xratio=1.;
		    break;
		}
	    }
	}
    }
    if(ii==0){
	for (int i=1; i<nxr; i++){
	    if(x>xpos(i) && x<xpos(i+1)){
		ii = i;
		xratio = (x-xpos(i))/(xpos(i+1)-xpos(i)); // or dx = x-xpos(i);
		break;
	    }
	}
    }
    
    if (y<ypos.front()){
	jj = 1; yratio = 0.;
    }else if (y>ypos.back()){
	jj = nyr-1; yratio = 1.;
    }else{
	for (int j=1; j<=nyr; j++){
	    if(fabs(y-ypos(j))<eps){
		if(j!=nyr){
		    jj=j; yratio=0.;
		    break;
		}else{
		    jj=nyr-1; yratio=1.;
		    break;
		}
	    }
	}
    }
    if(jj==0){
	for (int j=1; j<nyr; j++){
	    if(y>ypos(j) && y<ypos(j+1)){
		jj = j;
		yratio = (y-ypos(j))/(ypos(j+1)-ypos(j)); // or dy = y-ypos(j);
		break;
	    }
	}
    }
    
    return zpos(ii,jj)*(1-xratio-yratio+xratio*yratio)+zpos(ii+1,jj)*(xratio-xratio*yratio)
	+zpos(ii,jj+1)*(yratio-xratio*yratio)+zpos(ii+1,jj+1)*xratio*yratio;

    error("Interface3d::z - impossible!");
}

double Interface3d::dzdx(double x, double y) const
{
    double z1=z(x,y);
    int nxr=xpos.size();
    double zp,zn,slope_zxp,slope_zxn;
    if(x<xpos.front() || x>xpos.back()){
	return 0.0; // out of bounds
    }else if(fabs(x-xpos(nxr))<eps){
	zp=z(xpos(nxr-1),y);
	slope_zxp=(z1-zp)/(x-xpos(nxr-1));
	return slope_zxp;
    }else if(fabs(x-xpos(1))<eps){
	zn=z(xpos(2),y);
	slope_zxn=(zn-z1)/(xpos(2)-x);
	return slope_zxn;
    }else{
	for(int i=2; i<nxr; i++){
	    if(fabs(x-xpos(i))<eps){
		zp=z(xpos(i-1),y);
		zn=z(xpos(i+1),y);
		slope_zxp=(z1-zp)/(x-xpos(i-1));
		slope_zxn=(zn-z1)/(xpos(i+1)-x);
		return 0.5*(slope_zxp+slope_zxn);
	    }
	}
	for(int i=1; i<nxr; i++){
	    if(x>xpos(i) && x<xpos(i+1)){
		zp=z(xpos(i),y);
		zn=z(xpos(i+1),y);
		slope_zxp=(zn-zp)/(xpos(i+1)-xpos(i));
		return slope_zxp;
	    }
	}
    }
    error("Interface3d::dzdx - impossible!");
}

double Interface3d::dzdy(double x, double y) const
{
    double z1=z(x,y);
    int nyr=ypos.size();
    double zp,zn,slope_zyp,slope_zyn;
    if(y<ypos.front() || y>ypos.back()){
	return 0.0; // out of bounds
    }else if(fabs(y-ypos(nyr))<eps){
	zp=z(x,ypos(nyr-1));
	slope_zyp=(z1-zp)/(y-ypos(nyr-1));
	return slope_zyp;
    }else if(fabs(y-ypos(1))<eps){
	zn=z(x,ypos(2));
	slope_zyn=(zn-z1)/(ypos(2)-y);
	return slope_zyn;
    }else{
	for(int j=2; j<nyr; j++){
	    if(fabs(y-ypos(j))<eps){
		zp=z(x,ypos(j-1));
		zn=z(x,ypos(j+1));
		slope_zyp=(z1-zp)/(y-ypos(j-1));
		slope_zyn=(zn-z1)/(ypos(j+1)-y);
		return 0.5*(slope_zyp+slope_zyn);
	    }
	}
	for(int j=1; j<nyr; j++){
	    if(y>ypos(j) && y<ypos(j+1)){
		zp=z(x,ypos(j));
		zn=z(x,ypos(j+1));
		slope_zyp=(zn-zp)/(ypos(j+1)-ypos(j));
		return slope_zyp;
	    }
	}
    }
    error("Interface3d::dzdy - impossible!");
}

void Interface3d::locateInPlane(double x, double y, int& j1, int& j2, int& j3, int& j4) const
{
    int nxr=xpos.size();
    int nyr=ypos.size();
    int inode=0;
    int jnode=0;

    for (int i=1; i<=nxr; i++){
	if(fabs(x-xpos(i))<eps){
	    if(i!=nxr){
		inode=i;
		break;
	    }else{
		inode=nxr-1;
		break;
	    }
	}
    }

    if(inode==0){
	for (int i=1; i<nxr; i++){
	    double xnode=xpos(i);
	    double xnode2=xpos(i+1);
	    if (x >= xnode && x < xnode2){
		inode = i;
		break;
	    }
	}
    }

    if(inode==0){
        if (x<xpos.front() || x>xpos.back()){ // out of x-bounds
	    j1=j2=j3=j4=-1;
	    return;
	}else{
	    error("Interface3d::locateInPlane - unable to locate x");
	}
    }

    for (int j=1; j<=nyr; j++){
	if(fabs(y-ypos(j))<eps){
	    if(j!=nyr){
		jnode=j;
		break;
	    }else{
		jnode=nyr-1;
		break;
	    }
	}
    }

    if(jnode==0){
	for (int j=1; j<nyr; j++){
	    double ynode=ypos(j);
	    double ynode2=ypos(j+1);
	    if (y >= ynode && y < ynode2){
		jnode = j;
		break;
	    }
	}
    }

    if(jnode==0){
	if (y<ypos.front() || y>ypos.back()){ // out of y-bounds
	    j1=j2=j3=j4=-1;
	    return;
	}else{
	    error("Interface3d::locateInPlane - unable to locate y");
	}
    }
    // serial index array NEED TO DEFINE SERIAL INDEX AND NODE INDEX ARRAYS.
    j1=ser_index_r(inode,jnode);
    j2=ser_index_r(inode+1,jnode);
    j3=ser_index_r(inode,jnode+1);
    j4=ser_index_r(inode+1,jnode+1);
}

double Interface3d::x(int i) const
{
    Index2d node=node_index_r(i);
    int inode=node.i();
    if (inode>0 && inode<=xpos.size()) return xpos(inode);
    error("Interface3d::x - subscript out of range");
}

double Interface3d::y(int i) const
{
    Index2d node=node_index_r(i);
    int jnode=node.j();
    if (jnode>0 && jnode<=ypos.size()) return ypos(jnode);
    error("Interface3d::y - subscript out of range");
}

void Interface3d::cellNodes(int icell, int& j1, int& j2, int& j3, int& j4) const // Added
{
    j1 = cell_index_r(icell);
    Index2d index=node_index_r(j1);
    int i=index.i(), j=index.j();
    j2 = ser_index_r(i+1,j);
    j3 = ser_index_r(i,j+1);
    j4 = ser_index_r(i+1,j+1);
}

void Interface3d::cellNormKernel(int icell, double& fac) const // Added
{
    Index2d index=node_index_r(cell_index_r(icell));
    int i=index.i(), j=index.j();
    double dx=xpos(i+1)-xpos(i);
    double dy=ypos(j+1)-ypos(j);
    fac = sqrt(dx*dy);
}

void Interface3d::set(const Array1d<double>& a)
{
    if (a.size() != zpos.nRow()*zpos.nCol())
	error("Interface3d::set - size mismatch");
    int Nr=1;
    for (int i=1; i<=zpos.nRow(); i++){
	for (int j=1; j<=zpos.nCol(); j++){
	    zpos(i,j) = a(Nr);
	    Nr++;
	}
    }
    calc_slope();
}

void Interface3d::get(Array1d<double>& a) const
{
    if (a.size() != zpos.nRow()*zpos.nCol())
	error("Interface3d::get - size mismatch");
    int Nr=1;
    for (int i=1; i<=zpos.nRow(); i++){
	for (int j=1; j<=zpos.nCol(); j++){
	    a(Nr) = zpos(i,j);
	    Nr++;
	}
    }
}

ostream&
operator<<(ostream& os, const Interface3d& itf)
{
    os << itf.xpos.size() <<" "<< itf.ypos.size() <<'\n';
    for (int i=1; i<=itf.xpos.size(); i++){
	for (int j=1; j<=itf.ypos.size(); j++){
	    os << itf.xpos(i) << " " << itf.ypos(j) << " " <<  itf.zpos(i,j) << '\n';
	}
    }
    return os;
}
    
