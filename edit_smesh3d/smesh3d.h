/*
 * smesh3d.h - slowness mesh interface
 * based on smesh.h
 *
 * Adria Melendez and Jun Korenaga
 * Spring 2011
 */

#ifndef _TOMO_SMESH3D_H_
#define _TOMO_SMESH3D_H_

#include <list>
#include "array.h" // from mconv source
#include "geom3d.h"
#include "heap_deque.h"
#include "index3d.h"

class Interface3d; // forward declaration

class SlownessMesh3d {
 public:
  SlownessMesh3d(const char*); // file input

  void set(const Array1d<double>&);
  void get(Array1d<double>&) const;
  void vget(Array1d<double>&) const;
    
  // general bookkeeping functions
  int numNodes() const { return nnodes; }
  int Nx() const { return nx; }
  int Ny() const { return ny; }
  int Nz() const { return nz; }
  double xmin() const;
  double xmax() const;
  double ymin() const;
  double ymax() const;
  double zmin() const;
  double zmax() const;
  const Index3d& nodeIndex(int i) const;
  int nodeIndex(int i, int j, int k) const;
  Point3d nodePos(int i) const;

  // slowness interpolation
  double at(const Point3d& pos) const;
  double at(const Point3d& pos, Index3d& guess) const;
  double at(const Point3d& pos, Index3d& guess,
	    double& dudx, double& dudy,  double& dudz) const;
  double atWater() const { return p_water; }

  // for graph traveltime calculation
  double calc_ttime(int node_src, int node_rcv) const; 
									   
  // cell-oriented functions
  int numCells() const { return ncells; }
  void cellNodes(int, int&, int&, int&, int&, int&, int&, int&, int&) const;
  int cellIndex(int i, int j, int k) const { return index2cell(i,j,k); }
  void cellGradientKernel(int, Array2d<double>&, double, double) const;
  void cellNormKernel(int, Array2d<double>&) const;
  int locateInCell(const Point3d&, Index3d&,
		   int&, int&, int&, int&, int&, int&, int&, int&, 
		   double&, double&, double&, double&, double&, double&) const;
  void nodalCellVolume(Array1d<double>&,
		       Array1d<double>&, Array1d<double>&, Array1d<Point3d>&) const; // for gravity
    
  // miscellaneous
  int nearest(const Point3d& src) const;
  void nearest(const Interface3d& itf, Array1d<int>& inodes) const;
  bool inWater(const Point3d& pos) const;
  bool inAir(const Point3d& pos) const;

  // output functions
  void outMesh(ostream&) const;
  void printElements(ostream&) const;
  void printVGrid(ostream&,bool) const;
  void printVGrid(ostream&,
		  double,double,double,double,double,double,double,double,double) const;
  void printMaskGrid(ostream&, const Array1d<int>&) const;
  void printMaskGrid(ostream&, const Array1d<double>&) const;

  friend class Interface3d;
    
 private:
  void upperleft(const Point3d& pos, Index3d& guess) const;
  void calc_local(const Point3d& pos, int, int, int,
		  double&, double&, double&, double&, double&, double&) const;
  void commonGradientKernel();
  void commonNormKernel();
  bool in_water(const Point3d& pos, const Index3d& guess) const;
  bool in_air(const Point3d& pos, const Index3d& guess) const;
  //    double almost_exact_ttime(double v1, double v2,
  //			      double dpath) const;

  int nx, ny, nz, nnodes, ncells;
  double p_water; /* slowness of water column */
  double p_air; /* slowness of air column */
  Array3d<double> pgrid, vgrid;
  Array3d<int> ser_index, index2cell;
  Array1d<Index3d> node_index;
  Array1d<int> cell_index;
  Array1d<double> xpos, ypos, zpos;
  Array2d<double> topo, bx_vec, by_vec;
  Array1d<double> rdx_vec, rdy_vec, rdz_vec;
  Array1d<double> dx_vec, dy_vec, dz_vec;
  Array2d<double> Sm_H1, Sm_H2, Sm_V, T_common;
  const double eps;
};

#endif /* _TOMO_SMESH3D_H_ */
