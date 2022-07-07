/*
 * smesh3d.hpp slowness mesh interface
 * based on smesh.h
 *
 * Adria Melendez and Jun Korenaga
 * Spring 2011
 */

#ifndef _TOMO_SMESH3D_H_
#define _TOMO_SMESH3D_H_

#include <list>
#include "boost/optional.hpp"
#include "boost/shared_ptr.hpp" //ALOUREIRO 11/02/2022
#include "array.hpp"
#include "mesh.hpp"
#include "geom3d.hpp"
#include "heap_deque.hpp"
#include "index3d.hpp"
#include "ani3d.hpp"

class Interface3d; // forward declaration

class SlownessMesh3d : public Mesh3d {

public:
  SlownessMesh3d(std::string const& fname);
  SlownessMesh3d(SlownessMesh3d const& other);

  void set(const Array1d<double>&);
  void get(Array1d<double>&) const;
  void vget(Array1d<double>&) const;
    
  const Index3d& nodeIndex(int i) const;
  int nodeIndex(int i, int j, int k) const;
  Point3d nodePos(int i) const;

  // slowness interpolation
  double at(const Point3d& pos) const;
  double at(const Point3d& pos, Index3d& guess) const;
  double at(const Point3d& pos, Index3d& guess, Point3d& du) const;
  double atWater() const { return p_water; }

  // for graph traveltime calculation
  double calc_ttime(int node_src, int node_rcv, bool ani) const; 

  void tiltAngle(double);
  double Phi;
  boost::shared_ptr<AnisotropyMesh3d> anid, anie;

  void add_anisotropy(boost::shared_ptr<AnisotropyMesh3d> const& anidp,
		      boost::shared_ptr<AnisotropyMesh3d> const& aniep);

    void usingVh();
    bool useVh = false;

  // cell-oriented functions
  int numCells() const { return cell_index.size(); }
  void cellNodes(int cell, std::vector<int>& neighb) const;
  void cellNormKernel(int, Array2d<double>&) const;

  typedef boost::optional<int> opt_idx;    
  int locate_in_cell( const Point3d&, Index3d&,
		      opt_idx&, opt_idx&, opt_idx&, opt_idx&, opt_idx&, opt_idx&, opt_idx&, opt_idx&, 
		      Point3d& r) const;

  // miscellaneous
  int nearest(const Point3d& src) const;
  void nearest(const Interface3d& itf, std::vector<int>& inodes) const;
  bool inWater(const Point3d& pos) const;
  bool inAir(const Point3d& pos) const;

  // output functions
  void outMesh(ostream&) const;
  void printElements(ostream&) const;
  void printVGrid(ostream&,bool) const;
  void printVGrid(ostream&,
		  double,double,double,double,double,double,double,double,double) const;
  void printMaskGrid(ostream&, const std::vector<int>&) const;
  void printMaskGrid(ostream&, const Array1d<double>&) const;

  // Not sure
  boost::optional<double> xexp(Point3d const& p, int i, int j, int k, double lx) const;
  boost::optional<double> yexp(Point3d const& p, int i, int j, int k, double ly) const;
  boost::optional<double> zexp(Point3d const& p, int nodeIndex, double lz) const;

  friend class Interface3d;
    
private:
  void commonNormKernel();
  bool in_water(const Point3d& pos, const Index3d& guess) const;
  bool in_air(const Point3d& pos, const Index3d& guess) const;
  //    double almost_exact_ttime(double v1, double v2,
  //			      double dpath) const;

  //  double Phi;
  //  boost::shared_ptr<AnisotropyMesh3d> anid, anie;

  int ncells;
  double p_water; /* slowness of water column */
  double p_air; /* slowness of air column */
  Array3d<double> pgrid, vgrid;
  Array3d<int> ser_index, index2cell;
  Array1d<Index3d> node_index;
  std::vector<int> cell_index;
  Array2d<double> Sm_H1, Sm_H2, Sm_V, T_common;
};

#endif /* _TOMO_SMESH3D_H_ */
