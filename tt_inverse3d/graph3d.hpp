/*
 * graph3d.h based on graph.h - graph method related classes
 *
 * Adria Melendez and Jun Korenaga
 * Spring 2011
 */

#ifndef _TOMO_GRAPH3D_H_
#define _TOMO_GRAPH3D_H_

#include <list>
#include "array.hpp" // from mconv source
#include "geom3d.hpp"
#include "index3d.hpp"
#include "interface3d.hpp"

#include "heap_deque.hpp"

class GraphSolver3d {
public:
    GraphSolver3d(SlownessMesh3d const&, int xorder, int yorder, int zorder);
    GraphSolver3d(GraphSolver3d const& other);
    GraphSolver3d(GraphSolver3d const& other, SlownessMesh3d const&);

    int xOrder() const { return fs_xorder; }
    int yOrder() const { return fs_yorder; }
    int zOrder() const { return fs_zorder; }
  void solve(const Point3d& src, bool ani);
  void solve_refl(const Point3d& src, const Interface3d& itf, bool ani);
    void do_refl_downward(); 
    void limitRange(double xmin, double xmax, double ymin, double ymax);
    void delimitRange();
    void refine_if_long(double);
    double critLength() const { return crit_len; }
    
    void pickPath(const Point3d& rcv, Array1d<Point3d>& path) const;
    void pickPathRcvW(const Point3d& rcv, Array1d<Point3d>& path,
		      int&, int&) const;
    void pickReflPath(const Point3d& rcv, Array1d<Point3d>& path,
		      int&, int&) const;
    void pickReflPathRcvW(const Point3d& rcv, Array1d<Point3d>& path,
			  int&, int&, int&, int&) const;
    void pickPathSrcW(const Point3d& rcv, Array1d<Point3d>& path,
		      int&, int&) const;
    void pickPathSrcWRcvW(const Point3d& rcv, Array1d<Point3d>& path,
			  int&, int&, int&, int&) const;
    void pickReflPathSrcW(const Point3d& rcv, Array1d<Point3d>& path,
			  int&, int&, int&, int&) const;
    void pickReflPathSrcWRcvW(const Point3d& rcv, Array1d<Point3d>& path,
			      int&, int&, int&, int&, int&, int&) const;

    void printPath(ostream&) const;
    
    std::vector<int> const& previous_node() const { return my_prev_node;}
    int previous_node(int nid) const { return my_prev_node[nid-1];}

    std::vector<int> const& previous_node_up() const { return my_prev_node_up;}
    int previous_node_up(int nid) const { return my_prev_node_up[nid-1];}

    std::vector<int> const& previous_node_down() const { return my_prev_node_down;}
    int previous_node_down(int nid) const { return my_prev_node_down[nid-1];}

    Array1d<double> const& travel_time() const { return my_trav_time;}
    Array1d<double> const& travel_time_up() const { return my_trav_time_up;} 
    Array1d<double> const& travel_time_down() const { return my_trav_time_down;} 
    SlownessMesh3d const& slowness_mesh() const { return my_slowness_mesh; }

private:
    class trav_time_cmp {
    public:
        trav_time_cmp(Array1d<double> const& a) : data(a) { assert(&a); }
        static bool is_nan(double d) { return d!=d;}
        bool operator()(int i, int j) const { 
            double id = data(i);
            double jd = data(j);
            assert(!is_nan(id));
            assert(!is_nan(jd));
            return id > jd; 
        }
    private:
        Array1d<double> const& data;
    };

    class ForwardStar3d {
    public:
        ForwardStar3d(const Index3d& id, int ix, int iy, int iz);
        bool isIn(const Index3d&) const;
        
    private:
        const Index3d& orig;
        int xorder, yorder, zorder; //, dmin, dmid, dmax;
    };
 public:  
    typedef heap_deque<int,trav_time_cmp> work_list;
 private:

    void findRayEntryPoint(const Point3d& rcv, int& cur, 
			   const Array1d<double>& nodetime) const;
    void findRayExitPoint(const Point3d& src, int& cur) const;
    void process_travel_times(int src,
                              list<int>& nodes, 
                              Array1d<double>& trav_time, 
                              std::vector<int>& prev_node ) const;
    void process_travel_times(work_list& wl,
                              list<int>& nodes, 
                              Array1d<double>& trav_time, 
                              std::vector<int>& prev_node ) const;
    void crop_domain(int src, list<int>& points, Array1d<double>& trav_times, std::vector<int>& prev_nodes) const;
    void update_travel_times(ForwardStar3d const& fstar, int src, work_list& B,
                             Array1d<double>& trav_time, std::vector<int>& prev_node ) const;
    void compute_travel_times(ForwardStar3d const& fstar, int src, work_list& wl,  list<int>& nodes, 
                              Array1d<double>& trav_time, std::vector<int>& prev_node ) const;
    SlownessMesh3d const& my_slowness_mesh;
    int nnodes;
    bool is_solved, is_refl_solved, solve_refl_downward;
    Point3d src;
    int fs_xorder, fs_yorder, fs_zorder;
    std::vector<int> my_prev_node, my_prev_node_down, my_prev_node_up;
    Array1d<double> my_trav_time;
    Array1d<double> my_trav_time_up;
    Array1d<double> my_trav_time_down;
    const Interface3d* itfp;
    typedef list<int>::const_iterator nodeBrowser;
    typedef list<int>::iterator nodeIterator;
    const double local_search_radius;
    double xmin, xmax, ymin, ymax;
    bool my_refine_if_long;
    double crit_len;
  bool is_Ani;
};

#endif /* _TOMO_GRAPH3D_H_ */
