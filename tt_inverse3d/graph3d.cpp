/*
 * graph3d.cc based on graph.cc - graph solver implementation
 *
 * Adria Melendez and Jun Korenaga
 * Spring 2011
 */

#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <cassert>
#include "smesh3d.hpp"
#include "heap_deque.hpp"
#include "graph3d.hpp"

double const inf_trav_time = std::numeric_limits<double>::max();

GraphSolver3d::GraphSolver3d(SlownessMesh3d const& m, int xorder, int yorder, int zorder) // Modified
: my_slowness_mesh(m),
    nnodes(m.nb_nodes()),
    is_solved(false), is_refl_solved(false), solve_refl_downward(false),
    fs_xorder(xorder), fs_yorder(yorder), fs_zorder(zorder),
    my_prev_node(),
    my_trav_time(),
    my_trav_time_up(),
    my_trav_time_down(),
    local_search_radius(10.0),
    xmin(my_slowness_mesh.xmin()), xmax(my_slowness_mesh.xmax()),
    ymin(my_slowness_mesh.ymin()), ymax(my_slowness_mesh.ymax()),
    my_refine_if_long(false)
{}

GraphSolver3d::GraphSolver3d( GraphSolver3d const& g, SlownessMesh3d const& m)
: my_slowness_mesh(m),
    nnodes(g.nnodes),
    is_solved(g.is_solved),
    is_refl_solved(g.is_refl_solved),
    solve_refl_downward(g.solve_refl_downward),
    src(g.src),
    fs_xorder(g.fs_xorder),
    fs_yorder(g.fs_yorder),
    fs_zorder(g.fs_zorder),
    my_prev_node(g.my_prev_node),
    my_prev_node_down(g.my_prev_node_down),
    my_prev_node_up(g.my_prev_node_up),
    my_trav_time(g.my_trav_time),
    my_trav_time_up(g.my_trav_time_up),
    my_trav_time_down(g.my_trav_time_down),
    itfp(g.itfp),
    local_search_radius(g.local_search_radius),
    xmin(my_slowness_mesh.xmin()),
    xmax(my_slowness_mesh.xmax()),
    ymin(my_slowness_mesh.ymin()),
    ymax(my_slowness_mesh.ymax()),
    my_refine_if_long(g.my_refine_if_long),
    crit_len(g.crit_len) {}

void GraphSolver3d::limitRange(double x1, double x2, double y1, double y2)
{
    assert( x1 <= x2 );
    assert( y1 <= y2 );
    xmin = x1 - local_search_radius;
    xmax = x2 + local_search_radius;
    ymin = y1 - local_search_radius;
    ymax = y2 + local_search_radius;
}

void GraphSolver3d::delimitRange() // This function is not used.
{
    xmin = slowness_mesh().xmin();
    xmax = slowness_mesh().xmax();
}

void GraphSolver3d::refine_if_long(double x)
{
    my_refine_if_long = true;
    crit_len = x;
}

void
GraphSolver3d::crop_domain(int src, list<int>& points, Array1d<double>& trav_times, std::vector<int>& prev_nodes) const {
    assert(points.empty());
    assert(trav_times.size() == nnodes);
    assert(prev_nodes.size() == nnodes);
    
    for (int i=1; i<=nnodes; i++){
        if (i != src){
            double x = slowness_mesh().nodePos(i).x();
            double y = slowness_mesh().nodePos(i).y();
            if (x >= xmin && x <= xmax && y >= ymin && y <= ymax){
                points.push_back(i);
                trav_times(i) = inf_trav_time;
                prev_nodes[i-1] = src;
            }
        }
    }
}

void
GraphSolver3d::update_travel_times(ForwardStar3d const& fstar, 
                                   int src, work_list& wl, 
                                   Array1d<double>& trav_time, 
                                   std::vector<int>& prev_node ) const {
    assert(std::find(wl.begin(), wl.end(), src) == wl.end());
    // loop over [FS and B] to update traveltime
    //    cerr << "is_Ani = " << is_Ani << '\n';
    for(work_list::const_iterator it = wl.begin(); it != wl.end(); ++it) {            
        if (fstar.isIn(slowness_mesh().nodeIndex(*it))){
	  double tmp_trav_time=trav_time(src)+slowness_mesh().calc_ttime(src,*it,is_Ani);
            double orig_trav_time=trav_time(*it);
            if (orig_trav_time > tmp_trav_time){
                trav_time(*it) = tmp_trav_time;
                prev_node[*it-1] = src;
                wl.promote(it);
            }
        }
    } 
}

	
void
GraphSolver3d::compute_travel_times(ForwardStar3d const& fstar, 
                                    int src,
                                    work_list& wl, 
                                    list<int>& nodes, 
                                    Array1d<double>& trav_time, 
                                    std::vector<int>& prev_node ) const {
  //  cerr << "is_Ani = " << is_Ani << '\n';
    typedef list<int>::iterator liter;
    // loop over [FS and nodes] to calculate traveltime
    // and transfer them to wl
    std::vector<liter> to_erase;
    for( liter it = nodes.begin(); it != nodes.end(); ++it) {
        if (fstar.isIn(slowness_mesh().nodeIndex(*it))){
	  trav_time(*it) = trav_time(src)+slowness_mesh().calc_ttime(src,*it,is_Ani);
            prev_node[*it-1] = src;
            wl.push(*it);
            to_erase.push_back(it);
        }
    }
    for(std::vector<liter>::iterator it = to_erase.begin(); it != to_erase.end(); ++it) {
        nodes.erase(*it);
    }
}

void
GraphSolver3d::process_travel_times(work_list& wl,
                                    list<int>& nodes, 
                                    Array1d<double>& trav_time, 
                                    std::vector<int>& prev_node ) const {
    while(wl.size() != 0){
        // find a node with minimum traveltime from B,
        // (and transfer it to A)
        int N0 = wl.top();
        wl.pop();
        // form a forward star
        ForwardStar3d fstar(slowness_mesh().nodeIndex(N0), fs_xorder, fs_yorder, fs_zorder);
        update_travel_times(fstar, N0, wl, trav_time, prev_node );
	compute_travel_times(fstar, N0, wl, nodes, trav_time, prev_node );
    }
}

void
GraphSolver3d::process_travel_times(int src,
                                    list<int>& nodes, 
                                    Array1d<double>& trav_time, 
                                    std::vector<int>& prev_node ) const {
    work_list wl = work_list(trav_time_cmp(trav_time));
    wl.push(src); 
    process_travel_times(wl, nodes, trav_time, prev_node );
}

void GraphSolver3d::solve(const Point3d& s, bool ani)
{
    assert(!slowness_mesh().inAir(s));
    is_Ani = ani;
    //    cerr << "is_Ani = " << is_Ani << '\n';
    //    cerr << "ani = " << ani << '\n';
    src = s;

    if (slowness_mesh().inAir(src)){
      std::cout << "src(x,y,z)=(" << src.x() << "," << src.y() << "," << src.z() << ")\n";
      error("GraphSolver3d::solve - source in air: currently unsupported source configuration detected.\n");
    }

    int Nsrc;
    if (slowness_mesh().inWater(src)){
	int cur = slowness_mesh().nearest(src);
	findRayExitPoint(src,cur);
	Nsrc = cur;
    }else{
	Nsrc = slowness_mesh().nearest(src);
    }
    // initialization

    Array1d<double> trav_time(nnodes, std::numeric_limits<double>::quiet_NaN());
    std::vector<int> prev_node(nnodes, -1);    
    list<int> C;
    prev_node[Nsrc-1] = Nsrc;
    trav_time(Nsrc) = 0.0;

    crop_domain(Nsrc, C, trav_time, prev_node);
    process_travel_times(Nsrc, C, trav_time, prev_node );

    my_trav_time.stl().swap(trav_time.stl());
    my_prev_node.swap(prev_node);
    is_solved = true;
}

void GraphSolver3d::do_refl_downward(){ solve_refl_downward = true; }

void
GraphSolver3d::solve_refl(const Point3d& s, const Interface3d& itf, bool ani) {
    if (!is_solved && !solve_refl_downward) {
        error("GraphSolver3d::not ready for solve_refl()");
    }
    is_Ani = ani;
    //    cerr << "is_Ani = " << is_Ani << '\n';
    //    cerr << "ani = " << ani << '\n';
    std::vector<int> prev_node_down(nnodes, -1);
    std::vector<int> prev_node_up(nnodes, -1);
    Array1d<double> trav_time_down(nnodes, std::numeric_limits<double>::quiet_NaN());

    itfp = &itf;
    std::vector<int> itfnodes;
    slowness_mesh().nearest(itf,itfnodes);

    if (solve_refl_downward){
        // initialization (downgoing)
        src = s;
        if (slowness_mesh().inAir(src)){
	    std::cout << "src(x,y,z)=(" << src.x() << "," << src.y() << "," << src.z() << ")\n";
            error("GraphSolver3d::solve_refl - source in air: currently unsupported source configuration detected.\n");
        }

	int Nsrc;
	if (slowness_mesh().inWater(src)){
	    int cur = slowness_mesh().nearest(src);
	    findRayExitPoint(src,cur);
	    Nsrc = cur;
	}else{
	    Nsrc = slowness_mesh().nearest(src);
	}

        prev_node_down[Nsrc-1] = Nsrc;
        trav_time_down(Nsrc) = 0.0;
        
        list<int> C;
        for (int i=1; i<=nnodes; i++){
            if (i!=Nsrc){
                double x = slowness_mesh().nodePos(i).x();
                double y = slowness_mesh().nodePos(i).y();
                double z = slowness_mesh().nodePos(i).z();
                double itfz = itf.z(x,y);
                if (x>=xmin && x<=xmax && y>=ymin && y<=ymax){
                    bool is_ok = false;
                    if (z<=itfz){
                        is_ok = true;
                    }else{
                        for (int n=1; n<=itfnodes.size(); n++){
                            if (itfnodes[n-1]==i){
                                is_ok = true; break;
                            }
                        }
                    }
                    if (is_ok){
                        C.push_back(i);
                        trav_time_down(i) = inf_trav_time;
                        prev_node_down[i-1] = Nsrc;
                    }
                }
            }
        }
        process_travel_times(Nsrc, C, trav_time_down, prev_node_down );
    } else {
        // use refraction's graph solution
        prev_node_down = previous_node();
        trav_time_down = travel_time();
    }
    
    // initialization (upgoing)
    Array1d<double> trav_time_up(nnodes, std::numeric_limits<double>::quiet_NaN());
    trav_time_cmp B_comp(trav_time_up);
    work_list B_up(B_comp);
    list<int> C;
    for (int i=1; i<=nnodes; i++){
        bool is_itf = false;
        double x=slowness_mesh().nodePos(i).x();
        double y=slowness_mesh().nodePos(i).y();
        double z=slowness_mesh().nodePos(i).z();
        double itfz=itf.z(x,y);
        if (x>=xmin && x<=xmax && y>=ymin && y<=ymax){
            for (int n=1; n<=itfnodes.size(); n++){
                if (itfnodes[n-1]==i){
                    is_itf = true;
                    prev_node_up[i-1] = -1;
                    trav_time_up(i) = trav_time_down(i);
                    B_up.push(i);
                    break;
                }
            }
            if (!is_itf && z<=itfz){ // take nodes above the reflector
                trav_time_up(i) = inf_trav_time;
                prev_node_up[i-1] = 0;
                C.push_back(i);
            }
        }
    }
    process_travel_times(B_up, C, trav_time_up, prev_node_up );    

    my_prev_node_down.swap(prev_node_down);
    my_trav_time_down.swap(trav_time_down);
    my_prev_node_up.swap(prev_node_up);
    my_trav_time_up.swap(trav_time_up);
    is_refl_solved = true;
}

void GraphSolver3d::pickPath(const Point3d& rcv, Array1d<Point3d>& path) const
{
    assert(is_solved);
    
    list<Point3d> tmp_path;
    tmp_path.push_back(rcv);
    
    int cur = slowness_mesh().nearest(rcv);
    int prev;
    while ((prev = previous_node(cur)) != cur){
        Point3d p = slowness_mesh().nodePos(prev);
        tmp_path.push_back(p);
        cur = prev;
    }
    if (tmp_path.size()>1){
        tmp_path.pop_back();
    }
    if (tmp_path.size() == 1){
        // insert mid point for later bending refinement
        tmp_path.push_back(0.5*(src+rcv));
    }
    tmp_path.push_back(src);

    if (my_refine_if_long){
        list<Point3d>::iterator pt = tmp_path.begin();
        // skip the first point
        Point3d prev_p = *pt++;
        while(pt!=tmp_path.end()){
            double seg_len = prev_p.distance(*pt);
            if (seg_len>crit_len){
                int ndiv = int(seg_len/crit_len+1.0);
                double frac = 1.0/ndiv;
                for (int j=1; j<ndiv; j++){
                    double ratio = frac*j;
                    Point3d new_p = (1.0-ratio)*prev_p+ratio*(*pt);
                    tmp_path.insert(pt,new_p);
                }
            }
            prev_p = *pt++;
        }
    }
    
    path.resize(tmp_path.size());
    std::copy(tmp_path.begin(), tmp_path.end(), path.begin());
}

void GraphSolver3d::pickPathRcvW(const Point3d& rcv, Array1d<Point3d>& path,
				      int& i0, int& i1) const
{
    assert(is_solved);

    list<Point3d> tmp_path;
    tmp_path.push_back(rcv);

    int cur = slowness_mesh().nearest(rcv);
    findRayEntryPoint(rcv,cur,travel_time());
    Point3d entryp=slowness_mesh().nodePos(cur);
    tmp_path.push_back(entryp);
    tmp_path.push_back(entryp);
    tmp_path.push_back(entryp);
    i0 = 2; i1 = 4;
    
    int prev;
    int nadd=0;
    while ((prev=previous_node(cur)) != cur){
        Point3d p = slowness_mesh().nodePos(prev);
        tmp_path.push_back(p); nadd++;
        cur = prev;
    }
    
    if (nadd>0){
        tmp_path.pop_back();
    }
    if (nadd==1){
        // insert mid point for later bending refinement
        tmp_path.push_back(0.5*(src+entryp));
    }
    tmp_path.push_back(src);

    if (my_refine_if_long){
        list<Point3d>::iterator pt=tmp_path.begin();
        pt++; pt++; pt++; // now at the entry point
        Point3d prev_p = *pt;
        pt++; // skip the entry point
        while(pt!=tmp_path.end()){
            double seg_len = prev_p.distance(*pt);
            if (seg_len>crit_len){
                int ndiv = int(seg_len/crit_len+1.0);
                double frac = 1.0/ndiv;
                for (int j=1; j<ndiv; j++){
                    double ratio = frac*j;
                    Point3d new_p = (1.0-ratio)*prev_p+ratio*(*pt);
                    tmp_path.insert(pt,new_p);
                }
            }
            prev_p = *pt++;
        }
    }
    
    path.resize(tmp_path.size());
    std::copy(tmp_path.begin(), tmp_path.end(), path.begin());
}

void GraphSolver3d::pickReflPath(const Point3d& rcv, Array1d<Point3d>& path,
				 int& ir0, int& ir1) const
{

    if (!is_refl_solved) error("GraphSolver3d::not yet solved");

    list<Point3d> tmp_path;
    tmp_path.push_back(rcv);

    int cur = slowness_mesh().nearest(rcv);
    int prev;
    int nadd=0;
    while ((prev=previous_node_up(cur)) != -1){
        Point3d p = slowness_mesh().nodePos(prev);
        tmp_path.push_back(p); nadd++; 
        cur = prev;
    }
    double reflx = slowness_mesh().nodePos(cur).x();
    double refly = slowness_mesh().nodePos(cur).y();
    double reflz = itfp->z(reflx,refly);
    Point3d onrefl(reflx,refly,reflz);

    if (nadd>0){
        tmp_path.pop_back();
    }
    if (nadd==1){
        // insert mid point for later bending refinement
        tmp_path.push_back(0.5*(onrefl+rcv));
    }

    tmp_path.push_back(onrefl); ir0 = tmp_path.size();
    tmp_path.push_back(onrefl);
    tmp_path.push_back(onrefl); ir1 = tmp_path.size();
    
    nadd=0;
    while ((prev=previous_node_down(cur)) != cur){
        Point3d p = slowness_mesh().nodePos(prev);
        tmp_path.push_back(p); nadd++;
        cur = prev;
    }
    if (nadd>0){
        tmp_path.pop_back();
    }
    if (nadd==1){
        // insert mid point for later bending refinement
        tmp_path.push_back(0.5*(src+onrefl));
    }
    tmp_path.push_back(src);

    if (my_refine_if_long){
        int old_i=1, nadd_pre=0;
        list<Point3d>::iterator pt=tmp_path.begin();
        Point3d prev_p = *pt;
        pt++; old_i++; // skip the first point
        while(pt!=tmp_path.end()){
            double seg_len = prev_p.distance(*pt);
            if (seg_len>crit_len){
                int ndiv = int(seg_len/crit_len+1.0);
                double frac = 1.0/ndiv;
                for (int j=1; j<ndiv; j++){
                    double ratio = frac*j;
                    Point3d new_p = (1.0-ratio)*prev_p+ratio*(*pt);
                    tmp_path.insert(pt,new_p);
                    if (old_i<=ir0) nadd_pre++;
                }
            }
            prev_p = *pt++; old_i++;
        }
        ir0 += nadd_pre;
        ir1 += nadd_pre;
    }
    
    path.resize(tmp_path.size());
    std::copy(tmp_path.begin(), tmp_path.end(), path.begin());
}

void GraphSolver3d::pickReflPathRcvW(const Point3d& rcv, Array1d<Point3d>& path,
					  int& i0, int& i1, int& ir0, int& ir1) const
{
    if (!is_refl_solved) error("GraphSolver3d::not yet solved");

    list<Point3d> tmp_path;
    tmp_path.push_back(rcv);

    int cur = slowness_mesh().nearest(rcv);
    findRayEntryPoint(rcv,cur,travel_time_up());

    Point3d entryp=slowness_mesh().nodePos(cur);
    tmp_path.push_back(entryp);
    tmp_path.push_back(entryp);
    tmp_path.push_back(entryp);
    i0 = 2; i1 = 4;
    
    int prev;
    int nadd=0;
    while ((prev=previous_node_up(cur)) != -1){
        Point3d p = slowness_mesh().nodePos(prev);
        tmp_path.push_back(p); nadd++; 
        cur = prev;
    }
    double reflx = slowness_mesh().nodePos(cur).x();
    double refly = slowness_mesh().nodePos(cur).y();
    double reflz = itfp->z(reflx,refly);
    Point3d onrefl(reflx,refly,reflz);

    if (nadd>0){
        tmp_path.pop_back();
    }
    if (nadd==1){
        // insert mid point for later bending refinement
        tmp_path.push_back(0.5*(onrefl+entryp));
    }

    tmp_path.push_back(onrefl); ir0 = tmp_path.size();
    tmp_path.push_back(onrefl);
    tmp_path.push_back(onrefl); ir1 = tmp_path.size();
    
    nadd=0;
    while ((prev=previous_node_down(cur)) != cur){
        Point3d p = slowness_mesh().nodePos(prev);
        tmp_path.push_back(p); nadd++;
        cur = prev;
    }
    if (nadd>0){
        tmp_path.pop_back();
    }
    if (nadd==1){
        // insert mid point for later bending refinement
        tmp_path.push_back(0.5*(src+onrefl));
    }
    tmp_path.push_back(src);

    if (my_refine_if_long){
        int nadd_pre=0;
        list<Point3d>::iterator pt=tmp_path.begin();
        pt++; pt++; pt++; // now at the entry point
        Point3d prev_p = *pt;
        pt++; // skip the entry point
        int old_i = 5;
        while(pt!=tmp_path.end()){
            double seg_len = prev_p.distance(*pt);
            if (seg_len>crit_len){
                int ndiv = int(seg_len/crit_len+1.0);
                double frac = 1.0/ndiv;
                for (int j=1; j<ndiv; j++){
                    double ratio = frac*j;
                    Point3d new_p = (1.0-ratio)*prev_p+ratio*(*pt);
                    tmp_path.insert(pt,new_p);
                    if (old_i<=ir0) nadd_pre++;
                }
            }
            prev_p = *pt++; old_i++;
        }
        ir0 += nadd_pre;
        ir1 += nadd_pre;
    }

    path.resize(tmp_path.size());
    std::copy(tmp_path.begin(), tmp_path.end(), path.begin());
}

void GraphSolver3d::pickPathSrcW(const Point3d& rcv, Array1d<Point3d>& path,
				 int& is0, int& is1) const
{
    assert(is_solved);
    list<Point3d> tmp_path;
    tmp_path.push_back(rcv);
    
    int cur = slowness_mesh().nearest(rcv);
    int prev;
    while ((prev = previous_node(cur)) != cur){
        Point3d p = slowness_mesh().nodePos(prev);
        tmp_path.push_back(p);
        cur = prev;
    }
    if (tmp_path.size()>1){
        tmp_path.pop_back();
    }
    
    Point3d exitp=slowness_mesh().nodePos(cur);
    
    if (tmp_path.size() == 1){
        // insert mid point for later bending refinement
        tmp_path.push_back(0.5*(exitp+rcv));
    }
    // exit point: src on seafloor:
    tmp_path.push_back(exitp); is0 = tmp_path.size();
    tmp_path.push_back(exitp);
    tmp_path.push_back(exitp); is1 = tmp_path.size();
    // real source position:
    tmp_path.push_back(src);
    //     
    if (my_refine_if_long){
	int sold_i=1, snadd_pre=0;
	list<Point3d>::iterator pt=tmp_path.begin();
	Point3d prev_p = *pt;
	pt++; sold_i++; // skip the first point
	while(pt!=tmp_path.end()){
	    double seg_len = prev_p.distance(*pt);
	    if (seg_len>crit_len){
		int ndiv = int(seg_len/crit_len+1.0);
		double frac = 1.0/ndiv;
		for (int j=1; j<ndiv; j++){
		    double ratio = frac*j;
		    Point3d new_p = (1.0-ratio)*prev_p+ratio*(*pt);
		    tmp_path.insert(pt,new_p);
		    if (sold_i<=is0) snadd_pre++;
		}
	    }
	    prev_p = *pt++; sold_i++;
	}
	is0 += snadd_pre;
	is1 += snadd_pre;
    }
    path.resize(tmp_path.size());
    std::copy(tmp_path.begin(), tmp_path.end(), path.begin());
}

void GraphSolver3d::pickPathSrcWRcvW(const Point3d& rcv, Array1d<Point3d>& path,
				     int& i0, int& i1, int& is0, int& is1) const
{
    assert(is_solved);
    list<Point3d> tmp_path;
    tmp_path.push_back(rcv);
    
    int cur = slowness_mesh().nearest(rcv);
    findRayEntryPoint(rcv,cur,travel_time());
    Point3d entryp=slowness_mesh().nodePos(cur);
    
    tmp_path.push_back(entryp);
    tmp_path.push_back(entryp);
    tmp_path.push_back(entryp);
    i0 = 2; i1 = 4;
    
    int prev;
    int nadd=0;
    while ((prev=previous_node(cur)) != cur){
        Point3d p = slowness_mesh().nodePos(prev);
        tmp_path.push_back(p); nadd++;
        cur = prev;
    }
    
    if (nadd>0){
        tmp_path.pop_back();
    }
    
    Point3d exitp=slowness_mesh().nodePos(cur);
    
    if (nadd==1){
        // insert mid point for later bending refinement
        tmp_path.push_back(0.5*(exitp+entryp));
    }
    // exit point: src on seafloor:
    tmp_path.push_back(exitp); is0 = tmp_path.size();
    tmp_path.push_back(exitp);
    tmp_path.push_back(exitp); is1 = tmp_path.size();
    // real source position:
    tmp_path.push_back(src);
    
    if (my_refine_if_long){
	int snadd_pre=0;
	list<Point3d>::iterator pt=tmp_path.begin();
	pt++; pt++; pt++; // now at the entry point
	Point3d prev_p = *pt;
	pt++; // skip the entry point
	int sold_i = 5;
	while(pt!=tmp_path.end()){
	    double seg_len = prev_p.distance(*pt);
	    if (seg_len>crit_len){
		int ndiv = int(seg_len/crit_len+1.0);
		double frac = 1.0/ndiv;
		for (int j=1; j<ndiv; j++){
		    double ratio = frac*j;
		    Point3d new_p = (1.0-ratio)*prev_p+ratio*(*pt);
		    tmp_path.insert(pt,new_p);
		    if (sold_i<=is0) snadd_pre++;
		}
	    }
	    prev_p = *pt++; sold_i++;
	}
	is0 += snadd_pre;
	is1 += snadd_pre;
    }
    path.resize(tmp_path.size());
    std::copy(tmp_path.begin(), tmp_path.end(), path.begin());
}

void GraphSolver3d::pickReflPathSrcW(const Point3d& rcv, Array1d<Point3d>& path,
				     int& is0, int& is1, int& ir0, int& ir1) const
{
    if (!is_refl_solved) error("GraphSolver3d::not yet solved");
    list<Point3d> tmp_path;
    tmp_path.push_back(rcv);
    
    int cur = slowness_mesh().nearest(rcv);
    int prev;
    int nadd=0;
    while ((prev=previous_node_up(cur)) != -1){
        Point3d p = slowness_mesh().nodePos(prev);
        tmp_path.push_back(p); nadd++; 
        cur = prev;
    }
    double reflx = slowness_mesh().nodePos(cur).x();
    double refly = slowness_mesh().nodePos(cur).y();
    double reflz = itfp->z(reflx,refly);
    Point3d onrefl(reflx,refly,reflz);
    
    if (nadd>0){
        tmp_path.pop_back();
    }
    if (nadd==1){
        // insert mid point for later bending refinement
        tmp_path.push_back(0.5*(onrefl+rcv));
    }
    
    tmp_path.push_back(onrefl); ir0 = tmp_path.size();
    tmp_path.push_back(onrefl);
    tmp_path.push_back(onrefl); ir1 = tmp_path.size();
    
    nadd=0;
    while ((prev=previous_node_down(cur)) != cur){
        Point3d p = slowness_mesh().nodePos(prev);
        tmp_path.push_back(p); nadd++;
        cur = prev;
    }
    if (nadd>0){
        tmp_path.pop_back();
    }

    Point3d exitp=slowness_mesh().nodePos(cur);
    
    if (nadd==1){
        // insert mid point for later bending refinement
        tmp_path.push_back(0.5*(exitp+onrefl));
    }
    // exit point: src on seafloor:
    tmp_path.push_back(exitp); is0 = tmp_path.size();
    tmp_path.push_back(exitp);
    tmp_path.push_back(exitp); is1 = tmp_path.size();
    
    // real source position:
    tmp_path.push_back(src);
    
    if (my_refine_if_long){
	int old_i=1, sold_i=1, nadd_pre=0, snadd_pre=0;
	list<Point3d>::iterator pt=tmp_path.begin();
	Point3d prev_p = *pt;
	pt++; old_i++; sold_i++; // skip the first point
	while(pt!=tmp_path.end()){
	    double seg_len = prev_p.distance(*pt);
	    if (seg_len>crit_len){
		int ndiv = int(seg_len/crit_len+1.0);
		double frac = 1.0/ndiv;
		for (int j=1; j<ndiv; j++){
		    double ratio = frac*j;
		    Point3d new_p = (1.0-ratio)*prev_p+ratio*(*pt);
		    tmp_path.insert(pt,new_p);
		    if (old_i<=ir0) nadd_pre++;
		    if (sold_i<=is0) snadd_pre++;
		}
	    }
	    prev_p = *pt++; old_i++; sold_i++;
	}
	ir0 += nadd_pre;
	ir1 += nadd_pre;
	is0 += snadd_pre;
	is1 += snadd_pre;
    }
    path.resize(tmp_path.size());
    std::copy(tmp_path.begin(), tmp_path.end(), path.begin());
}

void GraphSolver3d::pickReflPathSrcWRcvW(const Point3d& rcv, Array1d<Point3d>& path,
					 int& i0, int& i1, int& is0, int& is1, int& ir0, int& ir1) const
{
    if (!is_refl_solved) error("GraphSolver3d::not yet solved");
    
    list<Point3d> tmp_path;
    tmp_path.push_back(rcv);
    
    int cur = slowness_mesh().nearest(rcv);
    findRayEntryPoint(rcv,cur,travel_time_up());
    Point3d entryp=slowness_mesh().nodePos(cur);
    
    tmp_path.push_back(entryp);
    tmp_path.push_back(entryp);
    tmp_path.push_back(entryp);
    i0 = 2; i1 = 4;
    
    int prev;
    int nadd=0;
    
    while ((prev=previous_node_up(cur)) != -1){
	Point3d p = slowness_mesh().nodePos(prev);
	tmp_path.push_back(p); nadd++;
	cur = prev;
    }
    
    double reflx = slowness_mesh().nodePos(cur).x();
    double refly = slowness_mesh().nodePos(cur).y();
    double reflz = itfp->z(reflx,refly);
    Point3d onrefl(reflx,refly,reflz);
    
    if (nadd>0){
	tmp_path.pop_back();
    }
    if (nadd==1){
	tmp_path.push_back(0.5*(onrefl+entryp));
    }
    
    tmp_path.push_back(onrefl); ir0 = tmp_path.size();
    tmp_path.push_back(onrefl);
    tmp_path.push_back(onrefl); ir1 = tmp_path.size();
    
    nadd=0;
    while ((prev=previous_node_down(cur)) != cur){
	Point3d p = slowness_mesh().nodePos(prev);
	tmp_path.push_back(p); nadd++;
	cur = prev;
    }	
    
    Point3d exitp=slowness_mesh().nodePos(cur);
    
    if (nadd>0){
	tmp_path.pop_back();
    }
    if (nadd==1){
	tmp_path.push_back(0.5*(exitp+onrefl));
    }
    
    tmp_path.push_back(exitp); is0 = tmp_path.size();
    tmp_path.push_back(exitp);
    tmp_path.push_back(exitp); is1 = tmp_path.size();
    
    tmp_path.push_back(src);
    
    if (my_refine_if_long){
	int nadd_pre=0, snadd_pre=0;
	list<Point3d>::iterator pt=tmp_path.begin();
	pt++; pt++; pt++; // now at the entry point
	Point3d prev_p = *pt;
	pt++; // skip the entry point
	int old_i = 5, sold_i = 5;
	while (pt!=tmp_path.end()){
	    double seg_len = prev_p.distance(*pt);
	    if (seg_len>crit_len){
		int ndiv = int(seg_len/crit_len+1.0);
		double frac = 1.0/ndiv;
		for (int j=1; j<ndiv; j++){
		    double ratio = frac*j;
		    Point3d new_p = (1.0-ratio)*prev_p+ratio*(*pt);
		    tmp_path.insert(pt,new_p);
		    if (old_i<=ir0) nadd_pre++;
		    if (sold_i<=is0) snadd_pre++;
		}
	    }
	    prev_p = *pt++; old_i++; sold_i++;
	}
	ir0 += nadd_pre;
	ir1 += nadd_pre;
	is0 += snadd_pre;
	is1 += snadd_pre;
    }
    path.resize(tmp_path.size());
    std::copy(tmp_path.begin(), tmp_path.end(), path.begin());
}

void GraphSolver3d::findRayEntryPoint(const Point3d& rcv, int& cur, const Array1d<double>& nodetime) const
{
    Point3d p0 = slowness_mesh().nodePos(cur);
    double pwater = slowness_mesh().atWater();
    double ttime0 = pwater*rcv.distance(p0)+nodetime(cur);
    Index3d rindex = slowness_mesh().nodeIndex(cur);
    int rnode = slowness_mesh().nodeIndex(rindex.i(),rindex.j(),1);
    Point3d rp = slowness_mesh().nodePos(rnode); // Projection of the rcv to the seafloor.

    // local grid search for best ray entry point (+/- 10km radius)
    Point3d pxmin(rcv.x()-local_search_radius,rcv.y(),rcv.z());
    Point3d pxmax(rcv.x()+local_search_radius,rcv.y(),rcv.z());
    Point3d pymin(rcv.x(),rcv.y()-local_search_radius,rcv.z());
    Point3d pymax(rcv.x(),rcv.y()+local_search_radius,rcv.z());
    Index3d indxmin = slowness_mesh().nodeIndex(slowness_mesh().nearest(pxmin));
    Index3d indxmax = slowness_mesh().nodeIndex(slowness_mesh().nearest(pxmax));
    Index3d indymin = slowness_mesh().nodeIndex(slowness_mesh().nearest(pymin));
    Index3d indymax = slowness_mesh().nodeIndex(slowness_mesh().nearest(pymax));

    for (int i=indxmin.i()+1; i<indxmax.i(); i++){
        for (int j=indymin.j()+1; j<indymax.j(); j++){
            int inode = slowness_mesh().nodeIndex(i,j,1);
            Point3d ip = slowness_mesh().nodePos(inode);
            if (rp.distance(ip)<=local_search_radius){ // Search within a circumference.
                double tmp_ttime = pwater*rcv.distance(ip)+nodetime(inode);
                if (tmp_ttime < ttime0){
                    ttime0 = tmp_ttime;
                    cur = inode; 
                }
            }
        }
    }
}

void GraphSolver3d::findRayExitPoint(const Point3d& src, int& cur) const
{
    Point3d p0 = slowness_mesh().nodePos(cur);
    double pwater = slowness_mesh().atWater();
    double ttime0 = pwater*src.distance(p0);
    Index3d sindex = slowness_mesh().nodeIndex(cur);
    int snode = slowness_mesh().nodeIndex(sindex.i(),sindex.j(),1);
    Point3d sp = slowness_mesh().nodePos(snode); // Projection of the src to the seafloor.
	//
    // local grid search for best ray exit point (+/- 10km radius)
    Point3d pxmin(src.x()-local_search_radius,src.y(),src.z());
    Point3d pxmax(src.x()+local_search_radius,src.y(),src.z());
    Point3d pymin(src.x(),src.y()-local_search_radius,src.z());
    Point3d pymax(src.x(),src.y()+local_search_radius,src.z());
    Index3d indxmin = slowness_mesh().nodeIndex(slowness_mesh().nearest(pxmin));
    Index3d indxmax = slowness_mesh().nodeIndex(slowness_mesh().nearest(pxmax));
    Index3d indymin = slowness_mesh().nodeIndex(slowness_mesh().nearest(pymin));
    Index3d indymax = slowness_mesh().nodeIndex(slowness_mesh().nearest(pymax));
    
    for (int i=indxmin.i()+1; i<indxmax.i(); i++){
        for (int j=indymin.j()+1; j<indymax.j(); j++){
            int inode = slowness_mesh().nodeIndex(i,j,1);
            Point3d ip = slowness_mesh().nodePos(inode);
            if (sp.distance(ip)<=local_search_radius){ // Search within a circumference.
                double tmp_ttime = pwater*src.distance(ip);
                if (tmp_ttime < ttime0){
                    ttime0 = tmp_ttime;
                    cur = inode;
                }
            }
        }
    }
}

void GraphSolver3d::printPath(ostream& os) const
{
    for (int cur=1; cur<=nnodes; cur++){
        int prev = previous_node(cur);
        Point3d cur_p = slowness_mesh().nodePos(cur);
        Point3d prev_p = slowness_mesh().nodePos(prev);
        os << ">\n"
           << cur_p.x() << " "  << cur_p.y() << " " << cur_p.z() << "\n"
           << prev_p.x() << " " << prev_p.y() << " " << prev_p.z() << '\n';
    }
}

//
// forward star
//
GraphSolver3d::ForwardStar3d::ForwardStar3d(const Index3d& id, int ix, int iy, int iz)
    : orig(id), xorder(ix), yorder(iy), zorder(iz) {
    assert(ix>0 && iy>0 && iz>0);
}

bool
GraphSolver3d::ForwardStar3d::isIn(const Index3d& node) const
{
    int d[3];
    d[0] = std::abs(orig.i()-node.i());
    d[1] = std::abs(orig.j()-node.j());
    d[2] = std::abs(orig.k()-node.k());

    if (d[0]>xorder || d[1]>yorder || d[2]>zorder) return false;
    if (d[0]==1 || d[1]==1 || d[2]==1) return true;

    sort(d, d + 3);
    int dmin=d[0];
    int dmid=d[1];
    int dmax=d[2];

    if (dmin>0){
        int mod_midmin = dmid%dmin;
        int mod_maxmin = dmax%dmin;
        if (mod_maxmin==0 && mod_midmin==0){
            return false;
        }else{
            return true;
        }
    }else{
        if (dmid==dmin){
            return false;  
        }else{
            int mod_maxmid = dmax%dmid;
            if (mod_maxmid==0){
                return false;
            }else{
                return true;
            }
        }
    }
}

