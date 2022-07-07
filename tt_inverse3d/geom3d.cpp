/*
 * geom3d.cc based on geom.cc
 *
 * Adria Melendez and Jun Korenaga
 * Spring 2011
 */

#include <iostream>
#include <iterator>
#include "error.hpp"
#include "geom3d.hpp"

// Point3d
Point3d::Point3d(double x, double y, double z) //modified
{
    v[0] = x;
    v[1] = y;
    v[2] = z;
}

double Point3d::inner_product(const Point3d& p) const  //modified
{
    double ip = v[0]*p.v[0]+v[1]*p.v[1]+v[2]*p.v[2];
    return ip;
}

Point3d Point3d::outer_product(const Point3d& p) const  //new
{
    double cx = v[1]*p.v[2]-v[2]*p.v[1];
    double cy = v[2]*p.v[0]-v[0]*p.v[2];
    double cz = v[0]*p.v[1]-v[1]*p.v[0];

    Point3d c(cx,cy,cz);
    return c;
}

std::ostream&
operator<<(std::ostream& out, Point3d const& p) {
    out << '(' << p.x() << ',' << p.y() << ',' << p.z() << ')';
    return out;
}

std::ostream&
operator<<(std::ostream& out, std::vector<Point3d> const& p) {
    out << '[' << p.size() << ':';
    std::copy(p.begin(), p.end(), std::ostream_iterator<Point3d const&>(out, ""));
    out << ']';
    return out;
}

RegularDomain3d::RegularDomain3d(double x1, double x2, double y1, double y2)
{
    if (x1>x2 || y1>y2) error("RegularDomain3d::BadInput"); // i range = 0-1
    
    vmin[0] = x1;
    vmax[0] = x2;
    vmin[1] = y1;
    vmax[1] = y2;
}

// RegularBC3d
RegularBC3d::RegularBC3d(int i, double x1, double x2, double y1, double y2, double d)
    : idof(i), RegularDomain3d(x1,x2,y1,y2), val_(d)
{
    // the following is commented out for VectorMesh2d_plus class
//    if (idof <=0 || idof > 2) error("RegularBC3d::wrong input for degree of freedom");
    if (idof <=0) error("RegularBC2d::wrong input for degree of freedom");
}


