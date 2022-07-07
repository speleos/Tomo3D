/*
 * geom3d.h based on geom.h
 *
 * Adria Melendez and Jun Korenaga
 * Spring 2011
 */

#ifndef _TOMO_GEOM3D_H_
#define _TOMO_GEOM3D_H_

#include <iosfwd>
#include <vector>
#include <cmath>

// Point3d  //modified
class Point3d {
public:
    Point3d(){}
    Point3d(double, double, double);

    double x() const { return v[0]; }
    double y() const { return v[1]; }
    double z() const { return v[2]; }
    double& x() { return v[0]; }
    double& y() { return v[1]; }
    double& z() { return v[2]; }
    void   x(double d) { v[0] = d; }
    void   y(double d) { v[1] = d; }
    void   z(double d) { v[2] = d; }
    void   set(double dx, double dy, double dz) { v[0] = dx; v[1] = dy; v[2] = dz; }
    double distance(const Point3d& p) const;
    double norm() const { return sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); }
    double inner_product(const Point3d& p) const;
    Point3d outer_product(const Point3d& p) const;

    // unary operators
    Point3d operator-() const;
    
    // binary operators
    Point3d& operator+=(const Point3d&);
    Point3d& operator-=(const Point3d&);
    Point3d& operator*=(double);
    Point3d& operator/=(double);
    
private:
    double v[3];
};

inline double Point3d::distance(const Point3d& p) const
{
    double dx = p.v[0]-v[0];
    double dy = p.v[1]-v[1];
    double dz = p.v[2]-v[2];
    return sqrt(dx*dx+dy*dy+dz*dz);
}

inline
Point3d Point3d::operator-() const
{
    Point3d neg;

    neg.v[0] = -v[0];
    neg.v[1] = -v[1];
    neg.v[2] = -v[2];

    return neg;
}

inline
Point3d& Point3d::operator+=(const Point3d& a)
{
    v[0] += a.v[0];
    v[1] += a.v[1];
    v[2] += a.v[2];
    return *this;
}

inline
Point3d& Point3d::operator-=(const Point3d& a)
{
    v[0] -= a.v[0];
    v[1] -= a.v[1];
    v[2] -= a.v[2];
    return *this;
}

inline
Point3d& Point3d::operator*=(double a)
{
    v[0] *= a;
    v[1] *= a;
    v[2] *= a;
    return *this;
}

inline
Point3d& Point3d::operator/=(double a)
{
    v[0] /= a;
    v[1] /= a;
    v[2] /= a;
    return *this;
}

 // nonmember functions
inline
Point3d operator+(const Point3d& a, const Point3d& b)
{
    Point3d c=a;
    return c+=b;
}

inline
Point3d operator-(const Point3d& a, const Point3d& b)
{
    Point3d c=a;
    return c-=b;
}

inline
Point3d operator*(const Point3d& a, double val)
{
    Point3d c=a;
    return c*=val;
}

inline
Point3d operator*(double val, const Point3d& a)
{
    return operator*(a,val);
}

inline
Point3d operator/(const Point3d& a, double val)
{
    Point3d c=a;
    return c/=val;
}

std::ostream& operator<<(std::ostream& out, Point3d const& p);
std::ostream& operator<<(std::ostream& out, std::vector<Point3d> const& p);

// RegularDomain3d
class RegularDomain3d {
  public:
    RegularDomain3d(){ }
    RegularDomain3d(double, double, double, double); // 3-D 

    bool inXRange(double) const;
    bool inYRange(double) const;
    bool inDomain(const Point3d&) const;
    
    double x_min() const { return vmin[0]; }
    double x_max() const { return vmax[0]; }
    double y_min() const { return vmin[1]; }
    double y_max() const { return vmax[1]; }

    void x_min(double d) { vmin[0] = d; }
    void x_max(double d) { vmax[0] = d; }
    void y_min(double d) { vmin[1] = d; }
    void y_max(double d) { vmax[1] = d; }

  private:
    double vmin[2];
    double vmax[2];
};

class RegularBC3d : public RegularDomain3d {
  public:
    RegularBC3d(int, double, double, double, double, double); // 3-D
	
    int    iDegOfFreedom() const { return idof; }
    double val() const           { return val_; }
    void   set_iDegOfFreedom(int i) { idof = i; }
    void   set_val(double d)        { val_ = d; }

  private:
    int idof;			// index for deg of freedom
    double val_;
};

#endif /* _TOMO_GEOM3D_H_ */
