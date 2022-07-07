/*
 * array.h
 */

#ifndef _JK_NEW_ARRAY_H_
#define _JK_NEW_ARRAY_H_

#include <cassert>
#include <vector>
#include <algorithm>
#include <cstddef>
#include <iostream>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

#include "error.hpp"

//using namespace std;

template<class T>
class Array1d : public std::vector<T> {
public:
    typedef std::vector<T> super;

    Array1d(size_t n, T val) : vector<T>(n, val) {}
    Array1d(size_t n = 0) : vector<T>(n) {}
    Array1d(const T*, size_t);	// array initialization
    
    std::vector<T>& stl() { return *this; }
    std::vector<T> const& stl() const { return *this; }
    
    T& operator()(size_t i) {
        assert(i>0); assert(i<=size());
        return this->data()[i-1];
    }
    const T& operator()(size_t i) const {
        return const_cast<Array1d<T>&>(*this)(i);
    }
    size_t size() const { return vector<T>::size(); }
    T* begin() { return &(vector<T>::front()); }
    T* end() { return (&(vector<T>::back()))+1; }
    const T* begin() const { return &(vector<T>::front()); }
    const T* end() const { return (&(vector<T>::back()))+1; }


    void resize(size_t n, const T& val = T()) { vector<T>::resize(n, val); }
    void push_back (const T& x) { vector<T>::push_back(x); }

    // for Numerical recipes' vector
    T*       toRecipe() { return begin()-1; }
    
    // unary operators
    Array1d<T> operator-();
    
    // binary operators
    Array1d<T>& operator=(const T& val);
    Array1d<T>& operator+=(const Array1d<T>&);
    Array1d<T>& operator+=(const T& val);
    Array1d<T>& operator-=(const Array1d<T>&);
    Array1d<T>& operator-=(const T& val);
    Array1d<T>& operator*=(const T&);
    Array1d<T>& operator/=(const T&);
};

template<class T>
inline
Array1d<T>::Array1d(const T* init, size_t n)
: vector<T>(init, init+n)
{}
    
template<class T>
inline
Array1d<T> Array1d<T>::operator-()
{
    Array1d<T> negative(size());
    std::transform(begin(), end(), negative.begin(), std::negate<T>());
    return negative;
}

template<class T>
inline
Array1d<T>& Array1d<T>::operator=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest = val;
    return *this;
}

template<class T>
inline
Array1d<T>& Array1d<T>::operator+=(const Array1d<T>& a)
{
    if (size() != a.size()) error("Array1d::operator+= size mismatch");
    T* dest = end();
    const T* src = a.end();
    while (dest > begin()) *--dest += *--src;
    return *this;
}

template<class T>
inline
Array1d<T>& Array1d<T>::operator+=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest += val;
    return *this;
}

template<class T>
inline
Array1d<T>& Array1d<T>::operator-=(const Array1d<T>& a)
{
    if (size() != a.size()) error("Array1d::operator-= size mismatch");
    T* dest = end();
    const T* src = a.end();
    while (dest > begin()) *--dest -= *--src;
    return *this;
}

template<class T>
inline
Array1d<T>& Array1d<T>::operator-=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest -= val;
    return *this;
}

template<class T>
inline
Array1d<T>& Array1d<T>::operator*=(const T& a)
{
    T* dest = end();
    while (dest > begin()) *--dest *= a;
    return *this;
}
    
template<class T>
inline
Array1d<T>& Array1d<T>::operator/=(const T& a)
{
    T* dest = end();
    while (dest > begin()) *--dest /= a;
    return *this;
}


// nonmember functions

template<class T>
ostream& operator<<(ostream& s, const Array1d<T>& a)
{
    for (int i=1; i<a.size(); i++){
	s << a(i) << ", ";
    }
    s << a(a.size()) << '\n';
    return s;
}

// Array2d
template<class T>
class Array2d : public vector<T> {
public:
    Array2d() : vector<T>(0) {}
    Array2d(size_t n1, size_t n2)
	: vector<T>(n1*n2), nrow(n1), ncol(n2) {}
    ~Array2d();
    
    T& operator()(size_t i, size_t j){
        assert(i>0 && i <= nRow());
        assert(j>0 && j <= nCol());
        return vector<T>::operator[](offset(i,j)); 
    }
    T operator()(size_t i, size_t j) const{
        assert(i>0 && i <= nRow());
        assert(j>0 && j <= nCol());
        return vector<T>::operator[](offset(i,j)); 
    }
    T* begin() { return &(vector<T>::front()); }
    T* end() { return (&(vector<T>::back())+1); }
    const T* begin() const { return &(vector<T>::front()); }
    const T* end() const { return (&(vector<T>::back())+1); }
    
    void resize(size_t n1, size_t n2, const T& val = T())
	{ vector<T>::resize(n1*n2, val); nrow=n1, ncol=n2; }

    size_t nCol() const { return ncol; }
    size_t nRow() const { return nrow; }

    T** toRecipe();

    // unary operators
    Array2d<T> operator-();
    
    // binary operators
    Array2d<T>& operator=(const T& val);
    Array2d<T>& operator+=(const Array2d<T>&);
    Array2d<T>& operator+=(const T& val);
    Array2d<T>& operator-=(const Array2d<T>&);
    Array2d<T>& operator-=(const T& val);
    Array2d<T>& operator*=(const T&);
    Array2d<T>& operator/=(const T&);
    
    std::vector<T>& stl() { return *this; }
    std::vector<T> const& stl() const { return *this; }

private:
    size_t size() const { return vector<T>::size(); }

    size_t offset(size_t i, size_t j) const;

    size_t nrow;
    size_t ncol;

    vector<T*> datapp;			// for toRecipe();
};

/// \brief This base class encapsulate the 2d representation of a 
/// 2d NX x NY plane stored in a flatt 1d array
template<typename T>
class plane {
public:
    plane(std::vector<T> const& a, int nx, int ny) 
        : my_plane(a), my_nx(nx), my_ny(ny) {
        assert (a.size() == nx*ny);
        assert (nx>0);
        assert (ny>0);
    }
    std::vector<T> const& data() const { return my_plane; }
    int                   nx()    const { return my_nx; }
    int                   ny()    const { return my_ny; }
    
    bool in_range(int i, int j) const {
        return (i >= 0 && j >= 0
                && i < my_nx && j < my_ny );
    }
    
    T const& at(int i, int j) const {
        assert(in_range(i,j));
        return my_plane[i*my_ny + j];
    }
    T const& operator()(int i, int j) const {
        return at(i,j);
    }

    T const& at(int i, int j, T const& def) const {
        return in_range(i,j) ? my_plane[i*my_ny + j] : def;
    }
    T const& operator()(int i, int j, T const& def) const {
        return at(i,j,def);
    }
private:
    std::vector<T> const& my_plane;
    int                   my_nx;
    int                   my_ny;
};


/// \brief This class allow the traversal of a single slice of a 2d NX X NY plane 
/// assuming the first dimension is fixed.
/// 
/// The wrapped array is a flat vector representation of a nx x ny plane
template<typename T>
class frozenx : public plane<T> {
public:
  typedef plane<T> super;
    frozenx(std::vector<T> const& a, int nx, int ny, int x)
        : plane<T>(a,nx,ny), my_index(x) {
        assert (this->in_range(x,0));
    }
    int                   index() const { return my_index; }
    T const&              operator[](int i) const { return super::at(index(),i); }
    
private:
    int                   my_index;
};

/// \brief This class allow the traversal of a single slice of a 2d NX X NY plane 
/// assuming the second dimension is fixed.
/// 
/// The wrapped array is a flat vector representation of a nx x ny plane
template<typename T>
class frozeny : public plane<T> {
public:
    typedef plane<T> super;
    frozeny(std::vector<T> const& a, int nx, int ny, int y)
        : plane<T>(a,nx,ny), my_index(y) {
        assert (this->in_range(0,y));
    }
    int                   index() const { return my_index; }
    T const&              operator[](int i) const { return super::at(i,index());}
    
private:
    int                   my_index;
};

template<class T> 
inline size_t Array2d<T>::offset(size_t i, size_t j) const
{
//    if ((i-1)*ncol+j-1 >=size()){
//	cerr << "Array2d OoR " << i << " " << j << " " << nrow << " " << ncol << '\n';
//    }
    return (i-1)*ncol+j-1;
}

template<class T>
Array2d<T>::~Array2d()
{
}

template<class T>
T** Array2d<T>::toRecipe()
{
    datapp.resize(nrow+1);
    datapp[1] = begin()-1;
    for (int i=2; i<=nrow; i++) datapp[i] = datapp[i-1]+ncol;
    return &(datapp.front());
}

template<class T>
inline
Array2d<T>& Array2d<T>::operator=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest = val;
    return *this;
}

template<class T>
inline
Array2d<T> Array2d<T>::operator-()
{
    Array2d<T> negative(size());

    T* dest = negative.end();
    const T* src = end();
    while (dest > negative.begin()) *--dest = (-1)*(*--src);
    return negative;
}

template<class T>
inline
Array2d<T>& Array2d<T>::operator+=(const Array2d<T>& a)
{
    if (size() != a.size()) error("Array2d::operator+= size mismatch");
    T* dest = end();
    const T* src = a.end();
    while (dest > begin()) *--dest += *--src;
    return *this;
}

template<class T>
inline
Array2d<T>& Array2d<T>::operator+=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest += val;
    return *this;
}

template<class T>
inline
Array2d<T>& Array2d<T>::operator-=(const Array2d<T>& a)
{
    if (size() != a.size()) error("Array2d::operator-= size mismatch");
    T* dest = end();
    const T* src = a.end();
    while (dest > begin()) *--dest -= *--src;
    return *this;
}

template<class T>
inline
Array2d<T>& Array2d<T>::operator-=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest -= val;
    return *this;
}

template<class T>
inline
Array2d<T>& Array2d<T>::operator*=(const T& a)
{
    T* dest = end();
    while (dest > begin()) *--dest *= a;
    return *this;
}
    
template<class T>
inline
Array2d<T>& Array2d<T>::operator/=(const T& a)
{
    T* dest = end();
    while (dest > begin()) *--dest /= a;
    return *this;
}

// nonmember functions
template<class T>
inline
Array2d<T> operator+(const Array2d<T>& a, const Array2d<T>& b)
{
    Array2d<T> r=a;
    return r+=b;
}

template<class T>
inline
Array2d<T> operator-(const Array2d<T>& a, const Array2d<T>& b)
{
    Array2d<T> r=a;
    return r-=b;
}

template<class T>
inline
Array2d<T> operator*(const Array2d<T>& a, const T& val)
{
    Array2d<T> r=a;
    return r*=val;
}

template<class T>
inline
Array2d<T> operator*(const T& val, const Array2d<T>& a)
{
    return operator*(a,val);
}

template<class T>
inline
Array2d<T> operator/(const Array2d<T>& a, const T& val)
{
    Array2d<T> r=a;
    return r/=val;
}

template<class T>
ostream& operator<<(ostream& s, const Array2d<T>& a)
{
    for (int i=1; i<=a.nRow(); i++){
	for (int j=1; j<a.nCol(); j++){
	    s << a(i,j) << ", ";
	}
	s << a(i,a.nCol()) << '\n';
    }
    return s;
}

// Array3d
template<class T>
class Array3d : public vector<T> {
public:
    Array3d() : vector<T>(0) {}
    Array3d(size_t n1, size_t n2, size_t n3)
	: vector<T>(n1*n2*n3), nn1(n1), nn2(n2), nn3(n3) {}
    ~Array3d();
    
    T& operator()(size_t i, size_t j, size_t k) {
        assert(i>0 && i <= nSize1());
        assert(j>0 && j <= nSize2());
        assert(k>0 && k <= nSize3());
        return vector<T>::operator[](offset(i,j,k));
    }
    T operator()(size_t i, size_t j, size_t k) const {
        assert(i>0 && i <= nSize1());
        assert(j>0 && j <= nSize2());
        assert(k>0 && k <= nSize3());
        return vector<T>::operator[](offset(i,j,k));
    }
    T* begin() { return &(vector<T>::front()); }
    T* end() { return (&(vector<T>::back())+1); }
    const T* begin() const { return &(vector<T>::front()); }
    const T* end() const { return (&(vector<T>::back())+1); }
    
    void resize(size_t n1, size_t n2, size_t n3, const T& val = T()) {
        vector<T>::resize(n1*n2*n3, val);
        nn1=n1, nn2=n2; nn3=n3;
    }

    size_t nSize1() const { return nn1; }
    size_t nSize2() const { return nn2; }
    size_t nSize3() const { return nn3; }

    // unary operators
    Array3d<T> operator-();
    
    // binary operators
    Array3d<T>& operator=(const T& val);
    Array3d<T>& operator+=(const Array3d<T>&);
    Array3d<T>& operator+=(const T& val);
    Array3d<T>& operator-=(const Array3d<T>&);
    Array3d<T>& operator-=(const T& val);
    Array3d<T>& operator*=(const T&);
    Array3d<T>& operator/=(const T&);
    
private:
    size_t size() const { return vector<T>::size(); }
    size_t offset(size_t i, size_t j, size_t k) const;

    size_t nn1, nn2, nn3;
};

template<class T> 
inline size_t Array3d<T>::offset(size_t i, size_t j, size_t k) const
{
    return ((i-1)*nn2+(j-1))*nn3+(k-1);
}

template<class T>
Array3d<T>::~Array3d()
{
}

template<class T>
inline
Array3d<T>& Array3d<T>::operator=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest = val;
    return *this;
}

template<class T>
inline
Array3d<T> Array3d<T>::operator-()
{
    Array3d<T> negative(size());

    T* dest = negative.end();
    const T* src = end();
    while (dest > negative.begin()) *--dest = (-1)*(*--src);
    return negative;
}

template<class T>
inline
Array3d<T>& Array3d<T>::operator+=(const Array3d<T>& a)
{
    if (size() != a.size()) error("Array3d::operator+= size mismatch");
    T* dest = end();
    const T* src = a.end();
    while (dest > begin()) *--dest += *--src;
    return *this;
}

template<class T>
inline
Array3d<T>& Array3d<T>::operator+=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest += val;
    return *this;
}

template<class T>
inline
Array3d<T>& Array3d<T>::operator-=(const Array3d<T>& a)
{
    if (size() != a.size()) error("Array3d::operator-= size mismatch");
    T* dest = end();
    const T* src = a.end();
    while (dest > begin()) *--dest -= *--src;
    return *this;
}

template<class T>
inline
Array3d<T>& Array3d<T>::operator-=(const T& val)
{
    T* dest = end();
    while (dest > begin()) *--dest -= val;
    return *this;
}

template<class T>
inline
Array3d<T>& Array3d<T>::operator*=(const T& a)
{
    T* dest = end();
    while (dest > begin()) *--dest *= a;
    return *this;
}
    
template<class T>
inline
Array3d<T>& Array3d<T>::operator/=(const T& a)
{
    T* dest = end();
    while (dest > begin()) *--dest /= a;
    return *this;
}

// nonmember functions
template<class T>
inline
Array3d<T> operator+(const Array3d<T>& a, const Array3d<T>& b)
{
    Array3d<T> r=a;
    return r+=b;
}

template<class T>
inline
Array3d<T> operator-(const Array3d<T>& a, const Array3d<T>& b)
{
    Array3d<T> r=a;
    return r-=b;
}

template<class T>
inline
Array3d<T> operator*(const Array3d<T>& a, const T& val)
{
    Array3d<T> r=a;
    return r*=val;
}

template<class T>
inline
Array3d<T> operator*(const T& val, const Array3d<T>& a)
{
    return operator*(a,val);
}

template<class T>
inline
Array3d<T> operator/(const Array3d<T>& a, const T& val)
{
    Array3d<T> r=a;
    return r/=val;
}

template<class T>
ostream& operator<<(ostream& s, const Array3d<T>& a)
{
    for (int i=1; i<=a.nSize1(); i++){
	for (int j=1; j<=a.nSize2(); j++){
	    for (int k=1; k<=a.nSize3(); k++){
		s << a(i,j,k) << ", ";
	    }
	    s << '\n';
	}
	s << "\n\n";
    }
    return s;
}

#endif /* _JK_NEW_ARRAY_H_ */
