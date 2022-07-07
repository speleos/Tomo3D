/*
 * sparse_rect.h - sparse rectangular matrix interface
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#ifndef _TOMO_SPARSE_RECT_H_
#define _TOMO_SPARSE_RECT_H_

#include <map>
#include <limits>
#include <functional>
#include <memory>   //ALOUREIRO 11/02/2022


#include "boost/optional.hpp"
#include "boost/serialization/access.hpp"
#include "boost/serialization/vector.hpp"
#include "boost/align/aligned_allocator.hpp"

#include "array.hpp" // from mconv

class SparseRectangularData;

template<typename T> using aligned_vector = std::vector<T, boost::alignment::aligned_allocator<T>>;

template<typename T>
class explicit_filter {
public:
    explicit_filter(std::vector<int> const& a) : my_answers(a) {}
    explicit_filter() = default;

    T operator()(int i) const {
        return T(my_answers->at(i-1));
    }
private:
    boost::optional<std::vector<int> const&> my_answers;
};

class SparseMatAux {
public:
    SparseMatAux(double d,
                 std::function<bool(int)> const& r,
                 std::function<int(int)>  const& c)
	: scale_factor(d), 
          my_rfilter(r),
          my_cfilter(c) {}

    SparseMatAux(double d,
                 std::vector<int> const& r,
                 std::vector<int> const& c)
	: scale_factor(d),
          my_rfilter(explicit_filter<bool>(r)), 
          my_cfilter(explicit_filter<int>(c)) {}
    
    SparseMatAux() = default;
    
    SparseMatAux(SparseMatAux const&)   = default;
    SparseMatAux& operator=(SparseMatAux const&) = default;
    ~SparseMatAux() = default;

    double scale_factor;
    
private:
    friend class SparseRectangularData;
    bool row(int i) const { return my_rfilter(i); }
    int  col(int i) const { return my_cfilter(i); }
    
private:
    std::function<bool(int)> my_rfilter;
    std::function<int(int)>  my_cfilter;
};

typedef std::map<int,double> sparse_line;

class sparse_matrix 
    : public Array1d<sparse_line> {
public:
    typedef Array1d<sparse_line> super;
    using super::super;
};

bool has_nan(map<int,double> const& m);
bool has_nan(Array1d< map<int,double>> const& m);

void print_sparse_line(std::ostream& out, map<int,double> const& m);
void print_sparse_matrix(std::ostream& out, Array1d<map<int,double>> const& m);

void dump_sparse_line(map<int,double> const& m);
void dump_sparse_matrix(Array1d<map<int,double>> const& m);

class SparseRectangularData {
public:
    SparseRectangularData(){};
    SparseRectangularData(const sparse_matrix& A, int ncol);
    SparseRectangularData(const sparse_matrix& A,
                          const std::vector<int>& r, const std::vector<int>& c, 
                          int ncol=0);
    SparseRectangularData(const sparse_matrix& A,
                          std::function<bool(int)> const& r, 
                          std::function<int(int)>  const& c, 
                          int ncol=0);
    SparseRectangularData(const Array1d<const sparse_matrix*>& As,
                          const Array1d<SparseMatAux>& spec, int ncol);
    ~SparseRectangularData() {};

private:
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar << nrow;
        ar << ncol;
        ar << coef;
        ar << index_j;
        ar << first_k;
    }

    template<class Archive, typename T>
    void
    load_from_v0(Archive& ar, aligned_vector<T>& v) {
        std::vector<T> u;
        ar >> u;
        v.resize(u.size());
        std::copy(u.begin(), u.end(), v.begin());
    }
    
    template<class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        ar & nrow;
        ar & ncol;
        if (version == 0) {
            load_from_v0(ar, coef);
            load_from_v0(ar, index_j);
            load_from_v0(ar, first_k);
        } else {
            ar >> coef;
            ar >> index_j;
            ar >> first_k;
        }
    }
    
    BOOST_SERIALIZATION_SPLIT_MEMBER();
    
public:
    int nrow, ncol;
    aligned_vector<double> coef;
    aligned_vector<int>    index_j;
    aligned_vector<int>    first_k;
};

BOOST_CLASS_VERSION(SparseRectangularData, 1);

class SparseRectangular : private SparseRectangularData {
public:
    typedef SparseRectangularData super;
    using super::super;
    ~SparseRectangular() {};
    
    int nRow() const { return nrow; }
    int nCol() const { return ncol; }
    void Ax (std::vector<double> const&, std::vector<double>&, bool init = true) const;
    void omp_mul(aligned_vector<double> const& b, aligned_vector<double>& x, bool init = true) const;
    aligned_vector<double> omp_mul(aligned_vector<double> const& b) const;
    void omp_mul_transposed(aligned_vector<double> const& b, aligned_vector<double>& x, bool init = true) const;
    aligned_vector<double> omp_mul_transposed(aligned_vector<double> const& b) const;
    void Atx(std::vector<double> const&, std::vector<double>&, bool init = true) const;
    void dump(const char* fname = 0) const;
    
    SparseRectangularData const& raw() const { return *this; }

private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        // serialize base class information
        ar & boost::serialization::base_object<SparseRectangularData>(*this);
    }
};

#endif /* _TOMO_SPARSE_RECT_H_ */
