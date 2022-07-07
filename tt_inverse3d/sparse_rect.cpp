/*
 * sparse_rect.cc - sparse rectangular matrix implementation
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#include <iostream>
#include <fstream>

#include <omp.h>

#include "sparse_rect.hpp"
#include "error.hpp"

#include "config.hpp"


namespace {
    template< typename T>
    bool is_nan(T const& v) {
        return v != v;
    }
    bool
    has_nan(std::map<int,double> const& m) {
        for(auto const& p : m) {
            if (is_nan(p.second)) {
                return true;
            }
        }
        return false;
    }
    
    bool
    has_nan(std::vector< std::map<int,double>> const& m) {
        for(auto const& l : m) {
            if (has_nan(l)) {
                return true;
            }
        }
        return false;
    }
}

void
print_sparse_line(std::ostream& out, std::map<int,double> const& m) {
    out << '{';
    for( auto const& v : m) {
        out << '[' << v.first << ':' << v.second << "] ";
    }
    out << "}\n";
}

void
print_sparse_matrix(std::ostream& out, Array1d<std::map<int,double>> const& m) {
    for(auto const& l : m) {
        print_sparse_line(out, l);
    }
}

void
dump_sparse_line(std::map<int,double> const& l) {
    print_sparse_line(std::cout, l);
}

void
dump_sparse_matrix(Array1d<std::map<int,double>> const& m) {
    print_sparse_matrix(std::cout, m);
}

SparseRectangularData::SparseRectangularData(const sparse_matrix& A, int n)
    : ncol(n),
      nrow(A.size())
{
    int ntotal=0;
    for (auto const& a : A.stl()) {
	ntotal += a.size();
    }
    
    coef.resize(ntotal);
    index_j.resize(ntotal);
    first_k.resize(nrow+1);

    {
        int k=1;
        for (int i=0; i < nrow; ++i){
            first_k[i] = k;
            for (auto const& p : A[i]) {
                index_j[k-1] = p.first;
                coef[k-1] = p.second;
                assert(!is_nan(coef[k-1]));
                k++;
            }
        }
        first_k.back() = k;
    }
    
    if (ncol==0){
	// scan index_j to get the largest colmun index (= ncol)
	for (int i=2; i<=nrow+1; i++){
	    int last_k = first_k[i-1]-1;
	    int jmax   = index_j[last_k-1];
	    if (jmax>ncol) {
                ncol = jmax;
            }
	}
    }
}

SparseRectangularData::SparseRectangularData(sparse_matrix const& A,
                                             std::function<bool(int)> const& rfilter,
                                             std::function<int(int)>  const& cfilter,
                                             int n)
    : ncol(n),
      nrow(0)
{
    typedef std::map<int,double> row;

    int ntotal=0;
    for (int i=1; i<= A.size(); i++){
	if (rfilter(i)){
	    ++nrow;
            row const& r = A(i);
	    for (row::const_iterator p = r.begin(); p != r.end(); ++p) {
		int j = p->first;
		if (bool(cfilter(j))) {
                    ++ntotal;
                }
	    }
	}
    }
    
    coef.resize(ntotal);
    index_j.resize(ntotal);
    first_k.resize(nrow+1);
    
    int k=1, ii=1;
    for (int i=1; i<=A.size(); i++){
	if (rfilter(i)){
	    first_k[ii++ - 1] = k;
            row const& r = A(i);
	    for (row::const_iterator p = r.begin(); p != r.end(); ++p){
		int j = p->first;
		if (bool(cfilter(j))){
		    index_j[k-1] = cfilter(j);
		    coef[k-1] = p->second;
                    assert(!is_nan(coef[k-1]));
		    ++k;
		}
	    }
	}
    }
    first_k[nrow] = k;
    
    if (ncol==0){
	// scan index_j to get the largest column index (= ncol)
	for (int i=2; i<=nrow+1; i++){
	    int last_k = first_k[i-1]-1;
	    int jmax = index_j[last_k-1];
	    if (jmax>ncol) {
                ncol = jmax;
            }
	}
    }
}

SparseRectangularData::SparseRectangularData(const sparse_matrix& A,
                                             const std::vector<int>& row_index,
                                             const std::vector<int>& col_index,
                                             int n)
    : SparseRectangularData(A,
                            explicit_filter<bool>(row_index),
                            explicit_filter<int>(col_index),
                            n) {}

SparseRectangularData::SparseRectangularData(const Array1d<const sparse_matrix*>& As,
                                             const Array1d<SparseMatAux>& spec, int n)
    : ncol(n),
      nrow(0)
{
    typedef std::map<int,double>::const_iterator mapBrowser;

    assert(As.size() == spec.size());
    int ntotal=0;
    for (int iA=1; iA<=As.size(); iA++){
	for (int i=1; i<=As(iA)->size(); i++){
	    if (spec(iA).row(i)) {
		++nrow;
		for (mapBrowser p=(*As(iA))(i).begin(); p!=(*As(iA))(i).end(); p++){
		    int j = p->first;
		    if (bool(spec(iA).col(j))) {
                        ++ntotal;
                    }
		}
	    }
	}
    }

    coef.resize(ntotal);
    index_j.resize(ntotal);
    first_k.resize(nrow+1);
    int k=1, ii=1;
    for (int iA=1; iA<=As.size(); iA++){
        sparse_matrix const& as = *(As(iA));
	double factor=spec(iA).scale_factor;
        cerr << "sparse_rect " << iA << " " << factor << '\n';
	for (int i=1; i <= as.size(); i++){
	    if (spec(iA).row(i)) {
		first_k[ii++ - 1] = k;
                std::map<int,double> const& line = as(i);
		for (auto const& p : line ) {
		    int j = p.first;
		    if (spec(iA).col(j)>0) {
			index_j[k-1] = spec(iA).col(j);
			coef[k-1] = factor*p.second;
                        assert(!is_nan(coef[k-1]));
			k++;
		    }
		}
	    }
	}
    }
    first_k[nrow] = k;

    if (ncol==0){
	// scan index_j to get the largest colmun index (= ncol)
	for (int i=2; i<=nrow+1; i++){
	    int last_k = first_k[i-1]-1;
	    int jmax = index_j[last_k-1];
	    if (jmax>ncol) ncol = jmax;
	}
    }
}

// compute Ax = b
void
SparseRectangular::Ax(std::vector<double> const& x, std::vector<double>& b, bool init) const
{
    assert(x.size() == ncol);
    assert(b.empty() || b.size() == nrow);
    
    if (b.empty()) {
        b.resize(nrow, 0.0);
    } else {
        if (init) {
            std::fill(b.begin(), b.end(), 0);
        }
    }

    std::size_t const n = b.size();
    for (int i = 0; i < n; ++i) {
        int begin = first_k[i] - 1;
        int end   = first_k[i+1] - 1;
        double v = 0;
	for (int k = begin; k < end; ++k){
	    v += coef[k]*x[index_j[k]-1];
	}
        b[i] += v;
    }
}

aligned_vector<double>
SparseRectangular::omp_mul(aligned_vector<double> const& b) const {
    assert(b.size() == ncol);
    aligned_vector<double> prod(nrow, 0);
    int    const* restrict first = first_k.data();
    double const* restrict coef  = this->coef.data();
    int    const* restrict js    = index_j.data();
    double* sol = prod.data();
#pragma omp parallel for 
    for (int i = 0; i < nrow; ++i) {
        int begin = first[i] - 1;
        int end   = first[i+1] - 1;
        if (begin != end) {
            double v = 0;
            for (int k = begin; k < end; ++k){
                v += coef[k]*b[js[k]-1];
            }
            sol[i] = v;
        }
    }
    return prod;
}

namespace {
    void
    omp_sum(aligned_vector<double>& a, aligned_vector<double> const& b) {
        assert(a.size() == b.size());
        double*       aarr = a.data();
        double const* barr = b.data();
        int const nb = a.size();
#pragma omp parallel for 
        for(int i=0; i < nb; ++i) {
            aarr[i] += barr[i];
        }
    }
}


// Compute Ax=b
void
SparseRectangular::omp_mul(aligned_vector<double> const& b, aligned_vector<double>& x, bool init) const
{
    aligned_vector<double> prod = omp_mul(b);
    if (init) {
        x = std::move(prod);
    } else {
        omp_sum(x,prod);
    }
}

// calculate At*x = b
void
SparseRectangular::Atx(std::vector<double> const& x, std::vector<double>& b, bool init) const
{
    assert(x.size() == nrow);
    assert(b.empty() || b.size() == ncol);
    
    if (b.empty()) {
        b.resize(ncol, 0.0);
    } else {
        if (init){
            std::fill(b.begin(), b.end(), 0);
        }
    }
    
    auto const& firsts = first_k;
    assert(firsts.size() >= 2);

    double const* restrict coef = this->coef.data();
    double*       restrict rhs = &(b.front());
    for (int c = 0; c < nrow; ++c) {
        int begin = firsts[c];
        int end   = firsts[c+1];
        for (int k = begin; k < end; ++k){
            int j = index_j[k-1] - 1;
            rhs[j] += coef[k-1]*x[c];
            assert(!is_nan(rhs[j]));
        }
    }
}

aligned_vector<double>
SparseRectangular::omp_mul_transposed(aligned_vector<double> const& x) const {
    assert(x.size() == nrow);
    aligned_vector<double> prod(ncol, 0);
    assert(first_k.size() >= 2);
    
    int    const* restrict firsts = first_k.data();
    int    const* restrict js     = index_j.data();
    double const* restrict coef   = this->coef.data();
    double const* restrict vec    = x.data();

    std::vector<aligned_vector<double>> locbs;
#pragma omp parallel shared(locbs)
    {
        int threadid = omp_get_thread_num();
        int nbthread = omp_get_num_threads();
#pragma omp single
        {
            locbs.resize(nbthread);
        }
#pragma omp barrier
        auto& locb = locbs[threadid];
        locb.resize(ncol, 0);
        double* restrict partial = locb.data();
        int chunk = nrow/nbthread;
#pragma omp for schedule(dynamic, chunk)
        for (int c = 0; c < nrow; c++ ) {
            int begin = firsts[c] - 1;
            int end   = firsts[c+1] - 1;
            if (begin != end) {
                double v  = vec[c];
                for (int k = begin; k < end; ++k){
                    int j = js[k] - 1;
                    partial[j] += coef[k]*v;
                }
            }
        }
    }
    int size = locbs.size();
    while (size != 1) {
        int half = (size+1)/2;
        for (int h1 = 0; h1 < half; ++h1) {
            int h2 = h1 + half;
            if (h2 < size) {
                double *      restrict d1 = locbs[h1].data();
                double const* restrict d2 = locbs[h2].data();
#pragma omp parallel for
                for (int i = 0; i < ncol; ++i) {
                    d1[i] += d2[i];
                }
            }
        }
        size = half;
        
    }
    prod.swap(locbs[0]);
    return prod;
}

void
SparseRectangular::omp_mul_transposed(aligned_vector<double> const& x, aligned_vector<double>& b, bool clear) const
{
    aligned_vector<double> prod = omp_mul_transposed(x);
    if (clear) {
        b = std::move(prod);
    } else {
        omp_sum(b, prod);
    }
}

// dump out the contents (for debugging)
void SparseRectangular::dump(const char* fn) const
{
    std::unique_ptr<std::ostream> ofile;
    std::ostream* os_ptr;
    if (fn) {
        ofile.reset(new std::ofstream(fn));
        os_ptr = ofile.get();
    } else {
        os_ptr = &std::cout;
    }
    
    std::ostream& os(*os_ptr);

    os << "SparseRectangular dump\n";

    os << "index_j=\n";
    for (int i=1; i<=index_j.size(); i++){
	os << index_j[i-1] << " ";
	//if (i%10==0) os << '\n';
    }
    os << '\n';

    os << "first_k=\n";
    for (int i=1; i<=first_k.size(); i++){
	os << first_k[i-1] << " ";
	//if (i%10==0) os << '\n';
    }
    os << '\n';

    os << "val=\n";
    for (int i=1; i<=coef.size(); i++){
	os << coef[i-1] << " ";
    }
    os << "\n\n\n";
    int k = 1;
    for (int i=1; i <= (first_k.size()-1); i++){
        for( int j = first_k[i-1]; j < first_k[i]; ++j) {
            os << index_j[j-1] << ": " << coef[k++-1] << ", ";
        }
        os << '\n';
	//if (i%10==0) os << '\n';
    }
    os << '\n';
}

