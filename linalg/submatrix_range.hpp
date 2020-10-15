#ifndef SUBMATRIX_HPP_
#define SUBMATRIX_HPP_

/* This is only an abstract description of a submatrix range. It can be
 * used on any matrix type that has methods .nrows() and .ncols()
 */
#include "macros.h"

struct submatrix_range {
    unsigned int i0=0,j0=0;
    unsigned int i1=0,j1=0;
    unsigned int nrows() const { return i1-i0; }
    unsigned int ncols() const { return j1-j0; }
    size_t size() const { return (size_t) nrows() * (size_t) ncols(); }
    submatrix_range & range() { return *this; }
    submatrix_range const & range() const { return *this; }
    submatrix_range() = default;
    submatrix_range(unsigned int i0, unsigned int j0, unsigned int ni, unsigned int nj) : i0(i0), j0(j0), i1(i0+ni), j1(j0+nj) {}
    template<typename T>
    submatrix_range(T const & M) : i0(0), j0(0), i1(M.nrows()), j1(M.ncols()) {}
    template<typename T>
    inline bool valid(T const & a) const {
        return i0 <= i1 && i1 <= a.nrows() && j0 <= j1 && j1 <= a.ncols();
    }
    submatrix_range operator*(submatrix_range const & a) const {
        ASSERT_ALWAYS(ncols() == a.nrows());
        return submatrix_range(i0, nrows(), a.j0, a.ncols());
    }
};


#endif	/* SUBMATRIX_HPP_ */
