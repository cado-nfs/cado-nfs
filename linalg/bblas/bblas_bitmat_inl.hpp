#ifndef BBLAS_BITMAT_INL_HPP_
#define BBLAS_BITMAT_INL_HPP_

#include "bblas.hpp"
#include "bblas_bitmat.hpp"
#include "gmp_aux.h"

namespace bblas_bitmat_details {
    template<typename matrix>
    void bitmat_ops<matrix>::fill_random(matrix & w, gmp_randstate_t rstate)
    {
        memfill_random(w.data(), sizeof(matrix), rstate);
    }

    template<typename matrix>
    bool bitmat_ops<matrix>::is_lowertriangular(matrix const & A) {
        typename matrix::datatype mask = ~(typename matrix::datatype) 1;
        for(unsigned int k = 0 ; k < matrix::width ; k++, mask<<=1) {
            if (A[k] & mask) return 0;
        }
        return 1;
    }
    template<typename matrix>
    bool bitmat_ops<matrix>::is_uppertriangular(matrix const & A) {
        typename matrix::datatype mask = 1;
        for(unsigned int k = 0 ; k < matrix::width ; k++, mask<<=1) {
            if (A[k]&(mask-1)) return 0;
        }
        return 1;
    }
    template<typename matrix>
    bool bitmat_ops<matrix>::triangular_is_unit(matrix const & A) {
        typename matrix::datatype mask = 1;
        for(unsigned int k = 0 ; k < matrix::width ; k++, mask<<=1) {
            if (!(A[k]&mask)) return 0;
        }
        return 1;
    }
    template<typename matrix>
    void bitmat_ops<matrix>::extract_uppertriangular(matrix & a, matrix const & b) {
        typename matrix::datatype mask = 1;
        for(unsigned int k = 0 ; k < matrix::width ; k++, mask<<=1) {
            a[k] = b[k] & ~(mask-1);
        }
    }
    template<typename matrix>
    void bitmat_ops<matrix>::extract_lowertriangular(matrix & a, matrix const & b) {
        typename matrix::datatype mask = ~(typename matrix::datatype) 1;
        for(unsigned int k = 0 ; k < matrix::width ; k++, mask<<=1) {
            a[k] = b[k] & ~mask;
        }
    }
    template<typename matrix>
    void bitmat_ops<matrix>::extract_LU(matrix & L, matrix & U) {
        typename matrix::datatype mask = 1;
        for(unsigned int k = 0 ; k < matrix::width ; k++, mask<<=1) {
            L[k] = (U[k] & (mask-1));
            U[k] ^= L[k];
            L[k] ^= mask;
        }
    }
    template<typename matrix>
    void bitmat_ops<matrix>::make_uppertriangular(matrix & u) {
        extract_uppertriangular(u, u);
    }
    template<typename matrix>
    void bitmat_ops<matrix>::make_lowertriangular(matrix & u) {
        extract_lowertriangular(u, u);
    }
    template<typename matrix>
    void bitmat_ops<matrix>::make_unit_uppertriangular(matrix & u) {
        make_uppertriangular(u);
        triangular_make_unit(u);
    }
    template<typename matrix>
    void bitmat_ops<matrix>::make_unit_lowertriangular(matrix & u) {
        make_lowertriangular(u);
        triangular_make_unit(u);
    }
    template<typename matrix>
    void bitmat_ops<matrix>::triangular_make_unit(matrix & u) {
        typename matrix::datatype mask = 1;
        for(unsigned int k = 0 ; k < matrix::width ; k++, mask<<=1) {
            u[k] |= mask;
        }
    }
    template<typename matrix>
    void bitmat_ops<matrix>::trsm(matrix const& L, matrix& U)
    {
        trsm(L, U, 0, matrix::width);
    }
    template<typename matrix>
    void bitmat_ops<matrix>::addmul(matrix & C, matrix const & A, matrix const & B)
    {
        addmul(C, A, B, 0, matrix::width, 0, matrix::width);
    }
}


#endif	/* BBLAS_BITMAT_INL_HPP_ */
