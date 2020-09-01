#ifndef BBLAS_BITMAT_INL_HPP_
#define BBLAS_BITMAT_INL_HPP_

// IWYU pragma: private, include "bblas.hpp"
// IWYU pragma: friend ".*/bblas.*"

#include "bblas.hpp"
#include "bblas_bitmat.hpp"
#include "gmp_aux.h"

namespace bblas_bitmat_details {
    template<typename T>
    void bitmat_ops<T>::fill_random(bitmat<T> & w, gmp_randstate_t rstate)
    {
        memfill_random(w.data(), sizeof(bitmat<T>), rstate);
    }

    template<typename T>
    bool bitmat_ops<T>::is_lowertriangular(bitmat<T> const & A) {
        T mask = ~(T) 1;
        for(unsigned int k = 0 ; k < bitmat<T>::width ; k++, mask<<=1) {
            if (A[k] & mask) return 0;
        }
        return 1;
    }
    template<typename T>
    bool bitmat_ops<T>::is_uppertriangular(bitmat<T> const & A) {
        T mask = 1;
        for(unsigned int k = 0 ; k < bitmat<T>::width ; k++, mask<<=1) {
            if (A[k]&(mask-1)) return 0;
        }
        return 1;
    }
    template<typename T>
    bool bitmat_ops<T>::triangular_is_unit(bitmat<T> const & A) {
        T mask = 1;
        for(unsigned int k = 0 ; k < bitmat<T>::width ; k++, mask<<=1) {
            if (!(A[k]&mask)) return 0;
        }
        return 1;
    }
    template<typename T>
    void bitmat_ops<T>::extract_uppertriangular(bitmat<T> & a, bitmat<T> const & b) {
        T mask = 1;
        for(unsigned int k = 0 ; k < bitmat<T>::width ; k++, mask<<=1) {
            a[k] = b[k] & ~(mask-1);
        }
    }
    template<typename T>
    void bitmat_ops<T>::extract_lowertriangular(bitmat<T> & a, bitmat<T> const & b) {
        T mask = ~(T) 1;
        for(unsigned int k = 0 ; k < bitmat<T>::width ; k++, mask<<=1) {
            a[k] = b[k] & ~mask;
        }
    }
    template<typename T>
    void bitmat_ops<T>::extract_LU(bitmat<T> & L, bitmat<T> & U) {
        T mask = 1;
        for(unsigned int k = 0 ; k < bitmat<T>::width ; k++, mask<<=1) {
            L[k] = (U[k] & (mask-1));
            U[k] ^= L[k];
            L[k] ^= mask;
        }
    }
    template<typename T>
    void bitmat_ops<T>::make_uppertriangular(bitmat<T> & u) {
        extract_uppertriangular(u, u);
    }
    template<typename T>
    void bitmat_ops<T>::make_lowertriangular(bitmat<T> & u) {
        extract_lowertriangular(u, u);
    }
    template<typename T>
    void bitmat_ops<T>::make_unit_uppertriangular(bitmat<T> & u) {
        make_uppertriangular(u);
        triangular_make_unit(u);
    }
    template<typename T>
    void bitmat_ops<T>::make_unit_lowertriangular(bitmat<T> & u) {
        make_lowertriangular(u);
        triangular_make_unit(u);
    }
    template<typename T>
    void bitmat_ops<T>::triangular_make_unit(bitmat<T> & u) {
        T mask = 1;
        for(unsigned int k = 0 ; k < bitmat<T>::width ; k++, mask<<=1) {
            u[k] |= mask;
        }
    }
    template<typename T>
    void bitmat_ops<T>::trsm(bitmat<T> const& L, bitmat<T>& U)
    {
        trsm(L, U, 0, bitmat<T>::width);
    }
    template<typename T>
    void bitmat_ops<T>::addmul(bitmat<T> & C, bitmat<T> const & A, bitmat<T> const & B)
    {
        addmul(C, A, B, 0, bitmat<T>::width, 0, bitmat<T>::width);
    }
}


#endif	/* BBLAS_BITMAT_INL_HPP_ */
