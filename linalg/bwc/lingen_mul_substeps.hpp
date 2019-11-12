#ifndef LINGEN_MUL_SUBSTEPS_HPP_
#define LINGEN_MUL_SUBSTEPS_HPP_

/* Use matpoly as imported from lingen_matpoly_ft.hpp */
#include "lingen_matpoly_ft.hpp"

/* middle product and multiplication are really the same thing, so better
 * avoid code duplication */

template<typename fft_type>
struct op_mul {
    typedef fft_type FFT;
    fft_type fti;
    size_t csize = 0;
    static constexpr const char * name = "MUL";
    op_mul() = default;
#ifdef SELECT_MPFQ_LAYER_u64k1
    op_mul(size_t asize, size_t bsize, unsigned int adj MAYBE_UNUSED = UINT_MAX) : fti(fft_type::mul_info(asize, bsize))
    {
        csize = asize + bsize; csize -= (csize > 0);
        if (adj != UINT_MAX)
            fti.adjust(GF2X_FFT_ADJUST_DEPTH, adj);
    }
    template<typename T>
    op_mul(T const & a, T const & b, unsigned int adj = UINT_MAX) : op_mul(a.get_size(), b.get_size(), adj)
    {}
#else
    op_mul(mpz_srcptr p, size_t asize, size_t bsize, unsigned int nacc, unsigned int adj = UINT_MAX) : fti(fft_type::polynomial_mul_info(p, asize, bsize, nacc))
    {
        csize = asize + bsize; csize -= (csize > 0);
        if (adj != UINT_MAX) fti.adjust_depth(adj);
    }
    template<typename T>
    op_mul(T const & a, T const & b, unsigned int adj = UINT_MAX) : op_mul(abfield_characteristic_srcptr(a.ab), a.get_size(), b.get_size(), a.n, adj)
    {}
#endif
};

template<typename fft_type>
struct op_mp {
    typedef fft_type FFT;
    fft_type fti;
    size_t csize = 0;
    static constexpr const char * name = "MP";
    unsigned int shift = 0;
    op_mp() = default;
#ifdef SELECT_MPFQ_LAYER_u64k1
    op_mp(size_t asize, size_t bsize, unsigned int MAYBE_UNUSED adj = UINT_MAX) : fti(fft_type::mp_info(asize, bsize))
    {
        csize = MAX(asize, bsize) - MIN(asize, bsize) + 1;
        shift = MIN(asize, bsize) - 1;
        if (adj != UINT_MAX)
            fti.adjust(GF2X_FFT_ADJUST_DEPTH, adj);
    }
    template<typename T>
    op_mp(T const & a, T const & b, unsigned int adj = UINT_MAX) : op_mp(a.get_size(), b.get_size(), adj)
    {}
#else
    op_mp(mpz_srcptr p, size_t asize, size_t bsize, unsigned int nacc, unsigned int adj = UINT_MAX) : fti(fft_type::polynomial_mp_info(p, asize, bsize, nacc))
    {
        csize = MAX(asize, bsize) - MIN(asize, bsize) + 1;
        shift = MIN(asize, bsize) - 1;
        if (adj != UINT_MAX) fti.adjust_depth(adj);
    }
    template<typename T>
    op_mp(T const & a, T const & b, unsigned int adj = UINT_MAX)
        : op_mp(abfield_characteristic_srcptr(a.ab), a.get_size(), b.get_size(), a.n, adj)
    {}
#endif
};

/* specializations for the no-fft case */
template<>
struct op_mul<void> {
    typedef void FFT;
    size_t csize = 0;
    static constexpr const char * name = "MUL";
    op_mul() = default;
#ifdef SELECT_MPFQ_LAYER_u64k1
    op_mul(size_t asize, size_t bsize, unsigned int adj MAYBE_UNUSED = UINT_MAX)
    {
        csize = asize + bsize; csize -= (csize > 0);
    }
    template<typename T>
    op_mul(T const & a, T const & b, unsigned int adj = UINT_MAX) : op_mul(a.get_size(), b.get_size(), adj)
    {}
#else
    mpz_srcptr p;
    op_mul(mpz_srcptr p, size_t asize, size_t bsize, unsigned int nacc MAYBE_UNUSED, unsigned int adj MAYBE_UNUSED = UINT_MAX) : p(p)
    {
        csize = asize + bsize; csize -= (csize > 0);
    }
    template<typename T>
    op_mul(T const & a, T const & b, unsigned int adj = UINT_MAX) : op_mul(abfield_characteristic_srcptr(a.ab), a.get_size(), b.get_size(), a.n, adj)
    {}
#endif
    static void addcompose(matpoly::view_t t, matpoly::const_view_t t0, matpoly::const_view_t t1) {
        matpoly::addmul(t, t0, t1);
    }
};

template<>
struct op_mp<void> {
    typedef void FFT;
    size_t csize = 0;
    static constexpr const char * name = "MP";
    unsigned int shift = 0;
    op_mp() = default;
#ifdef SELECT_MPFQ_LAYER_u64k1
    op_mp(size_t asize, size_t bsize, unsigned int MAYBE_UNUSED adj = UINT_MAX)
    {
        csize = MAX(asize, bsize) - MIN(asize, bsize) + 1;
        shift = MIN(asize, bsize) - 1;
    }
    template<typename T>
    op_mp(T const & a, T const & b, unsigned int adj = UINT_MAX) : op_mp(a.get_size(), b.get_size(), adj)
    {}
#else
    mpz_srcptr p;
    op_mp(mpz_srcptr p, size_t asize, size_t bsize, unsigned int nacc MAYBE_UNUSED, unsigned int adj MAYBE_UNUSED = UINT_MAX) : p(p)
    {
        csize = MAX(asize, bsize) - MIN(asize, bsize) + 1;
        shift = MIN(asize, bsize) - 1;
    }
    template<typename T>
    op_mp(T const & a, T const & b, unsigned int adj = UINT_MAX)
        : op_mp(abfield_characteristic_srcptr(a.ab), a.get_size(), b.get_size(), a.n, adj)
    {}
#endif
    static void addcompose(matpoly::view_t t, matpoly::const_view_t t0, matpoly::const_view_t t1) {
        matpoly::addmp(t, t0, t1);
    }
};


#endif	/* LINGEN_MUL_SUBSTEPS_HPP_ */
