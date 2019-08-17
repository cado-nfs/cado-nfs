#ifndef LINGEN_MUL_SUBSTEPS_HPP_
#define LINGEN_MUL_SUBSTEPS_HPP_

#include "lingen_matpoly.hpp"
#include "lingen_matpoly_ft.hpp"
#include "flint-fft/transform_interface.h"

/* middle product and multiplication are really the same thing, so better
 * avoid code duplication */

struct op_mul {
    fft_transform_info fti[1];
    size_t csize = 0;
    static constexpr const char * name = "MUL";
    op_mul() = default;
    op_mul(mpz_srcptr p, size_t asize, size_t bsize, unsigned int nacc, unsigned int adj = UINT_MAX)
    {
        csize = asize + bsize; csize -= (csize > 0);
        fft_get_transform_info_fppol(fti, p, asize, bsize, nacc);
        if (adj != UINT_MAX)
            fft_transform_info_adjust_depth(fti, adj);
    }
    template<typename T>
    op_mul(T const & a, T const & b, unsigned int adj = UINT_MAX) : op_mul(a.ab->p, a.size, b.size, a.n, adj)
    {}
    // op_mul(bigmatpoly const & a, bigmatpoly const & b, fft_transform_info * fti, unsigned int adj = UINT_MAX) : op_mul(a.ab->p, a.size, b.size, a.n, fti, adj)
    // {}

    inline void ift(matpoly::view_t a, matpoly_ft::view_t t) const
    {
        ::ift(a, t);
    }

    size_t get_transform_ram() const {
        size_t fft_alloc_sizes[3];
        fft_get_transform_allocs(fft_alloc_sizes, fti);
        return fft_alloc_sizes[0];
    }
};

struct op_mp {
    fft_transform_info fti[1];
    size_t csize = 0;
    static constexpr const char * name = "MP";
    unsigned int shift = 0;
    op_mp() = default;
    op_mp(mpz_srcptr p, size_t asize, size_t bsize, unsigned int nacc, unsigned int adj = UINT_MAX)
    {
        csize = MAX(asize, bsize) - MIN(asize, bsize) + 1;
        shift = MIN(asize, bsize) - 1;
        fft_get_transform_info_fppol_mp(fti, p, MIN(asize, bsize), MAX(asize, bsize), nacc);
        if (adj != UINT_MAX)
            fft_transform_info_adjust_depth(fti, adj);
    }
    template<typename T>
    op_mp(T const & a, T const & b, unsigned int adj = UINT_MAX)
        : op_mp(a.ab->p, a.size, b.size, a.n, adj)
    {}
    // op_mp(bigmatpoly const & a, bigmatpoly const & b, fft_transform_info * fti, unsigned int adj = UINT_MAX) : op_mp(a.ab->p, a.size, b.size, a.n, fti, adj)
    // {}

    inline void ift(matpoly::view_t a, matpoly_ft::view_t t) const
    {
        ::ift_mp(a, t, shift);
    }

    size_t get_transform_ram() const {
        size_t fft_alloc_sizes[3];
        fft_get_transform_allocs(fft_alloc_sizes, fti);
        return fft_alloc_sizes[0];
    }
};


#endif	/* LINGEN_MUL_SUBSTEPS_HPP_ */
