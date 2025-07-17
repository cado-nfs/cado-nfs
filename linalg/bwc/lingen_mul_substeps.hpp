#ifndef CADO_LINGEN_MUL_SUBSTEPS_HPP
#define CADO_LINGEN_MUL_SUBSTEPS_HPP

#include <cstddef>
#include <climits>

#include <string>
#include <algorithm>
#include <array>

#include <gmp.h>

#include "lingen_matpoly_select.hpp"
#include "lingen_fft_select.hpp"
#include "lingen_mul_substeps_base.hpp"
#include "macros.h"

template<bool is_binary, typename fft_type>
struct op_mul;

template<typename fft_type>
struct op_mul_inner : public op_mul_or_mp_base {
    typedef fft_type FFT;
    fft_type fti;
    size_t csize = 0;
    std::array<size_t, 3> get_alloc_sizes() const override { return fti.get_alloc_sizes(); }
    std::string explain() const override { return fti.explain(); }
    op_mul_inner(size_t asize, size_t bsize, fft_type fti)
        : op_mul_or_mp_base(OP_MUL)
        , fti(std::move(fti))
        , csize(asize + bsize - (asize || bsize))
    {}
};

template<>
struct op_mul_inner<void> : public op_mul_or_mp_base {
    typedef void FFT;
    size_t csize = 0;
    op_mul_inner(size_t asize, size_t bsize)
        : op_mul_or_mp_base { OP_MUL }
        , csize(asize + bsize - (asize || bsize))
    {}
    std::array<size_t, 3> get_alloc_sizes() const override { return {{ 0, 0, 0 }}; }
    std::string explain() const override { return "plain product"; }
};

#ifdef LINGEN_BINARY
template<typename fft_type>
struct op_mul<true, fft_type> : public op_mul_inner<fft_type> {
    static constexpr bool is_binary = true;
    // static constexpr const char * name = "MUL";
    op_mul() = default;
    op_mul(size_t asize, size_t bsize, unsigned int adj MAYBE_UNUSED = UINT_MAX)
        : op_mul_inner<fft_type>(asize, bsize, fft_type::mul_info(asize, bsize))
    {
        if (adj != UINT_MAX)
            op_mul_inner<fft_type>::fti.adjust(GF2X_FFT_ADJUST_DEPTH, adj);
    }
    op_mul(mpz_srcptr, size_t asize, size_t bsize, unsigned int, unsigned int adj = UINT_MAX) : op_mul(asize, bsize, adj) {}
    template<typename T>
    op_mul(T const & a, T const & b, unsigned int adj = UINT_MAX) : op_mul(a.get_size(), b.get_size(), adj)
    {}
};
#endif

#ifndef LINGEN_BINARY
template<typename fft_type>
struct op_mul<false, fft_type> : public op_mul_inner<fft_type> {
    static constexpr bool is_binary = false;
    // static constexpr const char * name = "MUL";
    op_mul() = default;
    op_mul(mpz_srcptr p, size_t asize, size_t bsize, unsigned int nacc, unsigned int adj = UINT_MAX)
        : op_mul_inner<fft_type>(asize, bsize,
                fft_type::polynomial_mul_info(p, asize, bsize, nacc))
    {
        if (adj != UINT_MAX)
            op_mul_inner<fft_type>::fti.adjust_depth(adj);
    }
    template<typename T>
    op_mul(T const & a, T const & b, unsigned int adj = UINT_MAX)
        : op_mul(a.ab->characteristic(), a.get_size(), b.get_size(), a.n, adj)
    {}
};
#endif

template<bool is_binary, typename fft_type>
struct op_mp;

template<typename fft_type>
struct op_mp_inner : public op_mul_or_mp_base {
    typedef fft_type FFT;
    fft_type fti;
    size_t csize = 0;
    unsigned int shift = 0;
    op_mp_inner() = default;
    op_mp_inner(size_t asize, size_t bsize, fft_type fti)
        : op_mul_or_mp_base { OP_MP }
        , fti(std::move(fti))
        , csize(std::max(asize, bsize) - std::min(asize, bsize) + 1)
        , shift(std::min(asize, bsize) - 1)
    {}
    std::array<size_t, 3> get_alloc_sizes() const override { return fti.get_alloc_sizes(); }
    std::string explain() const override { return fti.explain(); }
};

template<>
struct op_mp_inner<void> : public op_mul_or_mp_base {
    typedef void FFT;
    size_t csize = 0;
    // static constexpr const char * name = "MP";
    unsigned int shift = 0;
    op_mp_inner(size_t asize, size_t bsize)
        : op_mul_or_mp_base { OP_MP }
        , csize(std::max(asize, bsize) - std::min(asize, bsize) + 1)
        , shift(std::min(asize, bsize) - 1)
    {
    }
    std::array<size_t, 3> get_alloc_sizes() const override { return {{ 0, 0, 0 }}; }
    std::string explain() const override { return "plain product"; }
};


#ifdef LINGEN_BINARY
template<typename fft_type>
struct op_mp<true, fft_type> : public op_mp_inner<fft_type>
{
    static_assert(is_binary_fft<fft_type>::value);
    static constexpr bool is_binary = true;
    op_mp(size_t asize, size_t bsize, unsigned int MAYBE_UNUSED adj = UINT_MAX)
        : op_mp_inner<fft_type>(asize, bsize, fft_type::mp_info(asize, bsize))
    {
        if (adj != UINT_MAX)
            op_mp_inner<fft_type>::fti.adjust(GF2X_FFT_ADJUST_DEPTH, adj);
    }
    op_mp(mpz_srcptr,
            size_t asize, size_t bsize,
            unsigned int, unsigned int adj = UINT_MAX)
        : op_mp(asize, bsize, adj)
    {}
    template<typename T>
    op_mp(T const & a, T const & b, unsigned int adj = UINT_MAX)
        : op_mp(a.get_size(), b.get_size(), adj)
    {}
};
#endif

#ifndef LINGEN_BINARY
template<typename fft_type>
struct op_mp<false, fft_type> : public op_mp_inner<fft_type>
{
    static_assert(!is_binary_fft<fft_type>::value);
    static constexpr bool is_binary = false;
    op_mp(mpz_srcptr p, size_t asize, size_t bsize,
            unsigned int nacc, unsigned int adj = UINT_MAX)
        : op_mp_inner<fft_type>(asize, bsize,
                fft_type::polynomial_mp_info(p, asize, bsize, nacc))
    {
        if (adj != UINT_MAX)
            op_mp_inner<fft_type>::fti.adjust_depth(adj);
    }
    template<typename T>
    op_mp(T const & a, T const & b, unsigned int adj = UINT_MAX)
        : op_mp(a.ab->characteristic(), a.get_size(), b.get_size(), a.n, adj)
    {}
};
#endif

#ifdef LINGEN_BINARY
template<>
struct op_mul<true, void> : public op_mul_inner<void> {
    static constexpr bool is_binary = true;
    op_mul(size_t asize, size_t bsize, unsigned int adj MAYBE_UNUSED = UINT_MAX)
        : op_mul_inner(asize, bsize)
    {
    }
    op_mul(mpz_srcptr, size_t asize, size_t bsize,
            unsigned int, unsigned int adj = UINT_MAX)
        : op_mul(asize, bsize, adj)
    {}
    template<typename T>
    op_mul(T const & a, T const & b, unsigned int adj = UINT_MAX)
        : op_mul(a.get_size(), b.get_size(), adj)
    {}
    static void addcompose(matpoly<true>::view_t t, matpoly<true>::const_view_t t0, matpoly<true>::const_view_t t1) {
        matpoly<true>::addmul(t, t0, t1);
    }
};
#endif

#ifndef LINGEN_BINARY
/* specializations for the no-fft case */
template<>
struct op_mul<false, void> : public op_mul_inner<void> {
    static constexpr bool is_binary = false;
    mpz_srcptr p;
    op_mul(mpz_srcptr p, size_t asize, size_t bsize,
            unsigned int nacc MAYBE_UNUSED,
            unsigned int adj MAYBE_UNUSED = UINT_MAX)
        : op_mul_inner(asize, bsize)
        , p(p)
    { }
    template<typename T>
    op_mul(T const & a, T const & b, unsigned int adj = UINT_MAX)
        : op_mul(a.ab->characteristic(), a.get_size(), b.get_size(), a.n, adj)
    {}
    static void addcompose(matpoly<false>::view_t t, matpoly<false>::const_view_t t0, matpoly<false>::const_view_t t1) {
        matpoly<false>::addmul(t, t0, t1);
    }
};
#endif

#ifdef LINGEN_BINARY
template<>
struct op_mp<true, void> : public op_mp_inner<void> {
    static constexpr bool is_binary = true;
    op_mp(size_t asize, size_t bsize, unsigned int MAYBE_UNUSED adj = UINT_MAX)
        : op_mp_inner<void>(asize, bsize)
    {}
    op_mp(mpz_srcptr, size_t asize, size_t bsize,
            unsigned int, unsigned int adj = UINT_MAX)
        : op_mp(asize, bsize, adj)
    {}
    template<typename T>
    op_mp(T const & a, T const & b, unsigned int adj = UINT_MAX)
        : op_mp(a.get_size(), b.get_size(), adj)
    {}
    static void addcompose(matpoly<true>::view_t t, matpoly<true>::const_view_t t0, matpoly<true>::const_view_t t1) {
        matpoly<true>::addmp(t, t0, t1);
    }
};
#endif

#ifndef LINGEN_BINARY
template<>
struct op_mp<false, void> : public op_mp_inner<void> {
    static constexpr bool is_binary = false;
    mpz_srcptr p;
    op_mp(mpz_srcptr p,
            size_t asize, size_t bsize,
            unsigned int nacc MAYBE_UNUSED,
            unsigned int adj MAYBE_UNUSED = UINT_MAX)
        : op_mp_inner<void>(asize, bsize)
        , p(p)
    {
    }
    template<typename T>
    op_mp(T const & a, T const & b, unsigned int adj = UINT_MAX)
        : op_mp(a.ab->characteristic(), a.get_size(), b.get_size(), a.n, adj)
    {}
    static void addcompose(matpoly<false>::view_t t, matpoly<false>::const_view_t t0, matpoly<false>::const_view_t t1) {
        matpoly<false>::addmp(t, t0, t1);
    }
};
#endif

#endif	/* LINGEN_MUL_SUBSTEPS_HPP_ */
