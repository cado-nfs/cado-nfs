#include "cado.h"
#include <algorithm>

#include "arith-cross.hpp"
#include "bblas_level3a.hpp"
#include "bblas_level3a1.hpp"
#include "bblas_level3c.hpp"

#ifdef HAVE_SSE2
#include <x86intrin.h>
#endif

/* code for add_dotprod {{{ */
/* Parameter L may be fixed at the template level, and used as a way to
 * special-case the code. The default version does not do that, though.
 */
template<unsigned int, unsigned int>
struct add_dotprod {
    void operator()(
            uint64_t * w,           // 64L at a time
            const uint64_t * u,     // 64K at a time
            const uint64_t * v,     // 64L at a time
            unsigned int n,
            unsigned int K,
            unsigned int L) const
    {
        /* u has n rows of 64K bits
         * v has n rows of 64L bits.
         * Compute (u|v) == tr(u)*v into the area pointed to by v: 64K
         * rows of 64L bits.
         */
        /* This is meant for L > 2 -- whether K == 1 or not does not
         * matter much, as usual. We're not using sse-2 here because of
         * lazyness -- and little expected returns.
         */
        /* Apparently we don't have this code in bblas at this point.
         */
        for(unsigned int i = 0 ; i < n ; i++) {
            uint64_t * w0 = (uint64_t *) w;
            for(unsigned int l = 0 ; l < L ; l++) {
                uint64_t b = *v++;
                uint64_t *sw = w0;
                const uint64_t * u0 = u;
                for(unsigned int k = 0 ; k < K ; k++) {
                    uint64_t a = *u0++;
                    for (unsigned int j = 0; j < 64; j ++) {
                        *sw ^= b & -(a & 1);
                        a >>= 1;
                        sw += L;
                    }
                }
                w0++;
            }
            u += K;
        }
    }
};

/* u has n rows of 64K bits
 * v has n rows of 64 bits.
 * Compute (u|v) == tr(u)*v into the area pointed to by v: 64K rows of 64 bits.
 */
template<>
struct add_dotprod<0,1> {
    inline void operator()(
            uint64_t * w,
            const uint64_t * u,
            const uint64_t * v, unsigned int n,
            unsigned int K,
            unsigned int /* L=1 is hard-coded */) const
    {
        addmul_TN64K_N64(w, u, v, n, K);
    }
};

#if    defined(HAVE_SSE2) && GMP_LIMB_BITS == 64
/* u has n rows of 64K bits
 * v has n rows of 128 bits.
 * Compute (u|v) == tr(u)*v into the area pointed to by v: 64K rows of 128 bits.
 */
template<>
struct add_dotprod<0,2> {
    inline void operator()(
            uint64_t * w,           // 128 at a time
            const uint64_t * u,     // 64K at a time
            const uint64_t * v,     // 128 at a time
            unsigned int n,
            unsigned int K,
            unsigned int /* L=2 is hard-coded */) const
    {
        /* okay, we converted the m128i* to u64* for interchange, and now
         * we're casting it back...  */
        /* Apparently we don't have this code in bblas at this point.
         */
        __m128i const * vv = reinterpret_cast<const __m128i*>(v);
        for(unsigned int i = 0 ; i < n ; i++) {
            __m128i * w0 = (__m128i*) w;
            __m128i mb[4][2] = {
                {_mm_setzero_si128(), _mm_setzero_si128()},
                {*vv, _mm_setzero_si128()},
                {_mm_setzero_si128(), *vv},
                {*vv, *vv},
            };
            vv++;
            __m128i *sw = w0;
            for(unsigned int k = 0 ; k < K ; k++) {
                uint64_t a = *u++;
                for (unsigned int j = 0; j < 64; j += 2) {
                    _mm_storeu_si128(sw, _mm_xor_si128(_mm_loadu_si128(sw), mb[a & 3][0]));
                    sw ++;
                    _mm_storeu_si128(sw, _mm_xor_si128(_mm_loadu_si128(sw), mb[a & 3][1]));
                    sw ++;
                    a >>= 2;
                }
            }
        }
    }
};
#endif

/* Here is code for another specialization. (for <K,2>)
 *
static inline void add_dotprod_64K_128(uint64_t * b, const uint64_t * A, const uint64_t * x, unsigned int ncol, unsigned int K)
{
    uint64_t idx, i, rA;
    uint64_t rx[2];

    abort();    // untested.

    for(idx = 0; idx < ncol; idx++) {
        rx[0] = *x++;
        rx[1] = *x++;
        uint64_t* pb = b;
        for(unsigned int j = 0 ; j < K ; j++) {
            rA = *A++;
            for(i = 0; i < 64; i++) {
                *pb++ ^= rx[0] & -(rA & 1);
                *pb++ ^= rx[1] & -(rA & 1);
                rA >>= 1;
            }
        }
    }
}
*/

/* }}} */

/* {{{ code for addmul_tiny */
template<unsigned int, unsigned int>
struct addmul_tiny {
    /* multiply the n times 64K-bit vector u by the 64K by
     * 64L matrix v -- n must be even. Result is put in w.
     */
    void operator()(
            uint64_t * w,           // 64L at a time
            const uint64_t * u,     // 64K at a time
            const uint64_t * v,     // 64L at a time
            unsigned int n,
            unsigned int K,
            unsigned int L) const
    {
#if     defined(HAVE_SSE2) && GMP_LIMB_BITS == 64
        ASSERT_ALWAYS((n & 1) == 0);
        /* We're doing addmul, so don't clear the output. */
        /* We'll read and write rows two at a time. So we maintain two
         * pointers for each, interlaced. */
        uint64_t * u0 = (uint64_t *) u;
        uint64_t * u1 = (uint64_t *) (u + K);
        uint64_t * w0 = (uint64_t *) w;
        uint64_t * w1 = (uint64_t *) (w + L);
        for (unsigned int j = 0; j < n; j += 2 ) {
            const uint64_t * v0 = v;
            for(unsigned int l = 0 ; l < L ; l++) {
                __m128i r = _mm_setzero_si128();
                const uint64_t * vv = v0;
                for(unsigned int k = 0 ; k < K ; k++) {
                    __m128i a = _mm_setr_epi64(_mm_cvtsi64_m64(u0[k]), _mm_cvtsi64_m64(u1[k]));
                    __m128i one = _mm_set1_epi64(_mm_cvtsi64_m64(1));
                    for (unsigned int i = 0; i < 64; i++) {
                        __m128i zw = _mm_set1_epi64(_mm_cvtsi64_m64(*vv));
                        r = _mm_xor_si128(r, _mm_and_si128(zw, _mm_sub_epi64(_mm_setzero_si128(), _mm_and_si128(a, one))));
                        a = _mm_srli_epi64(a, 1);
                        vv += L;
                    }
                }
                /* Because L > 1, the destination words are not contiguous */
                w0[l] ^= _mm_cvtsi128_si64(r);
                w1[l] ^= _mm_cvtsi128_si64(_mm_srli_si128(r, 8));
                v0++; /* next column in v */
            }
            u0 += 2 * K; u1 += 2 * K;
            w0 += 2 * L; w1 += 2 * L;
        }
#else   /* HAVE_SSE2 */
        /* This is just a direct translation of the sse version. Has been
         * tested once. */
        uint64_t * u0 = (uint64_t *) u;
        uint64_t * w0 = (uint64_t *) w;
        for (unsigned int j = 0; j < n; j ++ ) {
            const uint64_t * v0 = v;
            for(unsigned int l = 0 ; l < L ; l++) {
                uint64_t rx = 0;
                const uint64_t * vv = v0;
                for(unsigned int k = 0 ; k < K ; k++) {
                    uint64_t a = u0[k];
                    for (unsigned int i = 0; i < 64; i++) {
                        rx ^= (*vv & -(a & (uint64_t) 1));
                        a >>= 1;
                        vv += L;
                    }
                }
                /* Because L > 1, the destination words are not contiguous */
                w0[l] ^= rx;
                v0++; /* next column in v */
            }
            u0 += K;
            w0 += L;
        }
#endif
    }
};
/* }}} */

/* Transpose a matrix of 64K rows of 64L bits */
struct transpose {
    inline void operator()(uint64_t * w, const uint64_t * u, unsigned int K, unsigned int L) const
    {
#if 1
        /* We have 1*(64*K)*L*1 64-bit words */
        /* We want 1*K*L*64*1 64-bit words */
        uint64_t* temp = new uint64_t[K*L*64];
        if (u != w)
            generic_transpose_words<uint64_t>(w, u, 1, 64, K*L, 1);
        else
            generic_transpose_words_inplace<uint64_t>(w, 1, 64, K*L, 1, temp);
        for (unsigned int k = 0; k < K * L; k++) {
            mat64& d(reinterpret_cast<mat64&>(w[k]));
            mat64_transpose(d, d);
        }
        /* We have 1*K*L*64*1 64-bit words */
        /* We want 1*L*(64*K)*1 64-bit words */
        generic_transpose_words_inplace<uint64_t>(w, 1, K, L*64, 1, temp);
        delete[] temp;
#else
        uint64_t md = 1UL;
        unsigned int od = 0;
        uint64_t ms = 1UL;
        unsigned int os = 0;
        uint64_t * dst = (uint64_t *) w;
        const uint64_t * src = (const uint64_t *) u;
        for(unsigned int l = 0 ; l < 64*L ; l++) {
            for(unsigned int k = 0 ; k < 64*K ; k++) {
                dst[od] |= md & -((src[os+k*L] & ms) != 0);
                md <<= 1;
                od += md == 0;
                md += md == 0;
            }
            ms <<= 1;
            os += ms == 0;
            ms += ms == 0;
        }
#endif
    }
};


template<unsigned int K, unsigned int L> struct arith_cross_gf2;

template<unsigned int L>
struct arith_cross_gf2<0, L>
: public arith_cross_generic
{
    unsigned int g0;
    unsigned int g1;
    arith_cross_gf2(arith_generic const * A0, arith_generic const * A1)
        : g0(A0->simd_groupsize())
        , g1(A1->simd_groupsize())
    {
        ASSERT_ALWAYS(L == 0 || g1 == L * 64);
    }
    virtual void add_dotprod(arith_generic::elt * w, arith_generic::elt const * u, arith_generic::elt const * v, unsigned int n) const override
    {
        auto xw = reinterpret_cast<uint64_t *>(w);
        auto xv = reinterpret_cast<const uint64_t *>(v);
        auto xu = reinterpret_cast<const uint64_t *>(u);
        ::add_dotprod<0,L>()(xw,xu,xv,n,g0/64,g1/64);
    }

    virtual void addmul_tiny(arith_generic::elt * w, arith_generic::elt const * u, arith_generic::elt const * v, unsigned int n) const override {
        auto xw = reinterpret_cast<uint64_t *>(w);
        auto xv = reinterpret_cast<const uint64_t *>(v);
        auto xu = reinterpret_cast<const uint64_t *>(u);
        ::addmul_tiny<0,L>()(xw,xu,xv,n,g0/64,g1/64);
    }

    virtual void transpose(arith_generic::elt * w, arith_generic::elt const * u) const override {
        auto xw = reinterpret_cast<uint64_t *>(w);
        auto xu = reinterpret_cast<const uint64_t *>(u);
        ::transpose()(xw,xu,g0/64,g1/64);
    }

    virtual ~arith_cross_gf2() override = default;
};

struct arith_cross_gfp : public arith_cross_generic
{
    arith_generic * A;
    arith_cross_gfp(arith_generic * A0, arith_generic const * A1)
        : A(A0)
    {
        ASSERT_ALWAYS(A1->simd_groupsize() == A->simd_groupsize());
        ASSERT_ALWAYS(A->simd_groupsize() == 1);
        ASSERT_ALWAYS(mpz_cmp(A->characteristic(), A1->characteristic()) == 0);
        ASSERT_ALWAYS(A->impl_name() == A1->impl_name());
    }
    virtual void add_dotprod(arith_generic::elt * w, arith_generic::elt const * u, arith_generic::elt const * v, unsigned int n) const override
    {
        A->vec_add_dotprod(*w, u, v, n);
    }

    virtual void addmul_tiny(arith_generic::elt * w, arith_generic::elt const * u, arith_generic::elt const * v, unsigned int n) const override
    {
        A->vec_addmul_and_reduce(w, u, *v, n);
    }

    virtual void transpose(arith_generic::elt * w, arith_generic::elt const * u) const override
    {
        A->set(*w, *u);
    }

    virtual ~arith_cross_gfp() override = default;
};

arith_cross_generic * arith_cross_generic::instance(arith_generic * A0, arith_generic * A1)
{
    if (A0->is_characteristic_two() && A1->is_characteristic_two()) {
        if (A1->simd_groupsize() == 64)
            return new arith_cross_gf2<0,1>(A0, A1);
        else if (A1->simd_groupsize() == 128)
            return new arith_cross_gf2<0,2>(A0, A1);
        else
            return new arith_cross_gf2<0,0>(A0, A1);
    } else if (mpz_cmp(A0->characteristic(), A1->characteristic()) == 0) {
        if (A0->simd_groupsize() != 1) return NULL;
        if (A1->simd_groupsize() != 1) return NULL;
        return new arith_cross_gfp(A0,A1);
    } else {
        return NULL;
    }
}
