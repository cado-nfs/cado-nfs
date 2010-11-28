#ifndef ABASE_BINARY_DOTPROD_BACKENDS_H_
#define ABASE_BINARY_DOTPROD_BACKENDS_H_

/* TODO: there's some duplicated code between here and matops.c */

#ifdef  need_dotprod_64K_64

#ifdef  HAVE_SSE2
#include <emmintrin.h>
/* u has n rows of 64K bits
 * v has n rows of 64 bits.
 * Compute (u|v) == tr(u)*v into the area pointed to by v: 64K rows of 64 bits.
 */
static inline void dotprod_64K_64(
        uint64_t * w,           // 64 at a time
        const uint64_t * u,     // 64K at a time
        const uint64_t * v,     // 64 at a time
        unsigned int n,
        unsigned int K)
{
    memset(w, 0, 64 * K * sizeof(uint64_t));
    for(unsigned int i = 0 ; i < n ; i++) {
        __v2di * w0 = (__v2di*) w;
        // TODO: It's possible to expand more, and use a __v2di
        // mb[4][2], or even [4]. This wouldn't change the code much
        // (see the u128 version), and is likely to speed things up a
        // wee bit maybe.
        __v2di mb[4] = {
            (__v2di) {0, 0},
            (__v2di) {*v, 0},
            (__v2di) {0, *v},
            (__v2di) {*v, *v},
        };
        v++;
        __v2di *sw = w0;
        for(unsigned int k = 0 ; k < K ; k++) {
            uint64_t a = *u++;
            for (unsigned int j = 0; j < 64; j += 2) {
                *sw ^= mb[a & 3];
                a >>= 2;
                sw ++;
            }
        }
    }
}
#else
static inline void dotprod_64K_64(uint64_t * b, const uint64_t * A, const uint64_t * x, unsigned int ncol, unsigned int K)
{
    uint64_t idx, i, rA;
    uint64_t rx;

    /* This has been tested once, seems to work ok */

    memset(b, 0, 64 * K * sizeof(uint64_t));
    for(idx = 0; idx < ncol; idx++) {
        rx = x[idx];
        uint64_t* pb = b;
        for(unsigned int j = 0 ; j < K ; j++) {
            rA = *A++;
            for(i = 0; i < 64; i++) {
                *pb++ ^= rx & -(rA & 1);
                rA >>= 1;
            }
        }
    }
}
#endif
#endif

#ifdef  need_dotprod_64K_128
#ifdef HAVE_SSE2
/* u has n rows of 64K bits
 * v has n rows of 128 bits.
 * Compute (u|v) == tr(u)*v into the area pointed to by v: 64K rows of 128 bits.
 */
#include <emmintrin.h>
static inline void dotprod_64K_128(
        uint64_t * w,           // 128 at a time
        const uint64_t * u,     // 128 at a time
        const uint64_t * v,     // 128 at a time
        unsigned int n,
        unsigned int K)
{
    memset(w, 0, 64 * K * sizeof(__v2di));
    /* okay, we've casted the v2di* to u64* for interchange, and now
     * we're casting it back...  */
    __v2di * vv = (__v2di*) v;
    for(unsigned int i = 0 ; i < n ; i++) {
        __v2di * w0 = (__v2di*) w;
#define SSE_0   ((__v2di) {0,0} )
        __v2di mb[4][2] = {
            {SSE_0, SSE_0},
            {*vv, SSE_0},
            {SSE_0, *vv},
            {*vv, *vv},
        };
#undef SSE_0
        vv++;
        __v2di *sw = w0;
        for(unsigned int k = 0 ; k < K ; k++) {
            uint64_t a = *u++;
            for (unsigned int j = 0; j < 64; j += 2) {
                *sw ^= mb[a & 3][0];
                sw ++;
                *sw ^= mb[a & 3][1];
                sw ++;
                a >>= 2;
            }
        }
    }
}
#else
static inline void dotprod_64K_128(uint64_t * b, const uint64_t * A, const uint64_t * x, unsigned int ncol, unsigned int K)
{
    uint64_t idx, i, rA;
    uint64_t rx[2];

    abort();    // untested.

    memset(b, 0, 64 * K * sizeof(uint64_t));
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
#endif
#endif

#ifdef  need_dotprod_64K_64L
/* u has n rows of 64K bits
 * v has n rows of 64L bits.
 * Compute (u|v) == tr(u)*v into the area pointed to by v: 64K rows of 64L bits.
 */
/* This is meant for L > 2 -- whether K == 1 or not does not matter much,
 * as usual. We're not using sse-2 here because of lazyness -- and little
 * expected returns.
 */
static void dotprod_64K_64L(
        uint64_t * w,           // 64L at a time
        const uint64_t * u,     // 64K at a time
        const uint64_t * v,     // 64L at a time
        unsigned int n,
        unsigned int K,
        unsigned int L)
{
    memset(w, 0, 64 * K * L * sizeof(uint64_t));
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
#endif

#ifdef need_vaddmul_tiny_64K_64L
/* multiply the n times 64K-bit vector u by the 64K by
 * 64L matrix v -- n must be even. Result is put in w.
 */
#ifdef HAVE_SSE2
#include <emmintrin.h>
static inline void vaddmul_tiny_64K_64L(
            uint64_t * w,
            const uint64_t * u,
            const uint64_t * v,
            unsigned int n,
            unsigned int K,
            unsigned int L)
{
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
            union { __v2di s; uint64_t x[2]; } r;
            r.s = (__v2di) {0,0};
            const uint64_t * vv = v0;
            for(unsigned int k = 0 ; k < K ; k++) {
                __v2di a = (__v2di) { u0[k], u1[k] };
                __v2di one = (__v2di) { 1, 1, };
                for (unsigned int i = 0; i < 64; i++) {
                    __v2di zw = { *vv, *vv, };
                    r.s ^= (zw & -(a & one));
                    a = _mm_srli_epi64(a, 1);
                    vv += L;
                }
            }
            /* Because L > 1, the destination words are not contiguous */
            w0[l] ^= r.x[0];
            w1[l] ^= r.x[1];
            v0++; /* next column in v */
        }
        u0 += 2 * K; u1 += 2 * K;
        w0 += 2 * L; w1 += 2 * L;
    }
}
#else   /* HAVE_SSE2 */
static inline void vaddmul_tiny_64K_64L(
            uint64_t * w,
            const uint64_t * u,
            const uint64_t * v,
            unsigned int n,
            unsigned int K,
            unsigned int L)
{
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
}
#endif
#endif

#ifdef need_vtranspose_64K_64L
/* Transpose a matrix of 64K rows of 64L bits */
static inline void vtranspose_64K_64L(
            uint64_t * w,
            const uint64_t * u,
            unsigned int K, unsigned int L)
{
    memset(w, 0, 64 * K * L * sizeof(uint64_t));
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
}
#endif

#endif	/* ABASE_BINARY_DOTPROD_BACKENDS_H_ */
