#ifndef TRANSFORM_INTERFACE_H_
#define TRANSFORM_INTERFACE_H_

#include <string.h>

#include "gmp.h"
#include "flint.h"

#define xxxDEBUG_FFT

struct fft_transform_info;

/* The transform data is provided as an opaque void* pointer ; in reality
 * its structure may be written as the following pseudo-C code.
 *  struct transform_data {
 *      union {
 *          mp_limb_t * x[4 << depth];
 *          mp_limb_t * xy[2 << depth2][1 << depth1];
 *          // 1 + depth = depth1 + depth2
 *      } p;
 *      mp_limb_t * extra[2];
 *      mp_limb_t data[0];
 *  };
 * with pointers in x[], or xy[][], or extra[], all pointing to areas
 * within the zone beginning at data[]. Those areas are coefficients in
 * R=Z/(2^(nw)+1), each occupying fti_rsize0(fti)+1 limbs. All these
 * coefficient areas are disjoint. xy[] and x[] are of course to ways of
 * accessing the same set of pointers.

 * A transform data object represents a sequence of 4n coefficients in R
 * (pointed to by x[]). Before a DFT, this may be intepreted as a
 * polynomial P modulo x^(4n)-1. After, it's a sequence of pointwise
 * evaluations of P.  The layout is to be interpreted as follows.

 * with fti->alg == 0:
 *     x[i] == P( sqrt(2)^(bitrev(i, depth+2)) )
 * with fti->alg == 1 (matrix algorithm), The meaning of xy[] differs in
 * the two halves of the array.  We use the notation n1 == 1<<depth1,
 * n2==1<<depth2. The contents of xy[][] are given as follows.  Let
 * i<n2. 
 
 * xy[bitrev(i,depth2)][j]             == P(sqrt(2)^(    2*(j + i<<depth1))
 *                                        P(             2^(j + i<<depth1))
 * xy[bitrev(i,depth2) + 1<<depth2][j] == P(sqrt(2)^(1 + 2*(j + i<<depth1)))
 *                                        P(   sqrt(2) * 2^(j + i<<depth1))

 * In case of truncation, only some of the entries are considered. For
 * fti->alg==0, these are the entries with i < trunc. For fti->alg==1, the
 * rule is instead (i+b<<depth2) < trunc, with b being the most significant
 * bit of the row index.
 */

#ifdef __cplusplus
extern "C" {
#endif

void fft_transform_info_init(struct fft_transform_info * fti, mp_bitcnt_t bits1, mp_bitcnt_t bits2, unsigned int nacc);
void fft_transform_info_init_mulmod(struct fft_transform_info * fti, mp_bitcnt_t xbits, mp_bitcnt_t ybits, unsigned int nacc, mp_bitcnt_t minwrap);
void fft_transform_info_adjust_depth(struct fft_transform_info * fti, unsigned int adj);
void fft_transform_info_set_first_guess(struct fft_transform_info * fti);
int fft_transform_info_check(const struct fft_transform_info * fti);
void fft_transform_info_get_alloc_sizes(const struct fft_transform_info * fti, size_t sizes[3]);

void fft_prepare(const struct fft_transform_info * fti, void * x);
void fft_dft(const struct fft_transform_info * fti, void * y, const mp_limb_t * x, mp_size_t nx, void * temp);
void fft_ift(const struct fft_transform_info * fti, mp_limb_t * x, mp_size_t nx, void * y, void * temp);
void fft_compose(const struct fft_transform_info * fti, void * z, const void * y0, const void * y1, void * temp);
void fft_addcompose(const struct fft_transform_info * fti, void * z, const void * y0, const void * y1, void * temp, void * qtemp);
void fft_add(const struct fft_transform_info * fti, void * z, const void * y0, const void * y1);
void fft_zero(const struct fft_transform_info * fti, void * z);
void fft_fill_random(const struct fft_transform_info * fti, void * z, gmp_randstate_t rstate);


void fft_transform_info_init_fppol(struct fft_transform_info * fti, mpz_srcptr p, mp_size_t n1, mp_size_t n2, unsigned int nacc);
void fft_transform_info_init_fppol_mp(struct fft_transform_info * fti, mpz_srcptr p, mp_size_t nmin, mp_size_t nmax, unsigned int nacc);
int fft_check(const struct fft_transform_info * fti, const void * x, int);

/* fft_transform_export modifies the transform area in x and makes it
 * position independent, so that the data may be moved, or transferred to
 * another machine.
 *             /=============================================\
 *             | This (reversibly) invalidates the data in x |
 *             \=============================================/
 * fft_transform_import must be called on x to revert the effect of
 * fft_transform_export (possibly after moving/transferring).
 */
void fft_export(const struct fft_transform_info * fti, void * x);
void fft_import(const struct fft_transform_info * fti, void * x);
/* indicates whether the integer returned is actually reduced modulo some
 * B^n-a, with a=\pm1. Returns n, and sets a. If the result is known to
 * be valid in Z, then n is returned as 0.
 */


struct fft_transform_info {
    mp_bitcnt_t bits1;
    mp_bitcnt_t bits2;
    unsigned int nacc;
    mp_size_t w;        /* Use sqrt(2)^w as a root of unity */
    mp_size_t depth;    /* Let n=2^depth. We work modulo 2^(wn)+1. Do a
                           transform length 4n. */
    mp_bitcnt_t bits;   /* Chunk sizes in bits */
    mp_size_t trunc0;   /* Number of coeffs of the transform computed.
                           This is not exactly the fourier transform
                           truncation point, because the truncation point
                           is also subject to a few extra criteria. */
    mp_bitcnt_t ks_coeff_bits;  /* This is used only for kronecker substitution */
    mp_bitcnt_t minwrap;        /* zero when no wraparound wanted */
    int alg;            /* alg==1: use matrix fourier algorithm */
    mpz_srcptr p;       /* non-NULL if we're doing Kronecker substitution
                           to multiply polynomials over GF(p). This is
                           the (pointer to) prime that got passed to the
                           ctor, therefore it is essential that this
                           prime remains alive throughout the lifetime of
                           the fft_transform_info object !!! */
    unsigned int mp_shift;      /* non-zero if we're doing a middle
                           product of polynomials over GF(p).  This is
                           the number of GF(p) coefficients that are
                           taken out from the result, and for which any
                           wraparound garbage is not important to us. */
#ifdef DEBUG_FFT
    char tmpdir[FILENAME_MAX];
#endif
#ifdef __cplusplus
    typedef void * ptr;
    typedef const void * srcptr;
    /* only for uniformity with callers that deal with other interfaces -- it's not even clear we'll need that. */
    fft_transform_info() { memset(this, 0, sizeof(*this)); }
    inline fft_transform_info(mp_bitcnt_t bits1, mp_bitcnt_t bits2, unsigned int nacc) {
        fft_transform_info_init(this, bits1, bits2, nacc);
    }
    static inline fft_transform_info mul_info(mp_bitcnt_t bits1, mp_bitcnt_t bits2, unsigned int nacc) {
        return fft_transform_info(bits1, bits2, nacc);
    }
    static inline fft_transform_info mulmod_info(mp_bitcnt_t xbits, mp_bitcnt_t ybits, unsigned int nacc, mp_bitcnt_t minwrap) {
        fft_transform_info fti;
        fft_transform_info_init_mulmod(&fti, xbits, ybits, nacc, minwrap);
        return fti;
    }
    static inline fft_transform_info polynomial_mul_info(mpz_srcptr p, mp_size_t n1, mp_size_t n2, unsigned int nacc)
    {
        fft_transform_info fti;
        fft_transform_info_init_fppol(&fti, p, n1, n2, nacc);
        return fti;
    }
    static inline fft_transform_info polynomial_mp_info(mpz_srcptr p, mp_size_t nmin, mp_size_t nmax, unsigned int nacc)
    {
        fft_transform_info fti;
        fft_transform_info_init_fppol_mp(&fti, p, nmin, nmax, nacc);
        return fti;
    }
    inline void adjust_depth(unsigned int adj) { fft_transform_info_adjust_depth(this, adj); }
    inline void get_alloc_sizes(size_t sizes[3]) const {
        fft_transform_info_get_alloc_sizes(this, sizes);
    }
    inline void prepare(ptr x) const { fft_prepare(this, x); }
    inline void dft(ptr y, const mp_limb_t * x, mp_size_t nx, ptr temp) const {
        fft_dft(this, y, x, nx, temp);
    }
    inline void ift(mp_limb_t * x, mp_size_t nx, ptr y, ptr temp) const {
        fft_ift(this, x, nx, y, temp);
    }
    inline void add(ptr z, srcptr y0, srcptr y1) const {
        return fft_add(this, z, y0, y1);
    }
    inline void compose(ptr z, srcptr y0, srcptr y1, ptr temp) const {
        fft_compose(this, z, y0, y1, temp);
    }
    inline void addcompose(ptr z, srcptr y0, srcptr y1, ptr temp, ptr qtemp) const {
        fft_addcompose(this, z, y0, y1, temp, qtemp);
    }
    inline void zero(ptr x) const { fft_zero(this, x); }
    inline void fill_random(ptr x, gmp_randstate_t rstate) const { fft_fill_random(this, x, rstate); }
    inline int check(srcptr x, int c) const { return fft_check(this, x, c); }
    inline void to_export(ptr x) const { fft_export(this, x); }
    inline void to_import(ptr x) const { fft_import(this, x); }
#endif
};

static inline mp_bitcnt_t fft_get_mulmod(const struct fft_transform_info * fti, int * a)
{
    *a=1;
    return fti->minwrap ? (4<<fti->depth)*fti->bits : 0;
}

static inline mp_size_t fft_get_mulmod_output_minlimbs(const struct fft_transform_info * fti)
{
    if (!fti->minwrap) return 0;
    mp_size_t w = fti->w;
    mp_size_t n = 1 << fti->depth;
    mp_bitcnt_t need = (4*n-1)*fti->bits+n*w;
    return (need + FLINT_BITS - 1) / FLINT_BITS;
}

#ifdef __cplusplus
}
#endif

#endif	/* TRANSFORM_INTERFACE_H_ */
