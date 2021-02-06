#ifndef TRANSFORM_INTERFACE_H_
#define TRANSFORM_INTERFACE_H_

#include <string.h>

#include <gmp.h>
#include "flint.h"
#ifdef __cplusplus
#include <array>
#include <string>
#endif

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

/* Initialize a transform type that is suitable to hold the sum of (nacc)
 * products of integers of size (bits1) and (bits2)
 */
void fft_transform_info_init(struct fft_transform_info * fti, mp_bitcnt_t bits1, mp_bitcnt_t bits2, unsigned int nacc);

/* Initialize a transform type that is suitable to hold the sum of (nacc)
 * products of integers of size (bits1) and (bits2), modulo some 2^n-1
 * with n>=minwrap. (the actual value of n may be freely chosen by the
 * implementation).
 */
void fft_transform_info_init_mulmod(struct fft_transform_info * fti, mp_bitcnt_t xbits, mp_bitcnt_t ybits, unsigned int nacc, mp_bitcnt_t minwrap);

/* This can be used to provide an external means to tune the parameter
 * choices of the implementation. It's called "adjust depth", but what
 * this actually means is up to the choice of the implementation. The
 * value adj==UINT_MAX means that no adjustment is to be done.
 */
void fft_transform_info_adjust_depth(struct fft_transform_info * fti, unsigned int adj);

/* This can be used to somehow revert the effect of adjust_depth, and
 * return to a default choice as chosen by one of the info_init
 * functions. Note that it is mandatory that fti has already been
 * initialized by one of the info_init functions !
 */
void fft_transform_info_set_first_guess(struct fft_transform_info * fti);

/* Perform sanity checks on the fti structure */
int fft_transform_info_check(const struct fft_transform_info * fti);

/* This returns the needed temporary storage for the different steps of a
 * fft multiplication process whose size corresponds
 * to the info given in *fti.
 * [0]: space to be allocated (size_t) for each transform.
 * [1]: temp space to be passed alongside with each transform
 *      (needs be allocated only once). This same amount is also needed
 *      when calling fft_addcompose
 * [2]: temp space to be passed alongside with each fft_compose or
 *      fft_addcompose convolution (needs be allocated only once).
 *
 * spaces returned in [1] and [2] are independent and the same area may
 * be used for both, but of course the caller must then ensure that they
 * are not used concurrently.
 *
 * Note that the size returned in [0] for each transform entails slight
 * overallocation due to pointer swap tricks here and there.
 */
void fft_transform_info_get_alloc_sizes(const struct fft_transform_info * fti, size_t sizes[3]);

/* prepares the allocated memory area pointed by x to hold a transform.
 * This might be a noop.
 */
void fft_prepare(const struct fft_transform_info * fti, void * x);

/* Computes in y the DFT of the integer {x, nx}. temp must point to a
 * memory are whose size is sizes[1], as returned by
 * fft_transform_info_get_alloc_sizes
 */
void fft_dft(const struct fft_transform_info * fti, void * y, const mp_limb_t * x, mp_size_t nx, void * temp);

/* Computes in {x, nx} the IFT of the transform y. temp must point to a
 * memory whose size is sizes[1], as returned by
 * fft_transform_info_get_alloc_sizes. Any data fast limb nx is
 * truncated.
 */
void fft_ift(const struct fft_transform_info * fti, mp_limb_t * x, mp_size_t nx, void * y, void * temp);


/* Computes in z the convolution product of the two transforms pointed to
 * by y and y1.  temp must point to a memory area whose size is sizes[2],
 * as returned by fft_transform_info_get_alloc_sizes
 */
void fft_compose(const struct fft_transform_info * fti, void * z, const void * y0, const void * y1, void * temp);

/* Adds to z the convolution product of the two transforms pointed to
 * by y and y1.  temp must point to a memory area whose size is sizes[2],
 * as returned by fft_transform_info_get_alloc_sizes qtemp must point to
 * a memory area whose size is sizes[1], as returned by
 * fft_transform_info_get_alloc_sizes
 */
void fft_addcompose(const struct fft_transform_info * fti, void * z, const void * y0, const void * y1, void * temp, void * qtemp);

/* z = y0 + y1 */
void fft_add(const struct fft_transform_info * fti, void * z, const void * y0, const void * y1);

/* zeroes out z. The are pointed to by z may be any sort of arbitrary
 * garbage.
 */
void fft_zero(const struct fft_transform_info * fti, void * z);

/* Fill z with random data from the random state (note that since random
 * state is not MT-safe, it might be tricky to do this random filling in
 * parallel).
 */
void fft_fill_random(const struct fft_transform_info * fti, void * z, gmp_randstate_t rstate);

/* Provide the transform info necessary for accumulating nacc products of
 * polynomials having respectively n1 and n2 coefficients, over GF(p),
 * using Kronecker substitution.
 * (we hare talking *length* n1 and n2, hence degrees n1-1 and n2-1).
 *
 * Up to nacc accumulated products should be supported by the
 * returned transform type.
 */
void fft_transform_info_init_fppol(struct fft_transform_info * fti, mpz_srcptr p, mp_size_t n1, mp_size_t n2, unsigned int nacc);

/* middle product of two polynomials of length nmin and nmax, with nmin
 * <= nmax.
 * Up to nacc accumulated products should be supported by the
 * returned transform type.
 */
void fft_transform_info_init_fppol_mp(struct fft_transform_info * fti, mpz_srcptr p, mp_size_t nmin, mp_size_t nmax, unsigned int nacc);

/* check that the area pointed to by x is a valid transform. if (diag) is
 * true, print diagnostics on stderr.
 */
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

/* In the mulmod case, return the integer so that the product is computed
 * mod 2^n-1 (see fft_transform_info_init_mulmod)
 */
static inline mp_bitcnt_t fft_get_mulmod(const struct fft_transform_info * fti, int * a);

/* In the mulmod case, return the minimum number of limbs that are
 * required to store the computed integer.  (see
 * fft_transform_info_init_mulmod and fft_get_mulmod)
 */
static inline mp_size_t fft_get_mulmod_output_minlimbs(const struct fft_transform_info * fti);

/* Returns a malloc()ed string (or NULL) providing explanation of what
 * this transform type is doing.
 */
char * fft_transform_info_explain(const struct fft_transform_info * fti);

#ifdef __cplusplus
}
#endif

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
    inline std::array<size_t, 3> get_alloc_sizes() const {
        std::array<size_t, 3> sizes;
        fft_transform_info_get_alloc_sizes(this, &sizes[0]);
        return sizes;
    }
    inline size_t size0_bytes() const { return get_alloc_sizes()[0]; }
    inline size_t size1_bytes() const { return get_alloc_sizes()[1]; }
    inline size_t size2_bytes() const { return get_alloc_sizes()[2]; }
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
    std::string explain() const { char * x = fft_transform_info_explain(this); std::string s = x; free(x); return s; }
#endif
};

#ifdef __cplusplus
extern "C" {
#endif

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
