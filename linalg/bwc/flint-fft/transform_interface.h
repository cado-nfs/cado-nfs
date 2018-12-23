#ifndef TRANSFORM_INTERFACE_H_
#define TRANSFORM_INTERFACE_H_

#ifdef __cplusplus
extern "C" {
#endif

#define xxxDEBUG_FFT

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
};

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

void fft_get_transform_info(struct fft_transform_info * fti, mp_bitcnt_t bits1, mp_bitcnt_t bits2, unsigned int nacc);
void fft_get_transform_info_mulmod(struct fft_transform_info * fti, mp_bitcnt_t xbits, mp_bitcnt_t ybits, unsigned int nacc, mp_bitcnt_t minwrap);
void fft_transform_info_adjust_depth(struct fft_transform_info * fti, unsigned int adj);
void fft_transform_info_set_first_guess(struct fft_transform_info * fti);
int fft_transform_info_check(const struct fft_transform_info * fti);
void fft_get_transform_allocs(size_t sizes[3], const struct fft_transform_info * fti);
void fft_transform_prepare(void * x, const struct fft_transform_info * fti);
void fft_do_dft(void * y, const mp_limb_t * x, mp_size_t nx, void * temp, const struct fft_transform_info * fti);
void fft_do_ift(mp_limb_t * x, mp_size_t nx, void * y, void * temp, const struct fft_transform_info * fti);
void fft_mul(void * z, const void * y0, const void * y1, void * temp, const struct fft_transform_info * fti);
void fft_addmul(void * z, const void * y0, const void * y1, void * temp, void * qtemp, const struct fft_transform_info * fti);
void fft_add(void * z, void * y0, void * y1, const struct fft_transform_info * fti);
void fft_zero(void * z, const struct fft_transform_info * fti);
void fft_get_transform_info_fppol(struct fft_transform_info * fti, mpz_srcptr p, mp_size_t n1, mp_size_t n2, unsigned int nacc);
void fft_get_transform_info_fppol_mp(struct fft_transform_info * fti, mpz_srcptr p, mp_size_t nmin, mp_size_t nmax, unsigned int nacc);
void fft_do_dft_fppol(void * y, const mp_limb_t * x, mp_size_t nx, void * temp, const struct fft_transform_info * fti, mpz_srcptr p);
void fft_do_ift_fppol(mp_limb_t * x, mp_size_t nx, void * y, void * temp, const struct fft_transform_info * fti, mpz_srcptr p);
void fft_do_ift_fppol_mp(mp_limb_t * x, mp_size_t nx, void * y, void * temp, const struct fft_transform_info * fti, mpz_srcptr p, mp_size_t shift);

/* fft_transform_export modifies the transform area in x and makes it
 * position independent, so that the data may be moved, or transferred to
 * another machine.
 *             /=============================================\
 *             | This (reversibly) invalidates the data in x |
 *             \=============================================/
 * fft_transform_import must be called on x to revert the effect of
 * fft_transform_export (possibly after moving/transferring).
 */
void fft_transform_export(void * x, const struct fft_transform_info * fti);
void fft_transform_import(void * x, const struct fft_transform_info * fti);
/* indicates whether the integer returned is actually reduced modulo some
 * B^n-a, with a=\pm1. Returns n, and sets a. If the result is known to
 * be valid in Z, then n is returned as 0.
 */
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
