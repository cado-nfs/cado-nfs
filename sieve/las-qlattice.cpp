#include "cado.h"
#include <math.h>       /* rint etc */
#include "las-types.hpp"
#include "las-qlattice.hpp"
#include "las-arith.hpp"

/* check that the double x fits into an int32_t */
#define fits_int32_t(x) \
  ((double) INT32_MIN <= (x)) && ((x) <= (double) INT32_MAX)

int mpz_rdiv_q(mpz_ptr q, mpz_t a, mpz_t b)
{
    /* Return the relative integer q which is closest to a/b. We
     * guarantee -1/2<=a/b-q<1/2.
     * It's a pity that this is not supported directly by gmp... */
    mpz_t r;
    mpz_init(r);

    mpz_fdiv_qr(q,r,a,b);
    /* b>0: We want -b/2 <= a-bq < b/2 */
    /* b<0: We want  b/2 < a-bq <= -b/2 */
    mpz_mul_2exp(r,r,1);
    if (mpz_cmp(r,b) * mpz_sgn(b) >= 0) {
        mpz_add_ui(q, q, 1);
    }
    mpz_clear(r);

    return 1;
}

/* We work with two vectors v0=(a0,b0) and v1=(a1,b1). The quadratic form
 * is proportional to a0^2+skewness^2*b0^2 */
static int 
generic_skew_gauss(mpz_t a[2], mpz_t b[2], double skewness)
{
    // skewness^2 can be written as a one-word (actually, 53 bits)
    // integer times a power of two, which is presumably larger than
    // 2^-53 in the most common case, since we expect skewness > 1.
    double mantissa;
    int64_t mantissa_z;
    int exponent;
    /* These are provided by c99 */
    mantissa = frexp(skewness * skewness, &exponent);
    mantissa_z = ldexp(mantissa, 53);
    exponent -= 53;

    mpz_t N[2], S, q, tmp;
    mpz_init(N[0]);
    mpz_init(N[1]);
    mpz_init(S);
    mpz_init(q);
    mpz_init(tmp);

    /* Compute the two norms, and the dot products */
#define QUADFORM(dst, i0, i1) do {					\
    mpz_mul(dst, a[i0], a[i1]);						\
    mpz_mul(tmp, b[i0], b[i1]);						\
    if (exponent < 0) mpz_mul_2exp(dst, dst, -exponent);		\
    if (exponent > 0) mpz_mul_2exp(tmp, tmp, exponent);			\
    mpz_addmul_int64(dst, tmp, mantissa_z);				\
} while (0)
    
    QUADFORM(N[0], 0, 0);
    QUADFORM(N[1], 1, 1);
    QUADFORM(S,    0, 1);

    /* After a reduction step (e.g. v0-=q*v1), N[0], N[1], and S are
     * updated using the following algorithm.
     * new_N0 = old_N0 + q^2old_N1 - 2q old_S
     * new_S = old_S - q old_N1
     *
     * i.e.
     *
     * new_N0 = old_N0 - q old_S - q*(old_S-q old_N1)
     * new_N0 = old_N0 - q * (old_S + new_S)
     *
     * cost: two products only
     */
    for(;;) {
        /* reduce v0 with respect to v1 */
        mpz_rdiv_q(q, S, N[1]);
        if (mpz_cmp_ui(q, 0) == 0) break;
        mpz_submul(a[0], q, a[1]);
        mpz_submul(b[0], q, b[1]);
        /* update */
        mpz_set(tmp, S);
        mpz_submul(S, q, N[1]);
        mpz_add(tmp, tmp, S);
        mpz_submul(N[0], q, tmp);

        /* reduce v1 with respect to v0 */
        mpz_rdiv_q(q, S, N[0]);
        if (mpz_cmp_ui(q, 0) == 0) break;
        mpz_submul(a[1], q, a[0]);
        mpz_submul(b[1], q, b[0]);
        /* update */
        mpz_set(tmp, S);
        mpz_submul(S, q, N[0]);
        mpz_add(tmp, tmp, S);
        mpz_submul(N[1], q, tmp);
    }

    /* We don't care about the sign of b. Normalization should be done
     * down the road (e.g. with IJToAB) */

    /* However we do care about vector 0 being the ``smallest'' in some
     * sense. The trick is that the comparison criterion used previously
     * by the code was the skewed L^\infty norm max(a,b*skewness). We
     * want the reduction step here to be oblivious to this L^\infty /
     * L^2 mess, so we provide something which is L^2 minimal. Wrapper
     * code around this function is still guaranteeing L^\infinity
     * minimality for compatibility, though (see sieve_info_adjust_IJ)
     */
    if (mpz_cmp(N[0], N[1]) > 0) {
        mpz_swap(a[0], a[1]);
        mpz_swap(b[0], b[1]);
    }
    mpz_clear(N[0]);
    mpz_clear(N[1]);
    mpz_clear(S);
    mpz_clear(q);
    mpz_clear(tmp);
    return 1;
}

int
SkewGauss (qlattice_basis &basis, const mpz_t p, const mpz_t r,
           const double skewness)
{
    mpz_t a[2], b[2];
    int fits;

    mpz_init_set (a[0], p);
    mpz_init_set (a[1], r);
    mpz_neg (a[1], a[1]);
    mpz_init_set_ui (b[0], 0);
    mpz_init_set_ui (b[1], 1);
    generic_skew_gauss (a, b, skewness);
    fits = mpz_fits_int64_p (a[0]);
    fits = fits && mpz_fits_int64_p (b[0]);
    fits = fits && mpz_fits_int64_p (a[1]);
    fits = fits && mpz_fits_int64_p (b[1]);
    if (fits)
      {
        basis.a0 = mpz_get_int64 (a[0]);
        basis.a1 = mpz_get_int64 (a[1]);
        basis.b0 = mpz_get_int64 (b[0]);
        basis.b1 = mpz_get_int64 (b[1]);
      }
    mpz_clear (a[0]);
    mpz_clear (a[1]);
    mpz_clear (b[0]);
    mpz_clear (b[1]);
    return fits ? 1 : 0;
}
