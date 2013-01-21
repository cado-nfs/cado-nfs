#include <math.h>       /* rint etc */
#include "las-types.h"
#include "las-qlattice.h"
#include "las-arith.h"

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
 * is a0^2+sigma*b0^2 */
int generic_skew_gauss(mpz_t a[2], mpz_t b[2], mpz_t sigma)
{
    mpz_t N[2], S, q, tmp;
    mpz_init(N[0]);
    mpz_init(N[1]);
    mpz_init(S);
    mpz_init(q);
    mpz_init(tmp);

    /* Compute the two norms, and the dot products */
    mpz_mul(tmp, a[0], a[0]); mpz_set(N[0], tmp);
    mpz_mul(tmp, b[0], b[0]); mpz_addmul(N[0], tmp, sigma);

    mpz_mul(tmp, a[1], a[1]); mpz_set(N[1], tmp);
    mpz_mul(tmp, b[1], b[1]); mpz_addmul(N[1], tmp, sigma);

    mpz_mul(tmp, a[0], a[1]); mpz_set(S, tmp);
    mpz_mul(tmp, b[0], b[1]); mpz_addmul(S, tmp, sigma);

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

    /* We don't care about the sign of b. Down the road, there's an
     * IJToAB function which guarantees positive b */

    /* However we do care about vector 0 being the ``smallest'' in some
     * sense. The trick is that the comparison criterion used previously
     * by the code was the skewed L^\infty norm max(a,b*skewness). We
     * want the reduction step here to be oblivious to this L^\infty /
     * L^2 mess, so we provide something which is L^2 minimal. Wrapper
     * code around this function is still guaranteeing L^\infinity
     * minimality for compatibility, though (see qlattice_safetyguard)
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

void mpz_set_square_of_d(mpz_t sigma, double skewness)
{
    /* Truncate the square, not the square root... */
    mpq_t qsigma;
    mpq_init(qsigma);
    mpq_set_d(qsigma, skewness);
    mpq_mul(qsigma, qsigma, qsigma);
    mpz_rdiv_q(sigma, mpq_numref(qsigma), mpq_denref(qsigma));
    mpq_clear(qsigma);
}

/* FIXME: I have doubts that this code really belongs here */
int qlattice_safetyguard(sieve_info_ptr si, double skewness)
{
    /* compare skewed max-norms */
    double maxab1 = MAX(fabs(si->a1), fabs(si->b1) * skewness);
    double maxab0 = MAX(fabs(si->a0), fabs(si->b0) * skewness);
    if (maxab0 > maxab1) {
        int64_t oa[2] = { si->a0, si->a1 };
        int64_t ob[2] = { si->b0, si->b1 };
        si->a0 = oa[1]; si->a1 = oa[0];
        si->b0 = ob[1]; si->b1 = ob[0];
        maxab1 = maxab0;
    }
    /* make sure J does not exceed I/2 */
    if (maxab1 >= si->B)
        si->J = (uint32_t) (si->B * skewness / maxab1 * (double) (si->I >> 1));
    else
        si->J = si->I >> 1;

    /* Make sure the bucket region size divides the sieve region size, 
       partly covered bucket regions may lead to problems when 
       reconstructing p from half-empty buckets. */
    /* Compute number of i-lines per bucket region, must be integer */
    ASSERT_ALWAYS(LOG_BUCKET_REGION >= si->logI);
    uint32_t i = 1U << (LOG_BUCKET_REGION - si->logI);
    i *= si->nb_threads;  /* ensures nb of bucket regions divisible by nb_threads */
    si->J = ((si->J - 1U) / i + 1U) * i; /* Round up to multiple of i */

    return 0;
}

int SkewGauss(sieve_info_ptr si, double skewness)
{
    mpz_t sigma;
    mpz_init(sigma);
    mpz_set_square_of_d(sigma, skewness);

    mpz_t a[2], b[2];
    mpz_init_set(a[0], si->q);
    mpz_init_set(a[1], si->rho);
    mpz_init_set_ui(b[0], 0);
    mpz_init_set_ui(b[1], 1);
    generic_skew_gauss(a, b, sigma);
    int fits = 1;
    fits = fits && mpz_cmp_si(a[0], INT64_MIN) >= 0;
    fits = fits && mpz_cmp_si(b[0], INT64_MIN) >= 0;
    fits = fits && mpz_cmp_si(a[1], INT64_MIN) >= 0;
    fits = fits && mpz_cmp_si(b[1], INT64_MIN) >= 0;
    fits = fits && mpz_cmp_si(a[0], INT64_MAX) <= 0;
    fits = fits && mpz_cmp_si(b[0], INT64_MAX) <= 0;
    fits = fits && mpz_cmp_si(a[1], INT64_MAX) <= 0;
    fits = fits && mpz_cmp_si(b[1], INT64_MAX) <= 0;
    if (fits) {
        si->a0 = mpz_get_int64(a[0]);
        si->a1 = mpz_get_int64(a[1]);
        si->b0 = mpz_get_int64(b[0]);
        si->b1 = mpz_get_int64(b[1]);
        qlattice_safetyguard(si, skewness);
    }
    mpz_clear(a[0]);
    mpz_clear(a[1]);
    mpz_clear(b[0]);
    mpz_clear(b[1]);
    mpz_clear(sigma);
    /* FIXME: That error convention looks odd */
    return fits ? 0 : 1;
}
