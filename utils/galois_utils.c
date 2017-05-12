#include "cado.h"
#include <stdio.h>
#include "portability.h"
#include "utils.h"

struct galois_automorphism_s sigmas[] = {
    /* the homography H={a,b,c,d} encodes (ax+b)/(cx+d) */
    { "autom2.1", 2, 1, 1, {0,1,1,0} },       // 1/x
    { "autom2.2", 2, 1, 1, {-1,0,0,1} },      // -x
    { "autom3.1", 3, 1, 1, {1,-1,1,0} },      // (x-1)/x = 1-1/x
    { "autom3.2", 3, 1, 1, {-1,-1,1,0} },     // (-x-1)/x = -1-1/x
    { "autom4.1", 4, 2, 2, {-1,-1,1,-1} },    // (-x-1)/x = -1-1/x
    { "autom6.1", 6, 3, 2, {-2,-1,1,-1} },    // -(2*x+1)/(x-1)
    // horrible trick for using autom6.1^3 which has order 2
    { "autom6.1Z", 2, 1, 1, {-1,-2,2,1} },    // -(x+2)/(2*x+1)
    { "_1/y", 2, 1, 1, {0,-1,1,0} },          // -1/x
    // some old aliases.
    { "autom2.1g", 2, 1, 1, {0,1,1,0} },       // 1/x
    { "1/y", 2, 1, 1, {0,1,1,0} },             // 1/x
    { "autom2.2g", 2, 1, 1, {-1,0,0,1} },      // -x
    { "_y", 2, 1, 1, {-1,0,0,1} },             // -x
    { "autom3.1g", 3, 1, 1, {1,-1,1,0} },      // (x-1)/x = 1-1/x
    { 0, },
};

galois_automorphism_srcptr galois_automorphism_get(const char * name)
{
    if (name == NULL)
        return NULL;
    for(unsigned int i = 0 ; sigmas[i].name ; i++) {
        if (strcmp(sigmas[i].name, name) == 0) {
            galois_automorphism_srcptr sigma = &(sigmas[i]);
            ASSERT_ALWAYS(sigma->order % sigma->period == 0);
            return sigma;
        }
    }
    return NULL;
}

/* returns how many times sigma->factor must by added to the
 * corresponding relation */
int galois_automorphism_apply_ab(galois_automorphism_srcptr sigma, 
				 int64_t * a, int64_t * b
)
{
    mpz_poly ab;
    mpz_poly_init(ab, 1);
    mpz_poly_setcoeff_int64(ab, 0, *a);
    mpz_poly_setcoeff_int64(ab, 1, *b);
    mpz_poly_homography(ab, ab, sigma->H, 1);
    if (mpz_sgn(ab->coeff[1]) < 0)
        mpz_poly_neg(ab, ab);
    *a = mpz_get_int64(ab->coeff[0]);
    *b = mpz_get_int64(ab->coeff[1]);
    mpz_poly_clear(ab);

    if (sigma->factor <= 1) return 0;

    /* we add [factor] every [period] applications. So the norm gets
     * [factor^degree] every [order] applications, and the degree is equal to
     * the order (we're talking about cyclic galois actions). So each
     * application gains us [order/period] occurences of the prime factor
     * [factor].  We then compensate every time we have a and b both
     * divisible by [factor] (which we don't expect would occur more than
     * once, normally).
     *
     * Note: if the automorphism order is less than the polynomial
     * degree, then it generates a subgroup. Norm contributions should be
     * multiplied by [degree/order], but that is done outside this
     * function, because we're only dealing with the automorphism itself
     * here.
     */
    int e = sigma->order / sigma->period;
    if (*a % sigma->factor == 0 && *b % sigma->factor == 0) {
        e -= sigma->order;
        *a /= sigma->factor;
        *b /= sigma->factor;
    }
    ASSERT_ALWAYS (*a % sigma->factor != 0 || *b % sigma->factor != 0);

    return e;
}

/* Perhaps we want this to deal correctly with projective roots. Do we
 * also care about powers ?
 *
 * Should we do a "fast" version ? Doubt it.
 */
void galois_automorphism_apply_root(galois_automorphism_srcptr sigma, mpz_ptr r1, mpz_srcptr r0, mpz_srcptr p)
{
    // x -> (A*x+B)/(C*x+D)
    mpz_t u,v;
    mpz_init(u);
    mpz_init(v);
    ASSERT_ALWAYS(mpz_cmp(r0, p) <= 0);
    if (mpz_cmp(r0, p) == 0) {
        /* image of infinity is a/c */
        mpz_set_si(v, sigma->H[2]);
        if (mpz_invert(v, v, p)) {
            mpz_mul_si(u, v,  sigma->H[0]);
            mpz_mod(r1, u, p);
        } else {
            mpz_set(r1, p);
        }
    } else {
        mpz_mul_si(v, r0, sigma->H[2]); mpz_add_si(v, v, sigma->H[3]);
        if (mpz_invert(v, v, p)) {
            mpz_mul_si(u, r0, sigma->H[0]); mpz_add_si(u, u, sigma->H[1]);
            mpz_mul(u, u, v);
            mpz_mod(r1, u, p);
        } else {
            mpz_set(r1, p);
        }
    }
    mpz_clear(u);
    mpz_clear(v);
}
