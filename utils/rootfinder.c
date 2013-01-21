#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <gmp.h>

#include "macros.h"
#include "rootfinder.h"
#include "plain_poly.h"
#include "modul_poly.h"
#include "poly.h"

/* put in r[0], ..., r[n-1] the roots of f (of degree d) modulo p,
   and the return value n is the number of roots (without multiplicities) */
int
poly_roots_ulong (unsigned long *r, mpz_t *f, int d, unsigned long p)
{
    int n;

    if (plain_poly_fits (d, p)) {
        plain_poly_coeff_t * rr;
        int i;

        if (r == NULL)
          return plain_poly_roots(NULL, f, d, p);

        rr = (plain_poly_coeff_t *) malloc(d * sizeof(plain_poly_coeff_t));
        n = plain_poly_roots(rr, f, d, p);
        for(i = 0 ; i < n ; i++) {
            r[i] = rr[i] < 0 ? (rr[i]+(plain_poly_coeff_t) p) : rr[i];
        }
        free(rr);
    } else {
        residueul_t * rr;
        modulusul_t pp;
        modul_initmod_ul(pp, p);
        int i;

        if (r == NULL)
            return modul_poly_roots(NULL, f, d, pp);

        rr = (residueul_t *) malloc(d * sizeof(residueul_t));
        for(i = 0 ; i < d ; i++) {
            modul_init_noset0(rr[i], pp);
        }
        n = modul_poly_roots(rr, f, d, pp);
        for(int i = 0 ; i < n ; i++) {
            /* The assumption is that p fits within an unsigned long
             * anyway. So the roots do as well.
             */
            r[i] = modul_get_ul(rr[i], pp);
        }
        for(i = 0 ; i < d ; i++) {
            modul_clear(rr[i], pp);
        }
        free(rr);
        modul_clearmod(pp);
    }
    return n;
}

int poly_roots_long(long * r, mpz_t * f, int d, unsigned long p)
{
    FATAL_ERROR_CHECK(p >= LONG_MAX, "Cannot use poly_roots_long here");
    return poly_roots_ulong((unsigned long *) r, f, d, p);
}

int
poly_roots_uint64 (uint64_t * r, mpz_t * f, int d, uint64_t p)
{
    FATAL_ERROR_CHECK(p >= ULONG_MAX, "poly_roots_uint64 is a stub");
    /* This is glue around poly_roots_ulong, nothing more. When uint64
     * means asking a lot more than ulong, we miss code. */
    unsigned long *rr;
    int i, n;

    if (r == NULL)
      return poly_roots_ulong(NULL, f, d, p);

    if (sizeof (unsigned long) != sizeof (uint64_t))
      rr = (unsigned long *) malloc(d * sizeof(unsigned long));
    else
      rr = (unsigned long *) r;
    n = poly_roots_ulong (rr, f, d, p);
    if (sizeof (unsigned long) != sizeof (uint64_t))
      {
        for(i = 0 ; i < n ; i++)
          r[i] = rr[i];
        free (rr);
      }
    return n;
}

typedef int (*sortfunc_t) (const void *, const void *);

static int compare_mpz_ptrs(
        const mpz_t * a,
        const mpz_t * b)
{
    if (mpz_cmp(*a, *b) < 0) return -1;
    if (mpz_cmp(*b, *a) < 0) return 1;
    return 0;
}

int poly_roots(mpz_t * r, mpz_t * f, int d, mpz_t p)
{
    if (mpz_cmp_ui(p, ULONG_MAX) <= 0) {
        /* There's a chance of using one of our layers. */
        unsigned long pp = mpz_get_ui(p);
        unsigned long * rr;
        int i;
        int n;

        if (r == NULL)
            return poly_roots_ulong(NULL, f, d, pp);

        rr = (unsigned long *) malloc(d * sizeof(unsigned long));
        n = poly_roots_ulong(rr,  f, d, pp);

        for(i = 0 ; i < n ; i++) {
            /* The assumption is that p fits within an unsigned long
             * anyway. So the roots do as well.
             */
            mpz_set_ui(r[i], rr[i]);
        }
        free(rr);
        return n;
    } else {
      int n;
      n = poly_roots_mpz (r, f, d, p);
      if (r && n)
          qsort(r, n, sizeof(mpz_t), (sortfunc_t) &compare_mpz_ptrs);
      return n;
    }
}
