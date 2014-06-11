#include "cado.h"
#include <stdio.h>
#include "portability.h"
#include "utils.h"

/* Init polynomial rel and set it to a - b*x */
void
mpz_poly_init_set_ab (mpz_poly_ptr rel, int64_t a, uint64_t b)
{
  mpz_poly_init(rel, 1);
  mpz_poly_setcoeff_int64(rel, 0, a);
  mpz_poly_setcoeff_int64(rel, 1, -b);
}

/* compute the SM */
void
compute_sm (mpz_poly_t SM, mpz_poly_t num, const mpz_poly_t F, const mpz_t ell,
            const mpz_t smexp, const mpz_t ell2, const mpz_t invl2)
{
  mpz_poly_power_mod_f_mod_mpz_Barrett(SM, num, F, smexp, ell2, invl2);
  mpz_poly_sub_ui(SM, 1);

  for(int j = 0; j <= SM->deg; j++)
  {
    ASSERT_ALWAYS(mpz_divisible_p(SM->coeff[j], ell));
    mpz_divexact(SM->coeff[j], SM->coeff[j], ell);
    ASSERT_ALWAYS(mpz_cmp(ell, SM->coeff[j])>0);
  }
}

/* Assume nSM > 0 */
void
print_sm (FILE *f, mpz_poly_t SM, int nSM, int d)
{
  d--;
  if (d > SM->deg)
    fprintf(f, "0");
  else
    gmp_fprintf(f, "%Zu", SM->coeff[d]);
  for(int j = 1; j < nSM; j++)
  {
    d--;
    if (d > SM->deg)
      fprintf(f, " 0");
    else
      gmp_fprintf(f, " %Zu", SM->coeff[d]);
  }
  fprintf(f, "\n");
}

/* Computed the Shirokauer maps of a single pair (a,b). 
   SM must be allocated and is viewed as a polynomial of degree F->deg. */
void
sm_single_rel(mpz_poly_t SM, int64_t a, uint64_t b, mpz_poly_t F,
              const mpz_t eps, const mpz_t ell, const mpz_t ell2,
              const mpz_t invl2)
{
  mpz_poly_t rel;

  mpz_poly_init_set_ab(rel, a, b);
  compute_sm (SM, rel, F, ell, eps, ell2, invl2);
  mpz_poly_clear(rel);
}

void
sm_relset_init (sm_relset_t r, int d)
{
  mpz_poly_init (r->num, d);
  mpz_poly_init (r->denom, d);
}

void
sm_relset_clear (sm_relset_t r)
{
  mpz_poly_clear (r->num);
  mpz_poly_clear (r->denom);
}

/* Given an array of index of rows and an array abpolys, construct the polynomial
 * that correspond to the relation-set, i.e. 
 *          rel = prod(abpoly[rk]^ek)
 *       where rk is the index of a row, ek its exponent 
 *             and abpoly[rk] = bk*x - ak
 * rel is build as a fraction (sm_relset_t) and should be init
 */
void
sm_build_one_relset (sm_relset_ptr rel, uint64_t *r, int64_t *e, int len,
                 mpz_poly_t * abpolys, mpz_poly_t F, const mpz_t ell2)
{
  mpz_t ee;
  mpz_init(ee);  
  mpz_poly_t tmp;
  mpz_poly_init(tmp, F->deg);

  /* Set the initial fraction to 1 / 1 */
  mpz_poly_set_zero (rel->num);
  mpz_poly_setcoeff_si(rel->num, 0, 1); 
  mpz_poly_set_zero (rel->denom);
  mpz_poly_setcoeff_si(rel->denom, 0, 1); 

  for(int k = 0 ; k < len ; k++)
  {
    /* Should never happen! */
    ASSERT_ALWAYS(e[k] != 0);

    if (e[k] > 0)
    {
      mpz_set_si(ee, e[k]);
      /* TODO: mpz_poly_long_power_mod_f_mod_mpz */
      mpz_poly_power_mod_f_mod_mpz (tmp, abpolys[r[k]], F, ee, ell2);
      mpz_poly_mul_mod_f_mod_mpz (rel->num, rel->num, tmp, F, ell2, NULL);
    }
    else
    {
      mpz_set_si(ee, -e[k]);
      /* TODO: mpz_poly_long_power_mod_f_mod_mpz */
      mpz_poly_power_mod_f_mod_mpz(tmp, abpolys[r[k]], F, ee, ell2);
      mpz_poly_mul_mod_f_mod_mpz(rel->denom, rel->denom, tmp, F, ell2, NULL);
    }
  }
  mpz_poly_cleandeg(rel->num, F->deg);
  mpz_poly_cleandeg(rel->denom, F->deg);

  mpz_clear (ee);
  mpz_poly_clear(tmp);
}
