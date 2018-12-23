#include "cado.h"
#include <stdio.h>
#include "portability.h"
#include "utils.h"

/* compute the SM */
static void
compute_sm_lowlevel (mpz_poly SM, mpz_poly_srcptr num, const mpz_poly F,
                     const mpz_t ell, const mpz_t smexp, const mpz_t ell2)
{
    mpz_poly_pow_mod_f_mod_mpz (SM, num, F, smexp, ell2);
    mpz_poly_sub_ui (SM, SM, 1);
    mpz_poly_divexact_mpz (SM, SM, ell);
}

double m_seconds = 0;

/* u and dst may be equal */
void compute_sm_piecewise(mpz_poly_ptr dst, mpz_poly_srcptr u, sm_side_info_srcptr sm)
{
    /* now for a split-chunk compute_sm */

    if (sm->nsm == 0)
        return;

    int n = sm->f0->deg;
    mpz_poly_realloc(dst, n);
    mpz_poly_factor_list_srcptr fac = sm->fac;

    mpz_poly temp;
    mpz_poly_init(temp, n);

    /* we can afford calling malloc() here */
    mpz_poly * chunks = malloc(fac->size * sizeof(mpz_poly));
    for(int j = 0 ; j < fac->size ; j++) {
        mpz_poly_srcptr g = fac->factors[j]->f;
        mpz_poly_init(chunks[j], g->deg-1);
    }

    int s = 0;
    for(int j = 0; j < fac->size ; j++) {
        if (!sm->is_factor_used[j])
            continue;
        /* simply call the usual function here */
        mpz_poly_srcptr g = fac->factors[j]->f;
        /* it is tempting to reduce u mod g ; depending on the context,
         * it may even be a good idea. However, doing so automatically is
         * wrong, since this function may also be called with tiny data
         * such as a-bx ; there, multiplying by a-bx with word-size a and
         * b is better than multiplying by a-bx mod g, which may be a
         * much longer integer (in the case deg(g)=1 ; otherwise reducing
         * is a noop in such a case).
         */
#if 0
        if (g->deg > 1) {
            compute_sm_lowlevel (chunks[j], u,
                    g, sm->ell, sm->exponents[j], sm->ell2);
        } else {
            ASSERT_ALWAYS(mpz_cmp_ui(g->coeff[1], 1) == 0);
            mpz_poly_set(chunks[j], u);
            mpz_poly_mod_f_mod_mpz(chunks[j], g, sm->ell2, NULL);
            ASSERT_ALWAYS(chunks[j]->deg == 0);
            mpz_ptr c = chunks[j]->coeff[0];
            mpz_powm(c, c, sm->exponents[j], sm->ell2);
            mpz_sub_ui(c, c, 1);
            mpz_divexact(c, c, sm->ell);
        }
#else
        compute_sm_lowlevel (chunks[j], u,
                g, sm->ell, sm->exponents[j], sm->ell2);
#endif
        for(int k = 0 ; k < g->deg ; k++, s++) {
            mpz_swap(temp->coeff[s], chunks[j]->coeff[k]);
        }
    }

    temp->deg = s - 1;

    mpz_poly_set(dst, temp);

    for(int j = 0 ; j < fac->size ; j++) {
        mpz_poly_clear(chunks[j]);
    }
    free(chunks);
    mpz_poly_clear(temp);
}

/* Assume nSM > 0 */
void
print_sm2 (FILE *f, mpz_poly SM, int nSM, MAYBE_UNUSED int d, const char * delim)
{
  if (nSM == 0)
      return;
  for (int j = 0; j < nSM; ++j) {
      if (j > SM->deg)
          fprintf(f, "0");
      else
          gmp_fprintf(f, "%Zu", SM->coeff[j]);
      if (j != nSM-1)
          fprintf(f, "%s", delim);
  }
}

void
print_sm (FILE *f, mpz_poly SM, int nSM, int d)
{
    print_sm2(f, SM, nSM, d, " ");
}


void
sm_relset_init (sm_relset_t r, int *d, int nb_polys)
{
  r->nb_polys = nb_polys;
  for (int side = 0; side < nb_polys; side++) {
    mpz_poly_init (r->num[side], d[side]);
    mpz_poly_init (r->denom[side], d[side]);
  }
}

void
sm_relset_clear (sm_relset_t r, int nb_polys)
{
  for (int side = 0; side < nb_polys; side++) {
    mpz_poly_clear (r->num[side]);
    mpz_poly_clear (r->denom[side]);
  }
}

void
sm_relset_copy (sm_relset_t r, sm_relset_srcptr s)
{
  r->nb_polys = s->nb_polys;
  for (int side = 0; side < r->nb_polys; side++) {
    mpz_poly_set (r->num[side], s->num[side]);
    mpz_poly_set (r->denom[side], s->denom[side]);
  }
}


/* Given an array of index of rows and an array of abpolys,
 * construct the polynomial that corresponds to the relation-set, i.e. 
 *          rel = prod(abpoly[rk]^ek)
 * where rk is the index of a row, ek its exponent 
 * and abpoly[rk] = bk*x - ak.
 * rel is built as a fraction (sm_relset_t) and should be initialized.
 */
void
sm_build_one_relset (sm_relset_ptr rel, uint64_t *r, int64_t *e, int len,
		     mpz_poly * abpolys, mpz_poly_ptr *F, int nb_polys,
		     const mpz_t ell2)
{
  mpz_t ee;
  mpz_init(ee);  
  mpz_poly tmp[NB_POLYS_MAX];
  for (int side = 0; side < nb_polys; side++) {
    if (F[side] == NULL) continue;
    mpz_poly_init(tmp[side], F[side]->deg);
  }

  /* Set the initial fraction to 1 / 1 */
  for (int side = 0; side < nb_polys; side++) {
    mpz_poly_set_zero (rel->num[side]);
    mpz_poly_setcoeff_si(rel->num[side], 0, 1); 
    mpz_poly_set_zero (rel->denom[side]);
    mpz_poly_setcoeff_si(rel->denom[side], 0, 1); 
  }

  for(int k = 0 ; k < len ; k++)
  {
    /* Should never happen! */
    ASSERT_ALWAYS(e[k] != 0);

    if (e[k] > 0)
    {
      mpz_set_si(ee, e[k]);
      /* TODO: mpz_poly_long_power_mod_f_mod_mpz */
      for (int s = 0; s < nb_polys; ++s) {
        if (F[s] == NULL) continue;
        mpz_poly_pow_mod_f_mod_mpz (tmp[s], abpolys[r[k]], F[s], ee, ell2);
        mpz_poly_mul_mod_f_mod_mpz (rel->num[s], rel->num[s], tmp[s], F[s],
                                    ell2, NULL);
      }
    }
    else
    {
      mpz_set_si(ee, -e[k]);
      /* TODO: mpz_poly_long_power_mod_f_mod_mpz */
      for (int s = 0; s < nb_polys; ++s) {
        if (F[s] == NULL) continue;
        mpz_poly_pow_mod_f_mod_mpz(tmp[s], abpolys[r[k]], F[s], ee, ell2);
        mpz_poly_mul_mod_f_mod_mpz(rel->denom[s], rel->denom[s], tmp[s], F[s],
                                   ell2, NULL);
      }
    }
  }
  for (int s = 0; s < nb_polys; ++s) {
    if (F[s] == NULL) continue;
    mpz_poly_cleandeg(rel->num[s], F[s]->deg);
    mpz_poly_cleandeg(rel->denom[s], F[s]->deg);
    mpz_poly_clear(tmp[s]);
  }

  mpz_clear (ee);
}

static int compute_unit_rank(mpz_poly_srcptr f)
{
    int r1 = mpz_poly_number_of_real_roots(f);
    ASSERT_ALWAYS((f->deg - r1) % 2 == 0);
    int r2 = (f->deg - r1) / 2;
    int unitrank = r1 + r2 - 1;
    return unitrank;
}

void sm_side_info_print(FILE * out, sm_side_info_srcptr sm)
{
    if (sm->unit_rank == 0) {
        fprintf(out, "# unit rank is 0, no SMs to compute\n");
        return;
    }
    fprintf(out, "# unit rank is %d\n", sm->unit_rank);
    fprintf(out, "# lifted factors of f modulo ell^2\n");
    for(int i = 0 ; i < sm->fac->size ; i++) {
        gmp_fprintf(out, "# factor %d (used=%d), exponent ell^%d-1=%Zd:\n# ",
                i, sm->is_factor_used[i], sm->fac->factors[i]->f->deg,
                sm->exponents[i]);
        mpz_poly_fprintf(out, sm->fac->factors[i]->f);
    }
}


void sm_side_info_init(sm_side_info_ptr sm, mpz_poly_srcptr f0, mpz_srcptr ell)
{
    memset(sm, 0, sizeof(*sm));

    sm->unit_rank = compute_unit_rank(f0);
    sm->nsm = sm->unit_rank; /* By default, we compute 'unit_rank' SMs. It can
                                be modify by the user. */

    if (sm->unit_rank == 0)
        return;

    /* initialize all fields */
    mpz_init_set(sm->ell, ell);
    mpz_init(sm->ell2);
    mpz_mul(sm->ell2, sm->ell, sm->ell);

    mpz_poly_init(sm->f, -1);
    sm->f0 = f0;
    mpz_poly_makemonic_mod_mpz(sm->f, f0, sm->ell2);
    mpz_poly_factor_list_init(sm->fac);

    mpz_init(sm->exponent);
    /* note that sm->exponents is initialized at the very end of this
     * function */

    /* compute the lifted factors of f modulo ell^2 */
    {
        /* polynomial factorization is Las Vegas type */
        gmp_randstate_t rstate;
        gmp_randinit_default(rstate);
        mpz_poly_factor_and_lift_padically(sm->fac, sm->f, sm->ell, 2, rstate);
        gmp_randclear(rstate);
    }

    /* select a subset of factors such that the sum of the degrees is
     * at least the number of sm to compute */
    {
        int s = 0;
        int i = 0;
        while (s < sm->nsm) {
            ASSERT_ALWAYS(sm->fac->size > i);
            s += sm->fac->factors[i]->f->deg;
            sm->is_factor_used[i] = 1;
            i++;
        }
        i--;
        while (s > sm->nsm && i >= 0) {
            int di = sm->fac->factors[i]->f->deg;
            if (s - di >= sm->nsm) {
                s -= di;
                sm->is_factor_used[i] = 0;
            }
            i--;
        }
    }

    /* also compute the lcm of the ell^i-1 */
    sm->exponents = malloc(sm->fac->size * sizeof(mpz_t));
    mpz_set_ui(sm->exponent, 1);
    for(int i = 0 ; i < sm->fac->size ; i++) {
        mpz_init(sm->exponents[i]);
        mpz_pow_ui(sm->exponents[i], sm->ell, sm->fac->factors[i]->f->deg);
        mpz_sub_ui(sm->exponents[i], sm->exponents[i], 1);
        mpz_lcm(sm->exponent, sm->exponent, sm->exponents[i]);
    }
}

void sm_side_info_clear(sm_side_info_ptr sm)
{
    if (sm->unit_rank == 0)
        return;

    mpz_clear(sm->exponent);

    for(int i = 0 ; i < sm->fac->size ; i++) {
        mpz_clear(sm->exponents[i]);
    }
    free(sm->exponents);

    mpz_poly_factor_list_clear(sm->fac);
    mpz_poly_clear(sm->f);

    mpz_clear(sm->ell2);
    mpz_clear(sm->ell);
}

