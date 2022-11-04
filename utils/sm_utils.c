#include "cado.h" // IWYU pragma: keep
#include <stdio.h>      // FILE
#include <stdlib.h>     // free, malloc, abort
#include <string.h>     // strcmp memset
#include <gmp.h>
#include "macros.h"     // for ASSERT_ALWAYS
#include "mpz_poly.h"   // mpz_poly
#include "sm_utils.h"
#include "cado_poly.h"  // NB_POLYS_MAX

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
            mpz_poly_mod_f_mod_mpz(chunks[j], g, sm->ell2, NULL, NULL);
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
            if (sm->mode == SM_MODE_2019REV)
                mpz_swap(temp->coeff[s], chunks[j]->coeff[g->deg-1-k]);
            else
                mpz_swap(temp->coeff[s], chunks[j]->coeff[k]);
        }
    }

    if (sm->mode == SM_MODE_LEGACY_PRE2018) {
        temp->deg = n - 1;
        /* now apply the change of basis matrix */
        for(int s = 0 ; s < n ; s++) {
            mpz_set_ui(dst->coeff[s], 0);
            for(int k = 0 ; k < n ; k++) {
                mpz_addmul(dst->coeff[s],
                        temp->coeff[k],
                        sm->matrix[k * n + s]);
            }
            mpz_mod(dst->coeff[s], dst->coeff[s], sm->ell);
        }
        dst->deg = n - 1;
    } else {
        temp->deg = s - 1;
        mpz_poly_set(dst, temp);
    }

    for(int j = 0 ; j < fac->size ; j++) {
        mpz_poly_clear(chunks[j]);
    }
    free(chunks);
    mpz_poly_clear(temp);
}

void
print_sm2 (FILE *f, sm_side_info_srcptr S, mpz_poly_srcptr SM, const char * delim)
{
    if (S->nsm == 0)
        return;

    for (int j = 0; j < S->nsm; ++j) {

        int jx = S->mode == SM_MODE_LEGACY_PRE2018 ? S->f->deg-1-j : j;

        if (jx > SM->deg)
            fprintf(f, "0");
        else
            gmp_fprintf(f, "%Zu", SM->coeff[jx]);

        if (j != S->nsm-1)
            fputs(delim, f);
    }
}

void
print_sm (FILE *f, sm_side_info_srcptr S, mpz_poly_srcptr SM)
{
    print_sm2(f, S, SM, " ");
}


void
sm_relset_init (sm_relset_t r, const mpz_poly_srcptr * F, int nb_polys)
{
  r->nb_polys = nb_polys;
  r->num = malloc(nb_polys * sizeof(mpz_poly));
  r->denom = malloc(nb_polys * sizeof(mpz_poly));
  for (int side = 0; side < nb_polys; side++) {
    mpz_poly_init (r->num[side], F[side] ? mpz_poly_degree(F[side]) : -1);
    mpz_poly_init (r->denom[side],  F[side] ? mpz_poly_degree(F[side]) : -1);
  }
}

void
sm_relset_clear (sm_relset_t r)
{
  for (int side = 0; side < r->nb_polys; side++) {
    mpz_poly_clear (r->num[side]);
    mpz_poly_clear (r->denom[side]);
  }
  free(r->num);
  free(r->denom);
}

void
sm_relset_copy (sm_relset_t r, sm_relset_srcptr s)
{
  if (r == s)
      return;
  sm_relset_clear(r);
  r->nb_polys = s->nb_polys;
  r->num = malloc(s->nb_polys * sizeof(mpz_poly));
  r->denom = malloc(s->nb_polys * sizeof(mpz_poly));
  for (int side = 0; side < r->nb_polys; side++) {
    mpz_poly_init (r->num[side], -1);
    mpz_poly_init (r->denom[side], -1);
    mpz_poly_set (r->num[side], s->num[side]);
    mpz_poly_set (r->denom[side], s->denom[side]);
  }
}


/* Given an array of row indices and an array of abpolys,
 * construct the polynomial that corresponds to the relation-set, i.e. 
 *          rel = prod(abpoly[rk]^ek)
 * where rk is the index of a row, ek its exponent 
 * and abpoly[rk] = bk*x - ak.
 * rel is built as a fraction (sm_relset_t) and should be initialized.
 */
void
sm_build_one_relset(sm_relset_ptr rel,
                    const uint64_t *r, const int64_t *e, int len,
                    const pair_and_sides * ps,
                    const mpz_poly_srcptr * F, int nb_polys,
		    mpz_srcptr ell2)
{
  mpz_t ee;
  mpz_init(ee);  
  mpz_poly * tmp = malloc(nb_polys * sizeof(mpz_poly));
  memset(tmp, 0, nb_polys * sizeof(mpz_poly));
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
    mpz_poly_srcptr ab = ps[r[k]]->ab;
    const unsigned int * si = ps[r[k]]->active_sides;

    if (e[k] > 0)
    {
      mpz_set_si(ee, e[k]);
      /* TODO: mpz_poly_long_power_mod_f_mod_mpz */
      for (int i = 0 ; i < 2 ; i++) {
          unsigned int s = si[i];
          if (F[s] == NULL) continue;
          mpz_poly_pow_mod_f_mod_mpz(tmp[s], ab, F[s], ee, ell2);
          mpz_poly_mul_mod_f_mod_mpz(rel->num[s], rel->num[s], tmp[s],
                  F[s], ell2, NULL, NULL);
      }
    } else {
      mpz_set_si(ee, -e[k]);
      /* TODO: mpz_poly_long_power_mod_f_mod_mpz */
      for (int i = 0 ; i < 2 ; i++) {
          unsigned int s = si[i];
          if (F[s] == NULL) continue;
          mpz_poly_pow_mod_f_mod_mpz(tmp[s], ab, F[s], ee, ell2);
          mpz_poly_mul_mod_f_mod_mpz(rel->denom[s], rel->denom[s], tmp[s], F[s],
                  ell2, NULL, NULL);
      }
    }
  }
  for (int s = 0; s < nb_polys; ++s) {
    if (F[s] == NULL) continue;
    mpz_poly_cleandeg(rel->num[s], F[s]->deg);
    mpz_poly_cleandeg(rel->denom[s], F[s]->deg);
    mpz_poly_clear(tmp[s]);
  }
  free(tmp);
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

void compute_change_of_basis_matrix(mpz_t * matrix, mpz_poly_srcptr f, mpz_poly_factor_list_srcptr fac, mpz_srcptr ell)
{
    /* now compute the change of basis matrix. This is given simply
     * as the CRT matrix from the piecewise representation modulo
     * the factors of f to the full representation.
     *
     * We thus have full_repr(x) = \sum multiplier_g * repr_g(x)
     *
     * Where multiplier_g is a multiple of f/g, such that
     * f/g*multiplier_g is one modulo f.
     *
     * Note that the computation of this matrix is ok to do mod ell
     * only ; going to ell^2 is unnecessary.
     */

    for(int i = 0, s = 0 ; i < fac->size ; i++) {
        mpz_poly d, a, b, h;
        mpz_poly_srcptr g = fac->factors[i]->f;
        mpz_poly_init(h, f->deg - g->deg);
        mpz_poly_init(d, 0);
        mpz_poly_init(a, f->deg - g->deg - 1);
        mpz_poly_init(b, g->deg - 1);

        /* compute h = product of other factors */
        mpz_poly_divexact(h, f, g, ell);

        /* now invert h modulo g */

        /* say a*g + b*h = 1 */
        mpz_poly_xgcd_mpz(d, g, h, a, b, ell);

        /* we now have the complete cofactor */
        mpz_poly_mul_mod_f_mod_mpz(h, b, h, f, ell, NULL, NULL);
        for(int j = 0 ; j < g->deg ; j++, s++) {
            /* store into the matrix the coefficients of x^j*h
             * modulo f */
            for(int k = 0 ; k < f->deg ; k++) {
                if (k <= h->deg)
                    mpz_set(matrix[s * f->deg + k], h->coeff[k]);
            }
            mpz_poly_mul_xi(h, h, 1);
            mpz_poly_mod_f_mod_mpz(h, f, ell, NULL, NULL);
        }
        mpz_poly_clear(b);
        mpz_poly_clear(a);
        mpz_poly_clear(d);
        mpz_poly_clear(h);
    }
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
                i, sm->is_factor_used[i],
                sm->fac->factors[i]->f->deg,
                sm->exponents[i]);
        mpz_poly_fprintf(out, sm->fac->factors[i]->f);
    }
    if (sm->mode == SM_MODE_LEGACY_PRE2018) {
        fprintf(out, "# change of basis matrix to OK/ell*OK from piecewise representation\n");
        for(int i = 0 ; i < sm->f->deg ; i++) {
            fprintf(out, "# ");
            for(int j = 0 ; j < sm->f->deg ; j++) {
                gmp_fprintf(out, " %Zd", sm->matrix[i * sm->f->deg + j]);
            }
            fprintf(out, "\n");
        }
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
    sm->exponents = (mpz_t *) malloc(sm->fac->size * sizeof(mpz_t));
    mpz_set_ui(sm->exponent, 1);
    for(int i = 0 ; i < sm->fac->size ; i++) {
        mpz_init(sm->exponents[i]);
        mpz_pow_ui(sm->exponents[i], sm->ell, sm->fac->factors[i]->f->deg);
        mpz_sub_ui(sm->exponents[i], sm->exponents[i], 1);
        mpz_lcm(sm->exponent, sm->exponent, sm->exponents[i]);
    }
}

void sm_side_info_set_mode(sm_side_info_ptr sm, const char * mode_string)
{
    if (mode_string == NULL || strcmp(mode_string, "default") == 0) {
        sm->mode = SM_MODE_2018;
    } else if (strcmp(mode_string, "2018") == 0) {
        sm->mode = SM_MODE_2018;
    } else if (strcmp(mode_string, "legacy") == 0) {
        sm->mode = SM_MODE_LEGACY_PRE2018;
    } else if (strcmp(mode_string, "2019rev") == 0) {
        sm->mode = SM_MODE_2019REV;
    } else if (strcmp(mode_string, "pseudorandom-combination") == 0) {
        sm->mode = SM_MODE_PSEUDORANDOM_COMBINATION;
        abort(); /* not implemented yet */
    }

    if (sm->mode == SM_MODE_LEGACY_PRE2018) {
        sm->matrix = (mpz_t *) malloc(sm->f->deg * sm->f->deg * sizeof(mpz_t));
        for(int i = 0 ; i < sm->f->deg ; i++)
            for(int j = 0 ; j < sm->f->deg ; j++)
                mpz_init_set_ui(sm->matrix[i * sm->f->deg + j], 0);
        compute_change_of_basis_matrix(sm->matrix, sm->f, sm->fac, sm->ell);
        for(int j = 0; j < sm->fac->size ; j++)
            sm->is_factor_used[j] = 1;
    }
}


void sm_side_info_clear(sm_side_info_ptr sm)
{
    if (sm->unit_rank == 0)
        return;

    if (sm->matrix) {
        for(int i = 0 ; i < sm->f->deg ; i++)
            for(int j = 0 ; j < sm->f->deg ; j++)
                mpz_clear(sm->matrix[i * sm->f->deg + j]);
        free(sm->matrix);
    }

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

