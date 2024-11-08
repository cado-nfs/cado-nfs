#include "cado.h" // IWYU pragma: keep

#include <cstdio>      // FILE
#include <cstdlib>     // free, malloc, abort
#include <cstring>     // strcmp memset

#include <vector>

#include <gmp.h>

#include "gmp_aux.h"
#include "macros.h"     // for ASSERT_ALWAYS
#include "mpz_poly.h"   // mpz_poly
#include "sm_utils.hpp"
#include "cado_poly.h"  // NB_POLYS_MAX

/* compute the SM */
static void
compute_sm_lowlevel (cxx_mpz_poly & SM,
        cxx_mpz_poly const & num, cxx_mpz_poly const & F,
        cxx_mpz const & ell, cxx_mpz const & smexp, cxx_mpz const & ell2)
{
    mpz_poly_pow_mod_f_mod_mpz (SM, num, F, smexp, ell2);
    if (SM->deg == -1) {
        throw std::runtime_error(fmt::format(
                "Fatal error: {} has non-zero valuation at the ideal ({}, {}). The code is currently unable to deal with this correctly", num, ell, F));
    }

    mpz_poly_sub_ui (SM, SM, 1);
    mpz_poly_divexact_mpz (SM, SM, ell);
}

double m_seconds = 0;

#define xxxxWIP_SM

/* u and dst may be equal */
void sm_side_info::compute_piecewise(cxx_mpz_poly & dst, cxx_mpz_poly const & u) const
{
    /* now for a split-chunk compute_sm */

    if (nsm == 0)
        return;

    int const n = f0->deg;
    mpz_poly_realloc(dst, n);

    cxx_mpz_poly temp;

    std::vector<cxx_mpz_poly> chunks(pieces.size());

    int s = 0;
    for(size_t j = 0; j < pieces.size() ; j++) {
        if (!pieces[j].is_used)
            continue;
        /* simply call the usual function here */
        mpz_poly_srcptr g = pieces[j].g;
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
                    g, ell, exponents[j], ell2);
        } else {
            ASSERT_ALWAYS(mpz_cmp_ui(mpz_poly_coeff_const(g, 1), 1) == 0);
            mpz_poly_set(chunks[j], u);
            mpz_poly_mod_f_mod_mpz(chunks[j], g, ell2, NULL, NULL);
            ASSERT_ALWAYS(chunks[j]->deg == 0);
            mpz_ptr c = mpz_poly_coeff(chunks[j], 0);
            mpz_powm(c, c, exponents[j], ell2);
            mpz_sub_ui(c, c, 1);
            mpz_divexact(c, c, ell);
        }
#else
        compute_sm_lowlevel (chunks[j], u,
                g, ell, pieces[j].exponent, ell2);
#endif
        for(int k = 0 ; k < g->deg ; k++, s++) {
            if (mode == SM_MODE_2019REV)
                mpz_swap(mpz_poly_coeff(temp, s), mpz_poly_coeff(chunks[j], g->deg-1-k));
            else
                mpz_swap(mpz_poly_coeff(temp, s), mpz_poly_coeff(chunks[j], k));
        }
    }

    if (mode == SM_MODE_LEGACY_PRE2018) {
        temp->deg = n - 1;
        /* now apply the change of basis matrix */
        for(int s = 0 ; s < n ; s++) {
            mpz_ptr y = mpz_poly_coeff(dst, s);
            mpz_set_ui(y, 0);
            for(int k = 0 ; k < n ; k++) {
                mpz_addmul(y,
                        mpz_poly_coeff_const(temp, k),
                        matrix(k, s));
            }
            mpz_mod(y, y, ell);
        }
        dst->deg = n - 1;
    } else {
        temp->deg = s - 1;
        mpz_poly_set(dst, temp);
    }
}

void
print_sm2 (FILE *f, sm_side_info const & S, mpz_poly_srcptr SM, const char * delim)
{
    if (S.nsm == 0)
        return;

    for (int j = 0; j < S.nsm; ++j) {

        int const jx = S.mode == SM_MODE_LEGACY_PRE2018 ? S.f->deg-1-j : j;

        if (jx > SM->deg)
            fprintf(f, "0");
        else
            gmp_fprintf(f, "%Zu", mpz_poly_coeff_const(SM, jx));

        if (j != S.nsm-1)
            fputs(delim, f);
    }
}

void
print_sm (FILE *f, sm_side_info const & S, mpz_poly_srcptr SM)
{
    print_sm2(f, S, SM, " ");
}


void
sm_relset_init (sm_relset_t r, const mpz_poly_srcptr * F, int nb_polys)
{
  r->nb_polys = nb_polys;
  r->num = (mpz_poly *) malloc(nb_polys * sizeof(mpz_poly));
  r->denom = (mpz_poly *) malloc(nb_polys * sizeof(mpz_poly));
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
  r->num = (mpz_poly *) malloc(s->nb_polys * sizeof(mpz_poly));
  r->denom = (mpz_poly *) malloc(s->nb_polys * sizeof(mpz_poly));
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
                    std::vector<pair_and_sides> const & ps,
                    const mpz_poly_srcptr * F, int nb_polys,
		    mpz_srcptr ell2)
{
  mpz_t ee;
  mpz_init(ee);  
  mpz_poly * tmp = (mpz_poly *) malloc(nb_polys * sizeof(mpz_poly));
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
    cxx_mpz_poly const & ab = ps[r[k]].ab;
    const unsigned int * si = ps[r[k]].active_sides;

    if (e[k] > 0)
    {
      mpz_set_si(ee, e[k]);
      /* TODO: mpz_poly_long_power_mod_f_mod_mpz */
      for (int i = 0 ; i < 2 ; i++) {
          unsigned int const s = si[i];
          if (F[s] == NULL) continue;
          mpz_poly_pow_mod_f_mod_mpz(tmp[s], ab, F[s], ee, ell2);
          mpz_poly_mul_mod_f_mod_mpz(rel->num[s], rel->num[s], tmp[s],
                  F[s], ell2, NULL, NULL);
      }
    } else {
      mpz_set_si(ee, -e[k]);
      /* TODO: mpz_poly_long_power_mod_f_mod_mpz */
      for (int i = 0 ; i < 2 ; i++) {
          unsigned int const s = si[i];
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
    int const r1 = mpz_poly_number_of_real_roots(f);
    ASSERT_ALWAYS((f->deg - r1) % 2 == 0);
    int const r2 = (f->deg - r1) / 2;
    int const unitrank = r1 + r2 - 1;
    return unitrank;
}

void compute_change_of_basis_matrix(cxx_mpz_mat & matrix, mpz_poly_srcptr f, std::vector<sm_side_info::piece> const & pieces, cxx_mpz const & ell)
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

    int s = 0;
    for(auto const & pf : pieces) {
        cxx_mpz_poly d, a, b, h;
        mpz_poly_srcptr g = pf.g;

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
                    mpz_set(mpz_mat_entry(matrix, s, k), mpz_poly_coeff_const(h, k));
            }
            mpz_poly_mul_xi(h, h, 1);
            mpz_poly_mod_f_mod_mpz(h, f, ell, NULL, NULL);
        }
    }
}

void sm_side_info::print(FILE * out) const
{
    if (unit_rank == 0) {
        fprintf(out, "# unit rank is 0, no SMs to compute\n");
        return;
    }
    fprintf(out, "# unit rank is %d\n", unit_rank);
    fprintf(out, "# lifted factors of f modulo ell^2\n");
    for(size_t i = 0 ; i < pieces.size() ; i++) {
        gmp_fprintf(out, "# factor %d (used=%d), exponent ell^%d-1=%Zd:\n# ",
                i, pieces[i].is_used,
                pieces[i].g->deg,
                (mpz_srcptr) pieces[i].exponent);
        mpz_poly_fprintf(out, pieces[i].g);
    }
    if (mode == SM_MODE_LEGACY_PRE2018) {
        fprintf(out, "# change of basis matrix to OK/ell*OK from piecewise representation\n");
        for(int i = 0 ; i < f->deg ; i++) {
            fprintf(out, "# ");
            for(int j = 0 ; j < f->deg ; j++) {
                gmp_fprintf(out, " %Zd", matrix[i * f->deg + j]);
            }
            fprintf(out, "\n");
        }
    }
}


sm_side_info::sm_side_info(mpz_poly_srcptr f0, mpz_srcptr ell, bool handle_small_ell)
    : unit_rank (compute_unit_rank(f0))
    , nsm(unit_rank) /* By default, we compute 'unit_rank' SMs. It can
                        be modify by the user. */
{
    if (unit_rank == 0)
        return;

    mpz_set(this->ell, ell);
    mpz_mul(ell2, ell, ell);
    mpz_mul(ell3, ell2, ell);

    mpz_poly_set(this->f0, f0);

    mpz_poly_makemonic_mod_mpz(f, f0, ell2);

    /* note that sm->exponents is initialized at the very end of this
     * function */

    /* compute the lifted factors of f modulo ell^2 */
    {
        /* polynomial factorization is Las Vegas type */
        cxx_gmp_randstate rstate;
        auto fac = mpz_poly_factor_and_lift_padically(f, ell, 2 + handle_small_ell, rstate);
        for(auto const & fm : fac) {
            pieces.emplace_back(fm.first, fm.second, false, 1);
        }
    }


    /* select a subset of factors such that the sum of the degrees is
     * at least the number of sm to compute */
    {
        int s = 0;
        int i = 0;
        while (s < nsm) {
            ASSERT_ALWAYS(pieces.size() > (size_t) i);
            s += pieces[i].g->deg;
            pieces[i].is_used = true;
            i++;
        }
        i--;
        while (s > nsm && i >= 0) {
            int const di = pieces[i].g->deg;
            if (s - di >= nsm) {
                s -= di;
                pieces[i].is_used = false;
            }
            i--;
        }
    }

    /* also compute the lcm of the ell^i-1 ; note that it's only used in
     * the legacy case, afaik. */
    exponent = 1;

    for(auto & P : pieces) {
        mpz_pow_ui(P.exponent, ell, P.g->deg);
        mpz_sub_ui(P.exponent, P.exponent, 1);
        mpz_lcm(exponent, exponent, P.exponent);
    }

    /* We need to compute an ell-maximal order. Well, most probably
     * Z[\alpha\hat] is an ell-maximal order, of course, so it's very
     * probably going to be an easy computation
     */
    // we need to be cxx to do that.
    // cxx_mpq_mat p_maximal_order(cxx_mpz_poly const& f, cxx_mpz const& p)
}

void sm_side_info::set_mode(const char * mode_string)
{
    if (mode_string == NULL || strcmp(mode_string, "default") == 0) {
        mode = SM_MODE_2018;
    } else if (strcmp(mode_string, "2018") == 0) {
        mode = SM_MODE_2018;
    } else if (strcmp(mode_string, "legacy") == 0) {
        mode = SM_MODE_LEGACY_PRE2018;
    } else if (strcmp(mode_string, "2019rev") == 0) {
        mode = SM_MODE_2019REV;
    } else if (strcmp(mode_string, "pseudorandom-combination") == 0) {
        mode = SM_MODE_PSEUDORANDOM_COMBINATION;
        abort(); /* not implemented yet */
    }

    if (mode == SM_MODE_LEGACY_PRE2018) {
        ASSERT_ALWAYS(!matrix);
        matrix = cxx_mpz_mat(f->deg, f->deg);
        compute_change_of_basis_matrix(matrix, f, pieces, ell);
        for(auto & P : pieces)
            P.is_used = true;
    }
}
