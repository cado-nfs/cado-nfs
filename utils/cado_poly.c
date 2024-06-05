#include "cado.h" // IWYU pragma: keep
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "macros.h"  // for ASSERT_ALWAYS, ASSERT
#include "params.h"
#include "cado_poly.h"
#include "portability.h"

/* Be conservative and allocate two polynomials by default. */
void cado_poly_init(cado_poly_ptr cpoly)
{
    /* ALL fields are zero upon init, EXCEPT the degree field (which is -1) */
    memset(cpoly, 0, sizeof(cpoly[0]));

    /* By default allocate 2 polynomials */
    cpoly->nb_polys = 0;
    cpoly->pols = NULL;

    mpz_init_set_ui(cpoly->n, 0);
}

void cado_poly_provision_new_poly(cado_poly_ptr cpoly)
{
    cpoly->nb_polys++;
    cpoly->pols = realloc(cpoly->pols, cpoly->nb_polys * sizeof(mpz_poly));
    mpz_poly_init(cpoly->pols[cpoly->nb_polys-1], -1);
}

void cado_poly_reset(cado_poly_ptr cpoly)
{
    for(int side = 0 ; side < cpoly->nb_polys ; side++)
      mpz_poly_clear (cpoly->pols[side]);

    free(cpoly->pols);
    cpoly->nb_polys = 0;
    cpoly->pols = NULL;

    mpz_set_ui(cpoly->n, 0);
}

void cado_poly_clear(cado_poly_ptr cpoly)
{
    for(int side = 0 ; side < cpoly->nb_polys ; side++)
      mpz_poly_clear (cpoly->pols[side]);

    free(cpoly->pols);
    mpz_clear(cpoly->n);

    memset(cpoly, 0, sizeof(cpoly[0]));
}

/* p <- q */
void
cado_poly_set (cado_poly_ptr p, cado_poly_srcptr q)
{
    if (p == q) return;
    mpz_set (p->n, q->n);
    p->skew = q->skew;
    for( ; p->nb_polys > q->nb_polys ; ) {
        mpz_poly_clear(p->pols[--p->nb_polys]);
    }
    for(int side = 0 ; side < q->nb_polys ; side++) {
        if (side == p->nb_polys)
            cado_poly_provision_new_poly(p);
        mpz_poly_set (p->pols[side], q->pols[side]);
    }
}

void
cado_poly_swap (cado_poly_ptr p, cado_poly_ptr q)
{
    if (p == q) return;
    mpz_swap (p->n, q->n);
    { double t = p->skew; p->skew = q->skew; q->skew = t; }
    { unsigned int t = p->nb_polys; p->nb_polys = q->nb_polys; q->nb_polys = t; }
    { mpz_poly * t = p->pols; p->pols = q->pols; q->pols = t; }
}

// This function is no longer exported
#define BUF_MAX 10000

// returns 0 on failure, 1 on success.
int cado_poly_set_plist(cado_poly_ptr cpoly, param_list_ptr pl)
{
  int ret = 1;

  /* Parse skew value. Set to 0.0 to ensure that we get an invalid skewness
     in case it is not given */
  cpoly->skew = 0.0;
  param_list_parse_double(pl, "skew", &(cpoly->skew));

  /* Parse polynomials. Either in line format (poly%d=%Zd,%Zd,...) either given
   * by coefficients. */
  for( ;; ) {
      char tag[15];
      snprintf(tag, sizeof(tag), "poly%d", cpoly->nb_polys);
      if (!param_list_lookup_string(pl, tag))
          break;
      cado_poly_provision_new_poly(cpoly);
      /* sure, we could probably use the result of the string lookup.
       * Well, that's not how param_list_parse_mpz_poly works.
       * Note that param_list_parse_mpz_poly in itself supports two
       * distinct formats: comma-separated lists of coefficients, or
       * algebraic expressions.
       */
      param_list_parse_mpz_poly(pl, tag, cpoly->pols[cpoly->nb_polys-1]);
  }

  if (cpoly->nb_polys == 0) /* No polynomials so far. Try legacy format. */
  /* Convention: c or X means side 1 (often algebraic)
   *             Y means side 0 (often rational) */
  {
    cado_poly_provision_new_poly(cpoly);
    cado_poly_provision_new_poly(cpoly);
    mpz_t coeff;
    mpz_init(coeff);
    /* reading polynomials coefficient by coefficient */
    for (unsigned int i = 0; i <= MAX_DEGREE; i++)
    {
      char tc[4]; snprintf(tc, sizeof(tc), "c%d", i);
      char tX[4]; snprintf(tX, sizeof(tX), "X%d", i);
      char tY[4]; snprintf(tY, sizeof(tY), "Y%d", i);
      if (param_list_parse_mpz(pl, tc, coeff))
          mpz_poly_setcoeff(cpoly->pols[1], i, coeff);
      if (param_list_parse_mpz(pl, tX, coeff))
          mpz_poly_setcoeff(cpoly->pols[1], i, coeff);
      if (param_list_parse_mpz(pl, tY, coeff))
          mpz_poly_setcoeff(cpoly->pols[0], i, coeff);
    }
    mpz_clear(coeff);
  }

  /* Parse value of N. Two keys possible: n or None. Return 0 if not found. */
  if (!param_list_parse_mpz(pl, "n", cpoly->n) &&
      !param_list_parse_mpz(pl, NULL, cpoly->n))
  {
    fprintf (stderr, "Error, no value for N in cado_poly_set_plist\n");
    ret = 0;
  }

  for (int side = 0; side < cpoly->nb_polys; side++)
    if (cpoly->pols[side]->deg < 0)
    {
      fprintf (stderr, "Error, polynomial on side %u has degree < 0 in "
                       "cado_poly_set_plist\n", side);
      ret = 0;
    }

  for (int side = 0; side < cpoly->nb_polys; side++) {
      mpz_poly_ptr g = cpoly->pols[side];
      mpz_t c;
      mpz_init(c);
      mpz_poly_content(c, g);
      if (mpz_cmp_ui(c, 1) != 0) {
          gmp_fprintf (stderr, "# Warning: the polynomial on side %d"
                  " has non-trivial content (gcd of coeffs = %Zd).\n"
                  "# This is a very unexpected situation."
                  " For all useful cases, it seems that\n"
                  "# the correct thing to do is to divide the polynomial"
                  " by its content. Doing it now.\n",
                  side, c);
          fprintf(stderr, "# Old poly: ");
          mpz_poly_fprintf(stderr, g);
          for(int j = 0 ; j <= g->deg ; j++) {
              mpz_fdiv_q(g->coeff[j],g->coeff[j],c);
          }
          fprintf(stderr, "# Reduced poly: ");
          mpz_poly_fprintf(stderr, g);
      }
      mpz_clear(c);
  }

  /* check that N divides the resultant*/
  ASSERT_ALWAYS(cado_poly_check_mapping(NULL, cpoly, cpoly->n));

  return ret;
}

// returns 0 on failure, 1 on success.
int cado_poly_read_stream (cado_poly_ptr cpoly, FILE * f)
{
  param_list pl;
  param_list_init (pl);
  param_list_read_stream (pl, f, 0);
  int r = cado_poly_set_plist (cpoly, pl);
  param_list_clear (pl);
  return r;
}

// returns 0 on failure, 1 on success.
int cado_poly_read_next_poly_from_stream (cado_poly_ptr cpoly, FILE * f)
{
  int r;
  param_list pl;
  param_list_init (pl);
  r = param_list_read_stream (pl, f, 1);
  if (r && !param_list_empty(pl)) {
    cado_poly_reset(cpoly);
    r = cado_poly_set_plist (cpoly, pl);
  } else {
    r = 0;
  }
  param_list_clear (pl);
  return r;
}

// returns 0 on failure, 1 on success.
int cado_poly_read(cado_poly_ptr cpoly, const char *filename)
{
    FILE *file;
    int r;
    const char * magic = "inline-poly://";
    if (strncmp(filename, magic, strlen(magic)) == 0) {
        /* interpret the reset as param list items, and create the
         * polynomial like this.
         */
        param_list pl;
        param_list_init(pl);
        const char * p = filename + strlen(magic);
        for( ; *p ; ) {
            const char * q = strchr(p, '/');
            char * newitem;
            if (q == NULL) {
                newitem = strdup(p);
                p += strlen(p);
            } else {
                newitem = strndup(p, q - p);
                p = q + 1;
            }
            char * newkey;
            char * newvalue;
            q = strchr(newitem, '=');
            if (q == NULL) {
                fprintf(stderr, "wrong value in inline-poly file\n");
                exit (EXIT_FAILURE);
            }
            newitem[q-newitem] = '\0';
            newkey = newitem;
            newvalue = newitem + (q-newitem+1);
            param_list_add_key(pl, newkey, newvalue, PARAMETER_FROM_FILE);
            free(newitem);
        }
        r = cado_poly_set_plist (cpoly, pl);
        param_list_clear(pl);
    } else {
        file = fopen(filename, "r");
        if (file == NULL)
        {
          fprintf(stderr, "Error, in cado_poly_read: could not open %s\n", filename);
          return 0;
        }
        r = cado_poly_read_stream(cpoly, file);
        fclose(file);
    }
    return r;
}

/* check that there is a mapping from Z[x] to the target ring (Z/N,
 * GF(p), GF(p^k)) that factors through each of the polynomials.
 *
 * If we have only two sides, rational or not, this is clearly something
 * that we want, for NFS to make any sort of sense.
 *
 * The question of how to generalize this check to multiple polynomials
 * is not entirely clear.
 *
 * If there is a rational side (only one, see cado_poly_get_ratside),
 * having a resultant with the rational polynomial that is a multiple of
 * N is actually a transitive relation, so that it suffices to check the
 * condition with the rational side only.
 *
 * If we have multiple nonlinear polynomials, we have to make sure that
 * there is a compatible mapping. This is done by checking the pseudo-gcd
 * of the set of polynomials in cado_poly
 *
 * => The conclusion of all that is that what we really care about is the
 * computation of the gcd, and this gcd is returned in the polynomial G.
 * Note that if G is null, then we don't store it, but we do the check
 * nevertheless.
 *
 * This returns true (1) on success, and false (0) on error, after an
 * error message is printed to stderr.
 */
int cado_poly_check_mapping(mpz_poly_ptr G, cado_poly_srcptr cpoly,
        mpz_srcptr N)
{
    /* scratch space for mpz_poly_pseudogcd_mpz, which destroys its
     * input.
     */
    mpz_poly S[2];
    mpz_poly_init(S[0], -1);
    mpz_poly_init(S[1], -1);

    /* a polynomial for the gcd */
    mpz_poly G_loc;
    mpz_poly_init(G_loc, -1);
    mpz_poly_set_xi(G_loc, -1);

    /* factors of N, if any. This should never happen */
    mpz_t factors_of_N, ftmp;
    mpz_init_set_ui(factors_of_N, 1);
    mpz_init_set_ui(ftmp, 0);

    /* A copy of N with all trivial factors removed */
    mpz_t Nsmall;
    mpz_init(Nsmall);

    for(int ret = 0 ; ret == 0; ) {
        mpz_divexact(Nsmall, N, factors_of_N);

        /* compute all pseudo-gcds */
        for(int i = 1 ; i < cpoly->nb_polys ; i++) {
            mpz_poly_set(S[0], cpoly->pols[0]);
            mpz_poly_set(S[1], cpoly->pols[i]);
            ret = mpz_poly_pseudogcd_mpz (S[0], S[1], Nsmall, ftmp);
            if (ret == 0) {
                mpz_mul(factors_of_N, factors_of_N, ftmp);
                break;
            }
            ret = mpz_poly_pseudogcd_mpz (G_loc, S[0], Nsmall, ftmp);
            if (ret == 0) {
                mpz_mul(factors_of_N, factors_of_N, ftmp);
                break;
            }
        }

        /* If a factor was found, start over, otherwise we're done with
         * this silly "restart if small factors pop up" loop (which will
         * presumably never be activated anyway...)
         */
    }

    int found_mapping = G_loc->deg >= 1;

    if (G)
        mpz_poly_set(G, G_loc);

    if (mpz_cmp_ui(factors_of_N, 1) != 0) {
        /* I don't think that there's any reason to have a situation
         * where N has non-trivial factors. Well, it happens if you feed
         * random numbers to polyselect, though.
         */
        gmp_fprintf(stderr, "Warning: non-trivial factors (%Zd) of N were found while checking the cpoly file. This should not happen.\n", factors_of_N);
        // found_mapping = 0;
    }


    mpz_clear(Nsmall);
    mpz_clear(ftmp);
    mpz_clear(factors_of_N);
    mpz_poly_clear(G_loc);
    mpz_poly_clear(S[0]);
    mpz_poly_clear(S[1]);

    if (!found_mapping) {
        fprintf (stderr, "Error, we found no compatible mapping that factors through the given number fields. Your poly file is most probably wrong.\n");
    }

    return found_mapping;
}


/* If NULL is passed for m, then, just check that N divides the resultant.
 * Assume there at least two polynomials.
 */
int cado_poly_getm(mpz_ptr m, cado_poly_srcptr cpoly, mpz_srcptr N)
{
    mpz_poly G;
    mpz_poly_init(G, -1);

    ASSERT_ALWAYS(cado_poly_check_mapping(G, cpoly, N));

    /* It doesn't make sense to ask for m when the mapping goes to a
     * larger structure. Note that this is not necessarily equivalent to
     * having no rational side: Two algebraic sides may well have only a
     * rational mapping in common (case of Joux-Lercier polynomial
     * selection).
     */

    if (m != NULL && mpz_poly_degree(G) != 1) {
        fprintf(stderr, "Error, cannot determine m given that the common mapping is not linear\n");
        exit (EXIT_FAILURE);
    }

    if (m != NULL) {
        mpz_t inv;
        mpz_init(inv);
        int ret2 = mpz_invert(inv, G->coeff[1], N);
        // This inversion should always work.
        // If not, it means that N has a small factor (not sure we want
        // to be robust against that...). This should have been reported,
        // at least as a warning, by cado_poly_check_mapping.
        // Or maybe the polynomial selection was really bogus!
        ASSERT_ALWAYS(ret2);
        mpz_mul(inv, inv, G->coeff[0]);
        mpz_neg(inv, inv);
        mpz_mod(m, inv, N);
        mpz_clear(inv);
      }

    mpz_poly_clear(G);

    return 1;
}

/* Return the rational side or -1 if only algebraic sides */
/* Assume that there is at most 1 rational side (renumber.cpp does the same
 * assumption).
 *
 * Note that several distinct rational sides would not really make sense.
 * A rational side prescribes the mapping from the indeterminate (the X
 * in a-b*X) to the underlying ring (Z/N, GF(p), ...). So if we have
 * distinct rational sides, we would have distinct mappings, and I don't
 * think that we could do something useful in such a situation.
 */
int
cado_poly_get_ratside (cado_poly_srcptr cpoly)
{
  int number_of_rational_sides = 0;
  int ratside = -1;
  for(int side = 0; side < cpoly->nb_polys; side++)
    if(cpoly->pols[side]->deg == 1) {
        ratside = side;
        number_of_rational_sides++;
    }
  ASSERT_ALWAYS(number_of_rational_sides <= 1);
  return ratside;
}

void
cado_poly_fprintf (FILE *fp, cado_poly_srcptr cpoly, const char *prefix)
{
  if (prefix)
    fputs (prefix, fp);
  gmp_fprintf (fp, "n: %Zd\n", cpoly->n);

  if (cpoly->nb_polys == 2) {
    mpz_poly_fprintf_cado_format (fp, cpoly->pols[0], 'Y', prefix);
    mpz_poly_fprintf_cado_format (fp, cpoly->pols[1], 'c', prefix);
  } else {
    for (int side = 0; side < cpoly->nb_polys; side++) {
      if (prefix)
        fputs (prefix, fp);
      fprintf (fp, "poly%u=", side);
      mpz_poly_fprintf_coeffs (fp, cpoly->pols[side], ',');
    }
  }

  if (prefix)
    fputs (prefix, fp);
  fprintf (fp, "skew: %1.3f\n", cpoly->skew);
}

/* if exp_E = 0, print E = lognorm + alpha (root-optimized polynomial),
   otherwise print exp_E */
void
cado_poly_fprintf_info (FILE *fp, double lognorm, double exp_E, double alpha,
                        double alpha_proj, unsigned int nrroots,
                        const char *prefix)
{
  if (prefix)
    fputs (prefix, fp);
  /* Always print "# " after the prefix and before the info line. */
  fprintf (fp, "# lognorm %1.2f, %s %1.2f, alpha %1.2f (proj %1.2f),"
             " %u real root%s\n",
             lognorm, (exp_E == 0) ? "E" : "exp_E",
             (exp_E == 0) ? lognorm + alpha : exp_E, alpha, alpha_proj,
             nrroots, (nrroots <= 1) ? "" : "s");
}

void
cado_poly_fprintf_MurphyE (FILE *fp, double MurphyE, double bound_f,
                           double bound_g, double area, const char *prefix)
{
  if (prefix)
    fputs (prefix, fp);
  /* Always print "# " after the prefix and before the MurphyE line. */
  fprintf (fp, "# MurphyE(Bf=%.3e,Bg=%.3e,area=%.3e)=%.3e\n", bound_f, bound_g,
               area, MurphyE);
}
