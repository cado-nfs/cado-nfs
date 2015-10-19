#include "cado.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "params.h"
#include "cado_poly.h"
#include "portability.h"

// FIXME: make this depend on nb_polys...
const char * sidenames[NB_POLYS_MAX] = { "side-0", "side-1", "side-2", };

/* Be conservative and allocate two polynomials by default. */
void cado_poly_init(cado_poly poly)
{
    /* ALL fields are zero upon init, EXCEPT the degree field (which is -1) */
    memset(poly, 0, sizeof(poly[0]));

    /* By default allocate 2 polynomials */
    poly->nb_polys = 2;
    for(int side = 0 ; side < poly->nb_polys ; side++)
      mpz_poly_init (poly->pols[side], MAXDEGREE);

    mpz_init_set_ui(poly->n, 0);
}

void cado_poly_clear(cado_poly poly)
{
    for(int side = 0 ; side < poly->nb_polys ; side++)
      mpz_poly_clear (poly->pols[side]);

    mpz_clear(poly->n);
    memset(poly, 0, sizeof(poly[0]));
}

/* p <- q */
void
cado_poly_set (cado_poly p, cado_poly q)
{
    mpz_set (p->n, q->n);
    p->skew = q->skew;
    p->nb_polys = q->nb_polys;
    for(int side = 0 ; side < q->nb_polys ; side++)
      mpz_poly_set (p->pols[side], q->pols[side]);
}

// This function is no longer exported
#define BUF_MAX 10000

// returns 0 on failure, 1 on success.
static
int cado_poly_set_plist(cado_poly poly, param_list pl)
{
  int ret = 1;

  /* Parse skew value. Set to 0.0 to ensure that we get an invalid skewness
     in case it is not given */
  poly->skew = 0.0;
  param_list_parse_double(pl, "skew", &(poly->skew));

  /* Parse polynomials. Either in line format (poly%d=%Zd,%Zd,...) either given
   * by coefficients. */
  for(int i = 0; i < NB_POLYS_MAX; i++)
  {
    char tag[6], buf[BUF_MAX];
    snprintf(tag, sizeof(tag), "poly%d", i);
    if(param_list_parse_string(pl, tag, buf, BUF_MAX))
    {
      /* feeding poly->pols[i] from buf="c0,c1,..." */
      char *tmp = buf;
      if(i >= 2)
        mpz_poly_init (poly->pols[i], MAXDEGREE);
      for(unsigned int j = 0; ; j++)
      {
        int ret = gmp_sscanf(tmp, "%Zd", poly->pols[i]->coeff[j]);
        if (ret != 1)
        {
          fprintf (stderr, "Error while parsing coefficients of degree %u "
                           "of polynomial %d in polynomial file\n", j, i);
          abort();
        }
        for( ; *tmp != '\0' && *tmp != ','; tmp++);
        if(*tmp == '\0')
          break;
        tmp++;
      }
    }
    else
    {
      poly->nb_polys = i;
      break;
    }
  }
  if (poly->nb_polys == 0) /* No polynomials so far. Try other format. */
  /* Convention: c or X means side 1 (often algebraic)
   *             Y means side 0 (often rational) */
  {
    poly->nb_polys = 2;
    /* reading polynomials coefficient by coefficient */
    for (unsigned int i = 0; i <= MAXDEGREE; i++)
    {
      char tag[4];
      snprintf(tag, sizeof(tag), "c%d", i);
      int b = param_list_parse_mpz(pl, tag, poly->pols[1]->coeff[i]);
      if (!b) 
      {
        snprintf(tag, sizeof(tag), "X%d", i);
        param_list_parse_mpz(pl, tag, poly->pols[1]->coeff[i]);
      }
      snprintf(tag, sizeof(tag), "Y%d", i);
      param_list_parse_mpz(pl, tag, poly->pols[0]->coeff[i]);
    }
  }
  /* setting degrees */
  for(int side = 0 ; side < poly->nb_polys ; side++)
    mpz_poly_cleandeg(poly->pols[side], MAXDEGREE);

  /* Parse value of N. Two keys possible: n or None. Return 0 if not found. */
  if (!param_list_parse_mpz(pl, "n", poly->n) &&
      !param_list_parse_mpz(pl, NULL, poly->n))
  {
    fprintf (stderr, "Error, no value for N in cado_poly_set_plist\n");
    ret = 0;
  }

  for (int side = 0; side < poly->nb_polys; side++)
    if (poly->pols[side]->deg < 0)
    {
      fprintf (stderr, "Error, polynomial on side %u has degree < 0 in "
                       "cado_poly_set_plist\n", side);
      ret = 0;
    }

  /* check that N divides the resultant */
  cado_poly_getm (NULL, poly, poly->n);

  return ret;
}

// returns 0 on failure, 1 on success.
int cado_poly_read_stream (cado_poly_ptr poly, FILE * f)
{
  param_list pl;
  param_list_init (pl);
  param_list_read_stream (pl, f, 0);
  int r = cado_poly_set_plist (poly, pl);
  param_list_clear (pl);
  return r;
}

// returns 0 on failure, 1 on success.
int cado_poly_read_next_poly_from_stream (cado_poly_ptr poly, FILE * f)
{
  int r;
  param_list pl;
  param_list_init (pl);
  r = param_list_read_stream (pl, f, 1);
  if (r && pl->size > 0)
    r = cado_poly_set_plist (poly, pl);
  else
    r = 0;
  param_list_clear (pl);
  return r;
}

// returns 0 on failure, 1 on success.
int cado_poly_read(cado_poly poly, const char *filename)
{
    FILE *file;
    int r;
    file = fopen(filename, "r");
    if (file == NULL)
    {
      fprintf(stderr, "Error, in cado_poly_read: could not open %s\n", filename);
      return 0;
    }
    r = cado_poly_read_stream(poly, file);
    fclose(file);
    return r;
}


/* If NULL is passed for m, then, just check that N divides the resultant.
 * Assume there at least two polynomials.
 * TODO: adapt for more than 2 polynomials:
 *    compute for each pair (0,i) the corresponding common root m_i
 *    check all m_i are equal
 */
int cado_poly_getm(mpz_ptr m, cado_poly_ptr cpoly, mpz_ptr N)
{
    // have to work with copies, because pseudo_gcd destroys its input
    mpz_poly_t f[2];

    ASSERT_ALWAYS(cpoly->nb_polys >= 2);

    for (int i = 0; i < 2; ++i) {
        mpz_poly_init(f[i], cpoly->pols[i]->alloc);
        mpz_poly_set(f[i], cpoly->pols[i]);
    }
    int ret;
    ret = mpz_poly_pseudogcd_mpz (f[0], f[1], N, m);

    /* if ret=0, a factor of N was found during the pseudo-gcd,
       we assume N divides the resultant */

    if (m == NULL) {
        if (f[0]->deg < 1) {
            fprintf (stderr, "Error, N does not divide resultant of given polynomials\n");
            ASSERT_ALWAYS(0);
        }
        goto clear_and_return;
    }

    if (ret) /* pseudo-gcd was successful */
      {
        if (f[0]->deg != 1)
          {
            fprintf (stderr, "Error, N does not divide resultant of given polynomials, or there are multiplicities.\n");
            fprintf (stderr, "Can not compute m.\n");
            exit (EXIT_FAILURE);
          }

        mpz_t inv;
        mpz_init(inv);
        int ret2 = mpz_invert(inv, f[0]->coeff[1], N);
        // This inversion should always work.
        // If not, it means that N has a small factor (not sure we want
        // to be robust against that...)
        // Or maybe the polynomial selection was really bogus!
        ASSERT_ALWAYS(ret2);
        mpz_mul(inv, inv, f[0]->coeff[0]);
        mpz_neg(inv, inv);
        mpz_mod(m, inv, N);
        mpz_clear(inv);
      }

 clear_and_return:
    for (int i = 0; i < 2; ++i)
        mpz_poly_clear(f[i]);
    return ret;
}

/* Return the rational side or -1 if only algebraic sides */
/* Assume that there is at most 1 rational side (renumber.c does the same
 * assumption).
 * XXX: Is there a case where more than 1 rational side is needed ??
 */
int
cado_poly_get_ratside (cado_poly_ptr pol)
{
  for(int side = 0; side < pol->nb_polys; side++)
    if(pol->pols[side]->deg == 1)
      return side;
  return -1;
}

void
cado_poly_fprintf (FILE *fp, cado_poly_srcptr poly, const char *prefix)
{
  if (prefix)
    fputs (prefix, fp);
  gmp_fprintf (fp, "n: %Zd\n", poly->n);

  if (poly->nb_polys == 2)
  {
    mpz_poly_fprintf_cado_format (fp, poly->pols[0], 'Y', prefix);
    mpz_poly_fprintf_cado_format (fp, poly->pols[1], 'c', prefix);
  }
  else
  {
    for (int side = 0; side < poly->nb_polys; side++)
    {
      if (prefix)
        fputs (prefix, fp);
      fprintf (fp, "poly%u=", side);
      mpz_poly_fprintf_coeffs (fp, poly->pols[side], ',');
    }
  }

  if (prefix)
    fputs (prefix, fp);
  fprintf (fp, "skew: %1.3f\n", poly->skew);
}

void
cado_poly_fprintf_info (FILE *fp, double lognorm, double alpha,
                        double alpha_proj, unsigned int nrroots,
                        const char *prefix)
{
  if (prefix)
    fputs (prefix, fp);
  /* Always print "# " after the prefix and before the info line. */
  fprintf (fp, "# lognorm %1.2f, alpha %1.2f (proj %1.2f), E %1.2f, "
               "%u real root%s\n", lognorm, alpha, alpha_proj, lognorm + alpha,
               nrroots, (nrroots <= 1) ? "" : "s");
}

void
cado_poly_fprintf_MurphyE (FILE *fp, double MurphyE, double bound_f,
                           double bound_g, double area, const char *prefix)
{
  if (prefix)
    fputs (prefix, fp);
  /* Always print "# " after the prefix and before the MurphyE line. */
  fprintf (fp, "# MurphyE(Bf=%.2e,Bg=%.2e,area=%.2e)=%.2e\n", bound_f, bound_g,
               area, MurphyE);
}
