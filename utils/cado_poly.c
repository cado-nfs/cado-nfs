#include "cado.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "params.h"
#include "cado_poly.h"

void cado_poly_init(cado_poly poly)
{
    int i;

    /* ALL fields are zero upon init, EXCEPT the degree field (which is -1) */
    memset(poly, 0, sizeof(poly));
    poly->f = (mpz_t *) malloc((MAXDEGREE + 1) * sizeof(mpz_t));
    poly->g = (mpz_t *) malloc((MAXDEGREE + 1) * sizeof(mpz_t));
    /* mpzs as well are zero */
    for (i = 0; i < (MAXDEGREE + 1); i++) {
	mpz_init_set_ui(poly->f[i], 0);
	mpz_init_set_ui(poly->g[i], 0);
    }
    mpz_init_set_ui(poly->n, 0);
    poly->degree = -1;
    poly->degreeg = -1;
    poly->type[0] = '\0';
    mpz_init_set_ui(poly->m, 0);
}

void cado_poly_clear(cado_poly poly)
{
    int i;

    for (i = 0; i < (MAXDEGREE + 1); i++) {
	mpz_clear(poly->f[i]);
	mpz_clear(poly->g[i]);
    }
    mpz_clear(poly->n);
    free(poly->f);
    free(poly->g);
    mpz_clear(poly->m);
    memset(poly, 0, sizeof(poly));
}

int cado_poly_set_plist(cado_poly poly, param_list pl)
{
    int have_n = 0;
    int have_f[(MAXDEGREE + 1)] = { 0, };
    int have_g[(MAXDEGREE + 1)] = { 0, };
    int have_m = 0;
    int degf, degg;
    int i;
    mpz_t tmp;

    have_n = param_list_parse_mpz(pl, "n", poly->n) || param_list_parse_mpz(pl, NULL, poly->n);
    poly->skew = 0.0; /* to ensure that we get an invalid skewness in case
                         it is not given */
    param_list_parse_double(pl, "skew", &(poly->skew));
    for (i = 0; i < (MAXDEGREE + 1); i++) {
        char tag[4];
        snprintf(tag, sizeof(tag), "c%d", i);
        have_f[i] = param_list_parse_mpz(pl, tag, poly->f[i]);
        if (!have_f[i]) {
            snprintf(tag, sizeof(tag), "X%d", i);
            have_f[i] = param_list_parse_mpz(pl, tag, poly->f[i]);
        }
        snprintf(tag, sizeof(tag), "Y%d", i);
        have_g[i] = param_list_parse_mpz(pl, tag, poly->g[i]);
    }
    param_list_parse_string(pl, "type", poly->type, sizeof(poly->type));
    param_list_parse_ulong(pl, "rlim", &(poly->rlim));
    param_list_parse_ulong(pl, "alim", &(poly->alim));
    param_list_parse_int(pl, "lpbr", &(poly->lpbr));
    param_list_parse_int(pl, "lpba", &(poly->lpba));
    param_list_parse_int(pl, "mfbr", &(poly->mfbr));
    param_list_parse_int(pl, "mfba", &(poly->mfba));
    param_list_parse_double(pl, "rlambda", &(poly->rlambda));
    param_list_parse_double(pl, "alambda", &(poly->alambda));
    param_list_parse_int(pl, "qintsize", &(poly->qintsize));
    mpz_set_ui (poly->m, 0);
    param_list_parse_mpz(pl, "m", poly->m);
    have_m = mpz_cmp_ui (poly->m, 0);

    for (degf = MAXDEGREE; degf >= 0 && !have_f[degf]; degf--) {
	if (have_f[degf] && mpz_cmp_ui(poly->f[degf], 0) != 0)
	    break;
    }
    for (degg = MAXDEGREE; degg >= 0 && !have_g[degg]; degg--) {
	if (have_g[degg] && mpz_cmp_ui(poly->g[degg], 0) != 0)
	    break;
    }

    ASSERT_ALWAYS(have_n);

    poly->degree = degf;
    poly->degreeg = degg;
    // compute m, the common root of f and g mod n
    mpz_init (tmp);
    if (degg == 1)
      {
	mpz_invert (tmp, poly->g[1], poly->n);
	mpz_mul (tmp, tmp, poly->g[0]);
	mpz_mod (tmp, tmp, poly->n);
	mpz_sub (tmp, poly->n, tmp);
	mpz_mod (tmp, tmp, poly->n);

        if (have_m && (mpz_cmp (poly->m, tmp) != 0))
          {
            fprintf (stderr, "m is not a root of g mod N\n");
            exit (EXIT_FAILURE);
          }
	mpz_set (poly->m, tmp);
      }
    else
      {
        if (have_m == 0)
          {
            fprintf (stderr, "must provide m for non-linear polynomials\n");
            exit (EXIT_FAILURE);
          }

        /* check m is a root of f mod n */
        mpz_set (tmp, poly->f[degf]);
        for (i = degf - 1; i >= 0; i--)
          {
            mpz_mul (tmp, tmp, poly->m);
            mpz_add (tmp, tmp, poly->f[i]);
          }
        mpz_mod (tmp, tmp, poly->n);
        if (mpz_cmp_ui (tmp, 0) != 0)
          {
            fprintf (stderr, "m is not a root of f modulo n\n");
            exit (EXIT_FAILURE);
          }

        /* check m is a root of g mod n */
        mpz_set (tmp, poly->g[degg]);
        for (i = degg - 1; i >= 0; i--)
          {
            mpz_mul (tmp, tmp, poly->m);
            mpz_add (tmp, tmp, poly->g[i]);
          }
        mpz_mod (tmp, tmp, poly->n);
        if (mpz_cmp_ui (tmp, 0) != 0)
          {
            fprintf (stderr, "m is not a root of g modulo n\n");
            exit (EXIT_FAILURE);
          }
      }
    mpz_clear (tmp);

    return 1;
}

// returns 0 on failure, 1 on success.
int cado_poly_read_stream(cado_poly poly, FILE * f)
{
    param_list pl;
    param_list_init(pl);
    param_list_read_stream(pl, f);
    int r = cado_poly_set_plist(poly, pl);
    param_list_clear(pl);
    return r;
}

// returns 0 on failure, 1 on success.
int cado_poly_read(cado_poly poly, const char *filename)
{
    FILE *file;
    int r;
    file = fopen(filename, "r");
    if (file == NULL) {
	fprintf(stderr, "read_polynomial: could not open %s\n", filename);
	return 0;
    }
    r = cado_poly_read_stream(poly, file);
    fclose(file);
    return r;
}


void fprint_polynomial(FILE * fp, mpz_t * f, const int d)
{
    int i, s, first = 1;
    mpz_t c;

    mpz_init(c);
    for (i = d; i >= 0; i--) {
	s = mpz_cmp_ui(f[i], 0);
	if (s != 0) {
	    if (s > 0) {
		if (first == 0)
		    gmp_fprintf(fp, " + ");
		gmp_fprintf(fp, "%Zd", f[i]);
	    } else if (s < 0) {
		mpz_abs(c, f[i]);
		gmp_fprintf(fp, " - %Zd", c);
	    }
	    first = 0;
	    if (i >= 2)
		gmp_fprintf(fp, "*x^%d", i);
	    else if (i == 1)
		gmp_fprintf(fp, "*x");
	}
    }
    mpz_clear(c);
    fprintf(fp, "\n");
}

/* check that n divides b^d*f(m/d), where g = b*x-m */
void check_polynomials(cado_poly cpoly)
{
    mpz_t r, q;
    int k;
    int d = cpoly->degree;

    mpz_init_set(r, cpoly->f[d]);
    mpz_init_set_ui(q, 1);

    for (k = d - 1; k >= 0; k--) {
	mpz_mul(q, q, cpoly->g[1]);
	/* invariant: q = b^(d-k) */
	mpz_mul(r, r, cpoly->g[0]);
	mpz_neg(r, r);
	mpz_addmul(r, q, cpoly->f[k]);
    }

    if (mpz_divisible_p (r, cpoly->n) == 0)
      {
	fprintf (stderr, "Error, n does not divide Res(f,g)\n");
	exit (EXIT_FAILURE);
      }

    mpz_clear (r);
    mpz_clear (q);
}
