/* complete_rels: same as complete_rels, but completes small factors if omitted */
#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#include <time.h>
#include <gmp.h>
#include "utils.h"

void
norm (mpz_t *f, int deg, mpz_t r, long a, unsigned long b)
{
  int i;
  mpz_t tmp;

  mpz_init_set_ui (tmp, 1);
  mpz_set (r, f[deg]);
  for (i = deg - 1; i >= 0; i--)
    {
      mpz_mul_si (r, r, a);
      mpz_mul_ui (tmp, tmp, b);
      mpz_addmul (r, tmp, f[i]);
    }
  mpz_abs (r, r);
  mpz_clear (tmp);
}

void
add_alg (relation_t *rel, unsigned long p)
{
  relation_provision_for_primes(rel, 0, rel->nb_ap + 1);
  rel->ap[rel->nb_ap++] = (alg_prime_t) { .p=p, .e=1 };
}

void
add_rat (relation_t *rel, unsigned long p)
{
  relation_provision_for_primes(rel, rel->nb_rp + 1, 0);
  rel->rp[rel->nb_rp++] = (rat_prime_t) { .p=p, .e=1 };
}

int
complete_relation (relation_t *rel, cado_poly_ptr cpoly)
{
  mpz_t no;
  int i, j;
  unsigned long p;

  mpz_init (no);

  // algebraic side
  norm (cpoly->f, cpoly->degree, no, rel->a, rel->b);
  for (i = 0; i < rel->nb_ap; ++i)
    {
      for (j = 0; j < (rel->ap[i]).e; ++j)
	if (mpz_divisible_ui_p (no, (rel->ap[i]).p) == 0)
	  {
	    fprintf (stderr,
                     "Wrong algebraic side for (%" PRId64 ", %" PRIu64 ")\n",
		     rel->a, rel->b);
	    fprintf (stderr, "Given factor %lu does not divide norm\n",
		     (rel->ap[i]).p);
	    mpz_clear (no);
	    return 0;
	  }
	else
	  mpz_divexact_ui (no, no, (rel->ap[i]).p);
    }

  for (p = 2; mpz_cmp_ui (no, 1) != 0; p += 1 + (p != 2))
    {
      while (mpz_divisible_ui_p (no, p))
	{
	  add_alg (rel, p);
	  mpz_divexact_ui (no, no, p);
	}
    }

  // rational side
  norm (cpoly->g, 1, no, rel->a, rel->b);
  for (i = 0; i < rel->nb_rp; ++i)
    {
      for (j = 0; j < (rel->rp[i]).e; ++j)
	if (mpz_divisible_ui_p (no, (rel->rp[i]).p) == 0)
	  {
	    fprintf (stderr,
                     "Wrong rational side for (%" PRId64 ", %" PRIu64 ")\n",
		     rel->a, rel->b);
	    fprintf (stderr, "Given factor %lu does not divide norm\n",
		     (rel->rp[i]).p);
	    mpz_clear (no);
	    return 0;
	  }
	else
	  mpz_divexact_ui (no, no, (rel->rp[i]).p);
    }
  for (p = 2; mpz_cmp_ui (no, 1) != 0; p += 1 + (p != 2))
    {
      while (mpz_divisible_ui_p (no, p))
	{
	  add_rat (rel, p);
	  mpz_divexact_ui (no, no, p);
	}
    }

  // were printing our output, thus we don't care about tidyness of the
  // internal structure.
  // relation_compress_rat_primes(rel);
  // relation_compress_alg_primes(rel);

  return 1;
}

int complete_relation_files(char ** files, cado_poly_ptr cpoly)
{
    relation_stream rs;
    relation_stream_init(rs);
    unsigned long ok = 0;
    unsigned long bad = 0;
    int had_error = 0;
    for( ; *files ; files++) {
        relation_stream_openfile(rs, *files);
        char line[RELATION_MAX_BYTES];
        unsigned long l0 = rs->lnum;
        unsigned long ok0 = ok;
        unsigned long bad0 = bad;
        for( ; relation_stream_get(rs, line) >= 0 ; ) {
            unsigned long l = rs->lnum - l0;
            if (complete_relation(&rs->rel, cpoly)) {
                fprint_relation (stdout, rs->rel);
                ok++;
                continue;
            }
            fprintf(stderr, "Failed at line %lu in %s: %s",
                    l, *files, line);
            bad++;
            had_error = 1;
        }
        if (bad == bad0) {
            fprintf(stderr, "%s : %lu ok\n", *files, ok-ok0);
        } else {
            fprintf(stderr, "%s : %lu ok ; FOUND %lu ERRORS\n",
                    *files, ok-ok0, bad-bad0);
        }
        relation_stream_closefile(rs);
    }
    relation_stream_clear(rs);

    return had_error ? -1 : 0;
}

void usage_and_die(char *str) {
    fprintf(stderr, "usage: %s -poly <polyfile> <relfile1> <relfile2> ...\n",
            str);
    exit(3);
}

int main(int argc, char * argv[])
{
    cado_poly cpoly;
    int had_error = 0;

    if (argc < 4 || strcmp(argv[1], "-poly") != 0) {
        usage_and_die(argv[0]);
    }
    cado_poly_init(cpoly);
    if (!cado_poly_read(cpoly, argv[2])) 
        return 2;

    argv++,argc--;
    argv++,argc--;
    argv++,argc--;

    had_error = complete_relation_files(argv, cpoly) < 0;

    cado_poly_clear(cpoly);

    return had_error;
}
