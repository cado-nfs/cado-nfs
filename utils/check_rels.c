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

void norm(mpz_t *f, int deg, mpz_t r, long a, unsigned long b)
{
    int i;

    mpz_t tmp;
    mpz_init_set_ui(tmp, 1);
    mpz_set(r, f[deg]);
    for (i = deg - 1; i >= 0; i--) {
        mpz_mul_si(r, r, a);
        mpz_mul_ui(tmp, tmp, b);
        mpz_addmul(r, tmp, f[i]);
    }
    mpz_abs(r, r);
    mpz_clear(tmp);
}

int
check_relation (relation_t *rel, cado_poly_ptr cpoly)
{
  mpz_t no, acc;
  int i;

  mpz_init (no);
  mpz_init (acc);

  // algebraic side
  norm (cpoly->f, cpoly->degree, no, rel->a, rel->b);
  mpz_set_ui (acc, 1);
  for (i = 0; i < rel->nb_ap; ++i)
    {
      int j;
      for (j = 0; j < (rel->ap[i]).e; ++j) 
	mpz_mul_ui (acc, acc, (rel->ap[i]).p);
    }
  if (mpz_cmp (acc, no) != 0)
    {
      if (mpz_divisible_p (no, acc))
	{
	  mpz_divexact (acc, no, acc);
	  gmp_fprintf (stderr, "Missing factor %Zd on algebraic side for (%ld, %lu)\n", acc, rel->a, rel->b);
	}
      else
	{
	  mpz_t g;
	  mpz_init (g);
	  mpz_gcd (g, acc, no);
	  mpz_divexact (acc, acc, g);
	  fprintf (stderr,
		   "Wrong algebraic side for (%" PRId64 ", %" PRIu64 ")\n",
                   rel->a, rel->b);
	  gmp_fprintf (stderr, "Given factor %Zd does not divide norm\n", acc);
	  mpz_clear (g);
	}
      mpz_clear (no);
      mpz_clear (acc);
      return 0;
    }

    // rational side
    norm(cpoly->g, 1, no, rel->a, rel->b);
    mpz_set_ui(acc, 1);
    for(i = 0; i < rel->nb_rp; ++i) {
        int j;
        for (j = 0; j < (rel->rp[i]).e; ++j) 
            mpz_mul_ui(acc, acc, (rel->rp[i]).p);
    }
    if (mpz_cmp(acc, no) != 0) {
        fprintf (stderr,
                 "Wrong rational side for (%" PRId64 ", %" PRIu64 ")\n",
                 rel->a, rel->b);
        mpz_clear(no);
        mpz_clear(acc);
        return 0;
    }
    mpz_clear(no);
    mpz_clear(acc);
    return 1;
}

void usage_and_die(char *str) {
    fprintf(stderr, "usage: %s -poly <polyfile> <relfile1> <relfile2> ...\n",
            str);
    exit(3);
}

int check_relation_files(char ** files, cado_poly_ptr cpoly)
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
            if (check_relation(&rs->rel, cpoly)) {
                printf("%s", line);
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

    had_error = check_relation_files(argv, cpoly) < 0;

    cado_poly_clear(cpoly);

    return had_error;
}
