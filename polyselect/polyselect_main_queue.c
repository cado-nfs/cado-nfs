#include "cado.h" // IWYU pragma: keep

#include <float.h>
#include <limits.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <pthread.h>

#include <gmp.h>

#include "auxiliary.h"
#include "macros.h"
#include "mpz_poly.h"
#include "polyselect_main_queue.h"
#include "polyselect_norms.h"
#include "portability.h"


/* Read-Only */

#if 0
const double exp_rot[] = { 0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 0 };

unsigned int sopt_effort = SOPT_DEFAULT_EFFORT;	/* size optimization effort */
double best_E = 0.0;		/* Murphy's E (the larger the better) */

/* read-write global variables */
double potential_collisions = 0.0;
#endif

/* print poly info */
static size_t
snprintf_poly_info(char *buf,
		  size_t size,
		  mpz_poly_srcptr f,
		  mpz_poly_srcptr g,
		  mpz_srcptr n,
		  const int raw)
{
  size_t np = 0;

  np += snprintf(buf + np, size - np, "# %s polynomial:\n", raw ? "Raw" : "Size-optimized");
  ASSERT_ALWAYS(np <= size);

  char * str;

  np += gmp_snprintf(buf + np, size - np, "%sn: %Zd\n", raw ? "# " : "", n);
  ASSERT_ALWAYS(np <= size);

  mpz_poly_asprintf_cado_format(&str, g, 'Y', raw ? "# " : "");
  np += strlcat(buf + np, str, size - np);
  ASSERT_ALWAYS(np <= size);
  free(str);

  mpz_poly_asprintf_cado_format(&str, f, 'c', raw ? "# " : "");
  np += strlcat(buf + np, str, size - np);
  ASSERT_ALWAYS(np <= size);
  free(str);

  double skew = L2_skewness(f);
  unsigned int nroots = mpz_poly_number_of_real_roots(f);
  double logmu = L2_lognorm(f, skew);
  double exp_E = logmu + expected_rotation_gain(f, g);

  np += snprintf(buf + np, size - np, "# %sexp_E", raw ? "raw " : "");
  ASSERT_ALWAYS(np <= size);

  np += snprintf(buf + np, size - np,
	       " %1.2f, lognorm %1.2f, skew %1.2f, %u rroots\n", exp_E,
	       logmu, skew, nroots);
  ASSERT_ALWAYS(np <= size);

  return np;
}

void polyselect_fprintf_poly_pair(FILE * fp, mpz_srcptr N, mpz_poly_srcptr f, mpz_poly_srcptr g, int raw)
{
    if (!f && !g) return;
    size_t sz = mpz_sizeinbase(N, 10);
    int length = sz * 12;
    char *str = malloc(length);
    ASSERT_ALWAYS(str);
    snprintf_poly_info(str, length, f, g, N, raw);
    fprintf(fp, "%s", str);
    free(str);
}
