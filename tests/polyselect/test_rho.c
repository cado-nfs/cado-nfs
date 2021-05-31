#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tests_common.h"
#include "macros.h"
#include "rho.h"
#include "misc.h"

/* check relative error is less than emax */
static int
check_num (double x, double y, double emax)
{
  double e = fabs (x - y) / fabs (y);

  if (e > emax)
    {
      printf ("expected %.16e, got %.16e (rel. error %e)\n", y, x, e);
      return 0;
    }
  return 1;
}

/* check absolute error is less than emax */
static int
check_num_abs (double x, double y, double emax)
{
  double e = fabs (x - y);

  if (e > emax)
    {
      printf ("expected %.16e, got %.16e (abs. error %e)\n", y, x, e);
      return 0;
    }
  return 1;
}

void
test_rho (void)
{
  double x, y;

  y = dickman_rho (-0.5);
  ASSERT_ALWAYS (y == 0.0);

  y = dickman_rho (0.0);
  ASSERT_ALWAYS (y == 1.0);

  x = ldexp(u64_random(state), -64);
  ASSERT_ALWAYS (0 <= x && x < 1.0);
  y = dickman_rho (x);
  ASSERT_ALWAYS (y == 1.0);

  y = dickman_rho (1.1);
  check_num (y, 0.904689820195675, 1.2e-9);

  y = dickman_rho (2.2);
  check_num (y, 0.220357137908328, 8.9e-10);

  y = dickman_rho (3.3);
  check_num (y, 0.0254647238733285, 7.7e-10);

  y = dickman_rho (4.4);
  check_num (y, 0.00177994246481535, 8.0e-10);

  y = dickman_rho (5.5);
  check_num (y, 0.0000860186111205116, 1.1e-9);

  y = dickman_rho (6.6);
  check_num (y, 3.11012649979137e-6, 1.6e-9);

  y = dickman_rho (7.7);
  check_num (y, 8.85046647687321e-8, 2.6e-9);

  y = dickman_rho (8.8);
  check_num (y, 2.05443505293307e-9, 4.4e-9);

  y = dickman_rho (9.9);
  check_num (y, 3.99531836601083e-11, 7.5e-9);

  y = dickman_rho (10.1);
  check_num (y, 1.91826261797451e-11, 1.2e-8);

  y = dickman_rho (11.2);
  check_num (y, 3.10667427553834420e-13, 1.9e-8);

  y = dickman_rho (12.3);
  check_num (y, 4.38519652833446e-15, 2.7e-8);

  y = dickman_rho (13.4);
  check_num (y, 5.46466232309370e-17, 3.7e-8);

  y = dickman_rho (14.5);
  check_num (y, 6.07650960951011e-19, 4.9e-8);

  y = dickman_rho (15.6);
  check_num_abs (y, 6.08381226695129e-21, 3.5e-21);
}

int
main (int argc, const char *argv[])
{
  tests_common_cmdline (&argc, &argv, PARSE_SEED);
  test_rho ();
  tests_common_clear ();
  exit (EXIT_SUCCESS);
}

