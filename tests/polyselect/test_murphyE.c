#include "cado.h"
#include <math.h>
#include "utils.h"
#include "polyselect/murphyE.h"
#include "tests_common.h"

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

void
test_ncx2_pdf (void)
{
  double y, epsilon = 0.001;

  y = ncx2_pdf (1.0, 2.0, 1.0, epsilon);
  check_num (y, 0.232879803796820, epsilon);

  y = ncx2_pdf (1.0, 2.0, 2.0, epsilon);
  check_num (y, 0.17472016746112826, epsilon);

  y = ncx2_pdf (1.0, 2.0, 3.0, epsilon);
  check_num (y, 0.12876542477554598, epsilon);

  y = ncx2_pdf (1.0, 4.0, 1.0, epsilon);
  check_num (y, 0.10395520767485419, epsilon);

  y = ncx2_pdf (1.0, 4.0, 2.0, epsilon);
  check_num (y, 0.070939964617860438, epsilon);

  y = ncx2_pdf (1.0, 4.0, 3.0, epsilon);
  check_num (y, 0.048210398188727084, epsilon);

  y = ncx2_pdf (2.0, 4.0, 3.0, epsilon);
  check_num (y, 0.080557660618074192, epsilon);
}

int
main (int argc, const char *argv[])
{
  tests_common_cmdline (&argc, &argv, PARSE_SEED);
  test_ncx2_pdf ();
  tests_common_clear ();
  exit (EXIT_SUCCESS);
}

