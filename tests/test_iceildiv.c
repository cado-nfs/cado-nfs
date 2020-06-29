#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>
#include "macros.h"

static int
do_test_siceildiv(int i, int j, double ref)
{
  int sicd = siceildiv(i, j);
  if (ceil(ref) != sicd) {
    printf("siceildiv(%d, %d) = %d != ceil(%f)\n", i, j, sicd, ref);
    return 1;
  }
  return 0;
}

static int
do_test_iceildiv(unsigned int k, unsigned int l, double ref)
{
  unsigned int icd = iceildiv(k, l);
  if (ceil(ref) != icd) {
    printf("iceildiv(%u, %u) = %u != ceil(%f)\n", k, l, icd, ref);
    return 1;
  }
  return 0;
}

int main()
{
  int i, j, had_error = 0;
  unsigned int k, l;
  for (i=-3; i <= 3; i++)
    for (j=-3; j <= 3; j++)
      if (j != 0) {
        had_error |= do_test_siceildiv(i, j, (double) i / (double) j);
      }
  had_error |= do_test_siceildiv(INT_MAX, INT_MAX, 1.0);
  had_error |= do_test_siceildiv(INT_MIN, INT_MIN, 1.0);

  for (k=0; k <= 3; k++)
    for (l=1; l <= 3; l++)
      had_error |= do_test_iceildiv(k, l, (double) k / (double) l);
  had_error |= do_test_iceildiv(UINT_MAX, UINT_MAX, 1.0);

  if (had_error) {
    exit(EXIT_FAILURE);
  } else {
    exit(EXIT_SUCCESS);
  }
}
