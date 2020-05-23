#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include <stdio.h>
#include "macros.h"
#include "sieving_bound.h"

void sieving_bound_init(sieving_bound_ptr H, unsigned int t)
{
  ASSERT(t > 2);

  H->t = t;
  H->h = (unsigned int * ) malloc(sizeof(unsigned int) * t);
}

void sieving_bound_set_hi(sieving_bound_ptr H, unsigned int i,
    unsigned int value)

{
  ASSERT(i < H->t);
#ifndef NDEBUG
  if (i == H->t - 1) {
    ASSERT(value > 1);
  }
#endif // NDEBUG

  H->h[i] = value;
}

uint64_t sieving_bound_number_element(sieving_bound_srcptr H)
{
  uint64_t nb = 1;
  for (unsigned int i = 0; i < H->t - 1; i++) {
    nb = nb * (2 * (uint64_t)H->h[i]);
  }
  nb = nb * ((uint64_t)H->h[H->t - 1]);

  ASSERT(nb != 1);

  return nb;
}

void sieving_bound_clear(sieving_bound_ptr H)
{
  free(H->h);
  H->t = 0;
}

void sieving_bound_fprintf(FILE * filew, sieving_bound_srcptr H)
{
  ASSERT(H->t > 2);

  fprintf(filew, "[");
  for (unsigned int i = 0 ; i < H->t - 1; i++) {
    fprintf(filew, "%u, ", H->h[i]);
  }
  fprintf(filew, "%u]\n", H->h[H->t - 1]);
}

void sieving_bound_fprintf_detailed(FILE * filew, sieving_bound_srcptr H)
{
  ASSERT(H->t > 2);

  for (unsigned int i = 0 ; i < H->t - 1; i++) {
    fprintf(filew, "-H%u: -%u -- H%u: %u\n", i, H->h[i], i, H->h[i] - 1);
  }
  fprintf(filew, "H%u: 0 -- H%u: %u\n", H->t - 1, H->t - 1, H->h[H->t - 1] - 1);
}

void sieving_bound_fprintf_detailed_comment(FILE * filew,
    sieving_bound_srcptr H)
{
  ASSERT(H->t > 2);

  fprintf(filew, "# [-H%u: -%u -- H%u: %u\n", 0, H->h[0], 0, H->h[0] - 1);

  for (unsigned int i = 1 ; i < H->t - 1; i++) {
    fprintf(filew, "# -H%u: -%u -- H%u: %u\n", i, H->h[i], i, H->h[i] - 1);
  }
  fprintf(filew, "# H%u: 0 -- H%u: %u]\n", H->t - 1, H->t - 1,
      H->h[H->t - 1] - 1);
}
