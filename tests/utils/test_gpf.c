#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include <stdio.h>
#include "gpf.h"

/* Generated with Pari/GP
   for(i=2, 150, print1(vecmax(factorint(4)~[1,]), ", ")) */
const unsigned int
known_gpf[151] = {0, 1,
  2, 3, 2, 5, 3, 7, 2, 3, 5, 11, 3, 13, 7, 5, 2, 17, 3, 19, 5, 7, 11, 23, 3,
  5, 13, 3, 7, 29, 5, 31, 2, 11, 17, 7, 3, 37, 19, 13, 5, 41, 7, 43, 11, 5,
  23, 47, 3, 7, 5, 17, 13, 53, 3, 11, 7, 19, 29, 59, 5, 61, 31, 7, 2, 13, 11,
  67, 17, 23, 7, 71, 3, 73, 37, 5, 19, 11, 13, 79, 5, 3, 41, 83, 7, 17, 43,
  29, 11, 89, 5, 13, 23, 31, 47, 19, 3, 97, 7, 11, 5, 101, 17, 103, 13, 7,
  53, 107, 3, 109, 11, 37, 7, 113, 19, 23, 29, 13, 59, 17, 5, 11, 61, 41, 31,
  5, 7, 127, 2, 43, 13, 131, 11, 19, 67, 5, 17, 137, 23, 139, 7, 47, 71, 13,
  3, 29, 73, 7, 37, 149, 5
};


int test(unsigned int n)
{
  gpf_init(n);
  for (unsigned int i = 0; i <= n; i++) {
    if (gpf_get(i) != known_gpf[i])  {
      fprintf(stderr, "get_gpf(%u) = %u wrong, %u is correct\n",
              i, gpf_get(i), known_gpf[i]);
      return 0;
    }
  }
  return 1;
}

int main ()
{
  int rc = EXIT_SUCCESS;
  if (test(100) == 0) {
    rc = EXIT_FAILURE;
  }
  if (test(150) == 0)
    rc = EXIT_FAILURE;
  exit(rc);
}
