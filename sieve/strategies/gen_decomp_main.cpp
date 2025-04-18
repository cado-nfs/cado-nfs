#include "cado.h" // IWYU pragma: keep

#include <cstdlib>
#include <climits>

#include "macros.h"
#include "gen_decomp.hpp"

int main(int argc, char const * argv[])
{
  ASSERT_ALWAYS (argc == 3);
  char * p;
  long mfb = strtol(argv[1], &p, 0);
  ASSERT_ALWAYS(*p == '\0');
  ASSERT_ALWAYS(mfb <= INT_MAX);
  unsigned long lim = strtoul (argv[2], &p, 0);
  ASSERT_ALWAYS(*p == '\0');

  generate_all_decomp((int) mfb, lim);
  
  return EXIT_SUCCESS;
}
