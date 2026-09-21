#include "cado.h" // IWYU pragma: keep

#include <cstdlib>
#include <climits>

#include "macros.h"
#include "gen_decomp.hpp"
#include "cado_main.hpp"

static int main_(int argc, char const * argv[]);

int main(int argc, char const * argv[])
{
    return cado::main_wrapper(main_, argc, argv);
}

static int main_(int argc, char const * argv[])
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
