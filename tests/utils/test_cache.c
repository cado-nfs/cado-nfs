#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include "cachesize_cpuid.h"

int
main ()
{
  cachesize_cpuid (1);
  exit (EXIT_SUCCESS);
}
