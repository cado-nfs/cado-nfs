/* Main function to invoke cachesize_guess() or cachesize_cpuid() */
#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include "version_info.h"

int cachesize_cpuid(int verbose);
size_t cachesize_guess(int verbose);

int main(int argc, char const * argv[])
{
  int i;
  fprintf (stderr, "# %s.r%s", *argv, cado_revision_string);
  for (i = 1; i < argc; i++)
    fprintf (stderr, " %s", *(argv+i));
  fprintf (stderr, "\n");

  if (argc != 1) {
    printf ("Usage: %s\n", *argv);
    return EXIT_FAILURE;
  }

  {
      printf ("-- invoking cachesize_cpuid() --\n");
      int ret = cachesize_cpuid (1);
      printf ("-- cachesize_cpuid() returns %d --\n", ret);
  }

  {
      printf ("-- invoking cachesize_guess() --\n");
      size_t ret = cachesize_guess (1);
      printf ("-- cachesize_guess() returns %zu --\n", ret);
  }

  return EXIT_SUCCESS;
}
