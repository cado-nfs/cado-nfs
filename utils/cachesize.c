/* Main function to invoke cachesize_guess() or cachesize_cpuid() */
#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include "version_info.h"

int cachesize_cpuid(int verbose);
int cachesize_guess(int verbose);

int main (int argc, char **argv)
{
  int i, ret;
  fprintf (stderr, "# %s.r%s", *argv, cado_revision_string);
  for (i = 1; i < argc; i++)
    fprintf (stderr, " %s", *(argv+i));
  fprintf (stderr, "\n");

  if (argc != 1) {
    printf ("Usage: %s\n", *argv);
    exit (1);
  }

  printf ("-- invoking cachesize_cpuid() --\n");  
  ret = cachesize_cpuid (1);
  printf ("-- cachesize_cpuid() returns %d --\n", ret);  

  printf ("-- invoking cachesize_guess() --\n");  
  ret = cachesize_guess (1);
  printf ("-- cachesize_guess() returns %d --\n", ret);  

}
