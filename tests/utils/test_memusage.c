#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include "memusage.h"

int
main ()
{
  long l = Memusage2 ();
  printf ("Memusage2: %ld\n", l);
  exit (EXIT_SUCCESS);
}
