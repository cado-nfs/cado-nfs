#include "cado.h" // IWYU pragma: keep
#include <stdio.h>
#include <stdlib.h>
#include "memusage.h"

int
main ()
{
    size_t l = Memusage2 ();
    printf ("Memusage2: %zu\n", l);
    return EXIT_SUCCESS;
}
