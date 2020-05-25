/* gen_decomp mfb fbb computes approximations of all decompositions of
   mfb-bit integers into primes larger than fbb. Example:
$ gen_decomp 60 524288
60 20 40 4.700208e+14 # 1
60 21 39 4.596914e+14 # 2
60 22 38 4.502935e+14 # 3
60 23 37 4.420430e+14 # 4
60 24 36 4.351291e+14 # 5
60 25 35 4.297433e+14 # 6
60 26 34 4.254332e+14 # 7
60 27 33 4.214013e+14 # 8
60 28 32 4.197615e+14 # 9
60 29 31 4.183856e+14 # 10
60 30 30 2.088895e+14 # 11
60 20 41 5.930007e+14 # 12
60 21 40 5.788613e+14 # 13
60 22 39 5.664933e+14 # 14
60 23 38 5.560455e+14 # 15
60 24 37 5.460721e+14 # 16
60 25 36 5.393045e+14 # 17
60 26 35 5.322997e+14 # 18
60 27 34 5.281604e+14 # 19
60 28 33 5.249884e+14 # 20
60 29 32 5.234872e+14 # 21
60 30 31 5.197604e+14 # 22
60 20 20 20 2.477904e+12 # 23
60 20 20 21 3.529032e+13 # 24
60 20 20 22 1.034259e+13 # 25
60 20 21 21 1.030657e+13 # 26
*/

#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include "macros.h"
#include "gen_decomp.h"
#include "tab_decomp.h"

/*
  To test the function generate_all_decomp ()!!
*/

int
main (int argc, char *argv[])
{
  unsigned long lim;
  int mfb;
  ASSERT_ALWAYS (argc == 3);
  mfb = atoi (argv[1]);
  lim = atol (argv[2]);
  tabular_decomp_t* res = generate_all_decomp (mfb, lim);

  tabular_decomp_print (res);

  tabular_decomp_free (res);
  
  return EXIT_SUCCESS;
}
