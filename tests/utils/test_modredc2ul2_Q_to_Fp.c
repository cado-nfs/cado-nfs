#include "cado.h"
#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>
#include "tests_common.h"
#include "getprime.h"
#include "modredc_2ul2.h"

#define P_LEN 128

int
main (int argc, const char *argv[])
{
  unsigned long iter = 100000;
  tests_common_cmdline (&argc, &argv, PARSE_ITER | PARSE_VERBOSE);
  tests_common_get_iter (&iter);

  unsigned long p[P_LEN], r[P_LEN];
  unsigned long t;
  prime_info pi;

  prime_info_init (pi);
  t = getprime_mt (pi);
  while(t < 8192)
    t = getprime_mt (pi);
  for (size_t i = 0; i < P_LEN; i++) {
    p[i] = t;
    t = getprime_mt (pi);
  }
  prime_info_clear (pi);

  const unsigned long k = 0;
  const int neg = 0;
  modredc2ul2_batch_Q_to_Fp_context_t *context;

  unsigned long p1[2] = {451UL, 0UL}, p2[2] = {41UL, 0UL};

#if ULONG_BITS == 32
      p1[1] = 452984832UL; /* p1 = 27*2^56+451 */
      p2[1] = 184549376UL; /* p2 = 11*2^56+41 */
#elif ULONG_BITS == 64
      p1[1] = 1UL << 56; /* p1 = 2^120+415 */
      p2[1] = 1UL << 57; /* p2 = 2^121+41 */
#else
#error "uh?"
#endif

  modintredc2ul2_t num, den;
  modredc2ul2_intinit(num);
  modredc2ul2_intinit(den);
  modredc2ul2_intset_uls(num, p1, 2);
  modredc2ul2_intset_uls(den, p2, 2);
  context = modredc2ul2_batch_Q_to_Fp_init (num, den);
  modredc2ul2_intclear(num);
  modredc2ul2_intclear(den);

  for (unsigned long i = 0; i< iter; i++) {
    modredc2ul2_batch_Q_to_Fp (r, context, k, neg, p, P_LEN);
  }

  modredc2ul2_batch_Q_to_Fp_clear (context);

  tests_common_clear();
  exit (EXIT_SUCCESS);
}
