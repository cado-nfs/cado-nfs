#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "modredc_ul.h"
#include "modredc_ul_default.h"
#include "roots_mod.h"
#include "macros.h"

void omega (residue_t o, residue_t b, const unsigned long k, const modulus_t pp);

int main(int argc, char **argv) {
  unsigned long p, minp = 3, maxp = 5000UL;
  enumeratediv_t div;
  int use_gp = 0;
  
  if (argc > 1 && strcmp(argv[1], "-gp") == 0){
    use_gp = 1;
    argc--;
    argv++;
  }
  
  if (argc > 1)
    minp = strtoul (argv[1], NULL, 10);
  if (argc > 2)
    maxp = strtoul (argv[2], NULL, 10);

#ifdef __COVERITY__
  __coverity_mark_pointee_as_sanitized(&maxp, LOOP_BOUND);
#endif

  printf ("/* minp = %lu, maxp = %lu */\n", minp, maxp);

  for (p = minp; p <= maxp; p += 2UL) {
    modulus_t pp;
    residue_t o, b, pow;
    unsigned long k;

    if (p <= 1 || (p > 3 && p % 3 == 0))
      continue;
    
    mod_initmod_ul (pp, p);
    
    if (!mod_isprime(pp))
      continue;
    
    mod_init (b, pp);
    mod_init (o, pp);
    mod_init (pow, pp);

    if (use_gp) {
      enumeratediv_init (&div, p-1);
      printf ("divisors(%lu) == vecsort([1", p-1);
      while ((k = enumeratediv(&div)) != 0)
        if (k > 1)
          printf (", %lu", k);
      printf ("])\n");
    }

    enumeratediv_init (&div, p-1);
    while ((k = enumeratediv(&div)) != 0) {
      omega (o, b, k, pp);

      if (use_gp) {
        /* Test using output to be piped into Pari/GP */
        printf ("znorder(Mod(%lu, %lu)) == %lu\n", 
                mod_get_ul(o, pp), mod_getmod_ul (pp), k);
      } else {
        /* Slow but reliable test */
        mod_pow_ul (pow, o, k, pp);
        ASSERT_ALWAYS (mod_is1(pow, pp));
        for (unsigned long i = 1; i < k; i++) {
          mod_pow_ul (pow, o, i, pp);
          ASSERT_ALWAYS (!mod_is1(pow, pp));
        }
      }
    }

    mod_clear (b, pp);
    mod_clear (o, pp);
    mod_clear (pow, pp);
    mod_clearmod (pp);
  }
  exit (EXIT_SUCCESS);
}
