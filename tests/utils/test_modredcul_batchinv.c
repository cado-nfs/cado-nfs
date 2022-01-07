#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include <stdio.h>
#include "tests_common.h"
#include "modredc_ul.h"
#include "misc.h"

static int verbose = 0;

void
print_residues(const char*prefix, residueredcul_t *a, size_t len,
               const modulusredcul_t m)
{
  if (!verbose)
    return;
  printf("%s: ", prefix);
  for(size_t i = 0; i < len; i++) {
    printf("Mod(%lu,%lu)%s", 
           modredcul_get_ul(a[i], m), modredcul_getmod_ul(m),
           i + 1 < len ? ", " : "");
  }
  printf("\n");
}

/* Returns 1 if a1[i]*a2[i] == 1, for 0 <= i < len,
   returns 0 otherwise. */
int
check_product_1(residueredcul_t *a1, residueredcul_t *a2,
                const size_t len, const modulusredcul_t m)
{
   residueredcul_t p;
   int ok = 1;
   
   modredcul_init(p, m);
   for (size_t i = 0; i < len; i++) {
     modredcul_mul(p, a1[i], a2[i], m);
     if(!modredcul_is1(p, m)) {
       ok = 0;
       fprintf (stderr, "Inverse is incorrect: %lu * %lu = %lu != 1 (mod %lu)\n",
                modredcul_get_ul(a1[i], m),
                modredcul_get_ul(a2[i], m),
                modredcul_get_ul(p, m),
                modredcul_getmod_ul(m));
       break;
     }
   }
   modredcul_clear(p, m);
   return ok;
}

int
check_product_c(residueredcul_t *a1, residueredcul_t *a2,
                const size_t len, const residueredcul_t c,
                const modulusredcul_t m)
{
   residueredcul_t p;
   int ok = 1;
   
   modredcul_init(p, m);
   for (size_t i = 0; i < len; i++) {
     modredcul_mul(p, a1[i], a2[i], m);
     if(!modredcul_equal(p, c, m)) {
       ok = 0;
       fprintf (stderr, "Inverse is incorrect: %lu * %lu = %lu != %lu (mod %lu)\n",
                modredcul_get_ul(a1[i], m),
                modredcul_get_ul(a2[i], m),
                modredcul_get_ul(p, m),
                modredcul_get_ul(c, m),
                modredcul_getmod_ul(m));
       break;
     }
   }
   modredcul_clear(p, m);
   return ok;
}

int
test_modredc_batchinv (const size_t len, unsigned long uc)
{
  residueredcul_t *a, *r, c;
  modulusredcul_t m;
  int ok = 1;
  
  /* Random, odd modulus */
  modredcul_initmod_ul(m, u64_random(state) | 1);
  
  a = (residueredcul_t *) malloc(len * sizeof(residueredcul_t));
  r = (residueredcul_t *) malloc(len * sizeof(residueredcul_t));
  
  for (size_t i = 0; i < len; i++) {
    modredcul_init(a[i], m);
    modredcul_set_ul(a[i], u64_random(state), m);
    if (modredcul_is0(a[i], m)) {
      modredcul_set1(a[i], m);
    }
    modredcul_init(r[i], m);
  }
  modredcul_init(c, m);
  modredcul_set_ul(c, uc, m);
  
  print_residues("Input numbers", a, len, m);
  int rc = modredcul_batchinv(r, (const residueredcul_t *) a, len, 
                              uc == 0 ? NULL : c, m);
  if (rc != 0)
    print_residues("Inverses", r, len, m);
  
  /* Check that product of input and output is 1 */
  if (rc == 1) {
    if (uc == 0) {
      if (!check_product_1(r, a, len, m)) {
        fprintf (stderr, "Product of input and its inverse it not 1\n");
        ok = 0;
        goto finish;
      }
    } else {
      if (!check_product_c(r, a, len, c, m)) {
        fprintf (stderr, "Product of input and its inverse it not %lu\n", uc);
        ok = 0;
        goto finish;
      }
    }
  }
  
  int had_common_factor = 0;
  modintredcul_t g, q;
  modredcul_intinit(g);
  modredcul_intinit(q);
  for (size_t i = 0; i < len; i++) {
    modredcul_gcd(g, a[i], m);
    /* May have to divide out factors with greater multiplicity in a than
       they have in m */
    while (!modredcul_intequal_ul(g, 1UL)) {
      had_common_factor = 1;
      /* Divide out gcd so we can try again with a set of residues for which
         the inverses exist */
      modredcul_get_int(q,  a[i], m);
      modredcul_intdivexact(q, q, g);
      modredcul_set_int(a[i], q, m);
      modredcul_gcd(g, a[i], m);
    }
  }
  modredcul_intclear(g);
  modredcul_intclear(q);

  if (had_common_factor && rc != 0) {
    fprintf(stderr, "Input numbers had common factor but modredcul_batchinv()"
            " returned %d", rc);
    ok = 0;
    goto finish;
  }

  print_residues("Input numbers after dividing by gcd", a, len, m);
  rc = modredcul_batchinv(r, (const residueredcul_t *) a, len, NULL, m);
  if (rc == 0) {
    fprintf (stderr, "After dividing out common factors, modredcul_batchinv()"
             " still returned 0\n");
    ok = 0;
    goto finish;
  }
  print_residues("Inverses", r, len, m);
  
  /* Check that product of input and output is 1 */
  if (!check_product_1(r, a, len, m)) {
    fprintf (stderr, "Product of input and its inverse it not 1\n");
    ok = 0;
    goto finish;
  }
  
  finish:
  for (size_t i = 0; i < len; i++) {
    modredcul_clear(a[i], m);
    modredcul_clear(r[i], m);
  }
  modredcul_clear(c, m);
  free(a);
  free(r);
  return ok;
}

int
main (int argc, const char *argv[])
{
  unsigned long iter = 20;
  int ok = 1;
  tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_ITER | PARSE_VERBOSE);
  tests_common_get_iter (&iter);
  verbose = tests_common_get_verbose();
  
  for (unsigned long i = 0; ok && i < iter; i++) {
    ok = test_modredc_batchinv(i, 0);
    ok = test_modredc_batchinv(i, 1);
    ok = test_modredc_batchinv(i, u64_random(state));
  }
  
  tests_common_clear();
  if (ok)
    exit (EXIT_SUCCESS);
  else
    exit (EXIT_FAILURE);
}
