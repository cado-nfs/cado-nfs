#include "cado.h" // IWYU pragma: keep

#include <stdint.h>
#include <stdio.h>

#include "arith/modredc_ul.h"
#include "arith/modredc_15ul.h"
#include "arith/modredc_2ul2.h"
#include "timing.h"

#define NBMUL_1UL  1000000000
#define NBMUL_15UL 1000000000
#define NBMUL_2UL  1000000000


double test_sqr_modredcul (modulusredcul_t m)
{
  uint64_t starttime, endtime;
  residueredcul_t x;

  modredcul_init (x, m);
  /* random 64-bit integer */
  modredcul_set_ul (x, 13231468506376890262UL, m);
  printf ("Running %d modredcul_sqr (x, x, m)...\n", NBMUL_1UL);

  starttime = microseconds();
  for (unsigned long i = 0; i < NBMUL_1UL ; i++)
    {
      modredcul_sqr (x, x, m);
    }
  endtime = microseconds();

  fprintf (stderr, "x = %lu\n", x[0]);
  modredcul_clear (x, m);

  return endtime - starttime;
}


double test_mul_modredcul (modulusredcul_t m)
{
  uint64_t starttime, endtime;
  residueredcul_t x, y;

  modredcul_init (x, m);
  modredcul_init (y, m);
  /* random 64-bit integers */
  modredcul_set_ul (x, 13231468506376890262UL, m);
  modredcul_set_ul (y, 2217200224969896255UL, m);

  printf ("Running %d modredcul_mul(x, x, y, m)...\n", NBMUL_1UL);

  starttime = microseconds();
  for (unsigned long i = 0; i < NBMUL_1UL ; i++)
    {
      modredcul_mul (x, x, y, m);
    }
  endtime = microseconds();
  
  fprintf (stderr, "x = %lu\n", x[0]);
  modredcul_clear (x, m);
  modredcul_clear (y, m);

  return endtime - starttime;
}

/* ------------------------------------------------------------------------- */

double test_sqr_modredc15ul (modulusredc15ul_t m)
{
  uint64_t starttime, endtime;
  modintredc15ul_t x;
  modredc15ul_init (x, m);
  /* random 96-bit integer */
  unsigned long t[2] = {10453954462211070795UL, 959555781UL};
  modredc15ul_intset_uls (x, t, 2);

  printf ("Running %d modredc15ul_sqr (x, x, m)...\n", NBMUL_15UL);

  starttime = microseconds();
  for (unsigned long i = 0; i < NBMUL_15UL ; i++)
    {
      modredc15ul_sqr (x, x, m);
    }
  endtime = microseconds();
  
  fprintf (stderr, "x[0], x[1] = %lu, %lu\n", x[0], x[1]);
  modredc15ul_clear (x, m);

  return endtime - starttime;
}

double test_mul_modredc15ul (modulusredc15ul_t m)
{
  uint64_t starttime, endtime;
  residueredc15ul_t x, y;
  modredc15ul_init (x, m);
  modredc15ul_init (y, m);

  /* random 96-bit integers */
  unsigned long tx[2] = {11737700200122420177UL, 4135799180UL};
  unsigned long ty[2] = {14294739848942950433UL, 1799857174UL};
  modredc15ul_intset_uls (x, tx, 2);
  modredc15ul_intset_uls (y, ty, 2);

  printf ("Running %d modredc15ul_mul(x, x, y, m)...\n", NBMUL_15UL);

  starttime = microseconds();
  for (unsigned long i = 0; i < NBMUL_15UL ; i++)
    {
      modredc15ul_mul (x, x, y, m);
    }
  endtime = microseconds();

  fprintf (stderr, "x[0], x[1] = %lu, %lu\n", x[0], x[1]);
  modredc15ul_clear (x, m);
  modredc15ul_clear (y, m);

  return endtime - starttime;
}

/* ------------------------------------------------------------------------- */

double test_sqr_modredc2ul (modulusredc2ul2_t m)
{
  uint64_t starttime, endtime;
  modintredc2ul2_t x;
  modredc2ul2_init (x, m);
  
  /* random 126-bit integer */
  unsigned long t[2] = {9924822112879250552UL, 14313105974994231104UL};
  modredc2ul2_intset_uls (x, t, 2);

  printf ("Running %d modredc2ul_sqr (x, x, m)...\n", NBMUL_2UL);

  starttime = microseconds();
  for (unsigned long i = 0; i < NBMUL_2UL ; i++)
    {
      modredc2ul2_sqr (x, x, m);
    }
  endtime = microseconds();

  fprintf (stderr, "x[0], x[1] = %lu, %lu\n", x[0], x[1]);
  modredc2ul2_clear (x, m);

  return endtime - starttime;
}

double test_mul_modredc2ul (modulusredc2ul2_t m)
{
  uint64_t starttime, endtime;
  residueredc2ul2_t x, y;
  modredc2ul2_init (x, m);
  modredc2ul2_init (y, m);

  /* random 126-bit integers */
  unsigned long tx[2] = {17311386395127023723UL, 3078708599844013985UL};
  unsigned long ty[2] = {1430637974752199232UL, 8694567303182869329UL};
  modredc2ul2_intset_uls (x, tx, 2);
  modredc2ul2_intset_uls (y, ty, 2);

  printf ("Running %d modredc2ul_mul(x, x, y, m)...\n", NBMUL_2UL);

  starttime = microseconds();
  for (unsigned long i = 0; i < NBMUL_2UL ; i++)
    {
      modredc2ul2_mul (x, x, y, m);
    }
  endtime = microseconds();

  fprintf (stderr, "x[0], x[1] = %lu, %lu\n", x[0], x[1]);
  modredc2ul2_clear (x, m);
  modredc2ul2_clear (y, m);

  return endtime - starttime;
}


/* ------------------------------------------------------------------------- */

int main()
{
  uint64_t starttime, endtime;
  double usrtime;

  /* ul */

  modulusredcul_t m64;

  /* largest prime less than 2^64 */
  modredcul_initmod_ul (m64, 18446744073709551557UL);
  printf ("m = %lu\n", m64->m);
  
  usrtime = test_sqr_modredcul (m64);
  printf ("Total time: %.2f s\n", usrtime / 1000000);

  usrtime = test_mul_modredcul (m64);
  printf ("Total time: %.2f s\n", usrtime / 1000000);
  
  modredcul_clearmod (m64);

  printf ("========================================\n");
  
  /* 15ul */

  modintredc15ul_t t;
  modulusredc15ul_t m96;

  /* largest prime less than 2^96 */
  unsigned long s[2] = {18446744073709551599UL, 4294967295UL};
  modredc15ul_intset_uls (t, s, 2);
  modredc15ul_initmod_int (m96, t);
  
  printf ("m0, m1 = %lu, %lu\n", m96->m[0], m96->m[1]);
  
  starttime = microseconds();
  test_sqr_modredc15ul (m96);
  endtime = microseconds();

  usrtime = endtime - starttime;
  printf ("Total time: %.2f s\n", usrtime / 1000000);

  usrtime = test_mul_modredc15ul (m96);
  printf ("Total time: %.2f s\n", usrtime / 1000000);
  
  modredc15ul_clearmod (m96);

  printf ("========================================\n");
  
  /* 2ul */

  modintredc2ul2_t tt;
  modulusredc2ul2_t m126;

  /* largest prime less than 2^126 */
  unsigned long ss[2] = {18446744073709551479UL, 4611686018427387903UL};

  modredc2ul2_intset_uls (tt, ss, 2);
  modredc2ul2_initmod_int (m126, tt);
  
  printf ("m0, m1 = %lu, %lu\n", m126->m[0], m126->m[1]);
  
  usrtime = test_sqr_modredc2ul (m126);
  printf ("Total time: %.2f s\n", usrtime / 1000000);

  usrtime = test_mul_modredc2ul (m126);
  printf ("Total time: %.2f s\n", usrtime / 1000000);
  
  modredc2ul2_clearmod (m126);

  return 0;
}
