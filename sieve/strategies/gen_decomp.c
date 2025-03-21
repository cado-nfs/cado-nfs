#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include <math.h>
#include "decomp.h"
#include "gen_decomp.h"
#include "tab_decomp.h"

#define PREC 1000000

/* number of primes <= 2^n */
/* All values below fit in 53-bit mantissa, so that makes well-defined
 * floating point literals.
 */
double A7053[256] = {0,1,2,4,6,11,18,31,54,97,172,309,564,1028,1900,
                            3512,6542,12251,23000,43390,82025,155611,295947,
                            564163,1077871,2063689,3957809,7603553,14630843,
                            28192750,54400028,105097565,203280221,393615806,
                            762939111,1480206279,2874398515,5586502348,
                            10866266172,21151907950,41203088796,80316571436,
                            156661034233,305761713237,597116381732,
                            1166746786182,2280998753949,4461632979717,
                            8731188863470,17094432576778,33483379603407,
                            65612899915304,128625503610475,0,};

/* return prime_pi(2^i) */
static double
prime_pi (unsigned long i)
{
  if (A7053[i] != 0)
    return (double) A7053[i];
  else
    {
      double x = ldexp (1.0, i);
      return x / log (x);
    }
}

/* generate a random i-bit integer, distributed with density in 1/log(x). */
static double
gen_random (unsigned int i)
{
  double n0 = prime_pi (i - 1);
  double n1 = prime_pi (i);
  double a, b, n;

  /* we assume the n-th prime is in a*n*log(n)+b, thus we want:
     a*n0*log(n0) + b = 2^(i-1)
     a*n1*log(n1) + b = 2^i */
  a = ldexp (1.0, i - 1) / (n1 * log (n1) - n0 * log (n0));
  b = ldexp (1.0, i) - a * n1 * log (n1);
  n = n0 + (n1 - n0) * ((double)rand() / (RAND_MAX+1UL));
  return a * n * log (n) + b;
}


double T[256] = {0,};

/* estimates number of products of primes:
   l[0]-bit * l[1]-bit * ... * l[k-1]-bit
   with product of mfb bits, each prime should be >= lim */
static double
psi (unsigned int *l, unsigned int k, unsigned long mfb, unsigned long lim)
{
  unsigned long ok, tot;
  double S, r, rmin, rmax;

  ok = tot = 0;
  S = 1;
  for (unsigned int i = 0 ; i < k ; i++)
      S *= T[l[i]];

  if (S == 0.0)
    return 0;
  rmin = ldexp (1.0, (int) mfb - 1);
  rmax = ldexp (1.0, (int) mfb);
  while (tot < PREC)
    {
      r = 1;
      for (unsigned int i = 0 ; i < k ; i++)
        {
          double p = gen_random (l[i]);
          if (p >= (double) lim)
            r *= p;
          else
            r = 0.0;
        }
      tot ++;
      if (rmin <= r && r <= rmax)
        ok ++;
    }
  //printf ("%lu", mfb);
  for (unsigned int i = 0; i < k; i++)
    {
      unsigned int j;
      for (j = 1; i >= j && l[i-j] == l[i]; j++);
      S /= (double) j;
    }
  /* for (i = 0; i < k; i++) */
  /*   printf (" %d", l[i]); */
  /* printf (" %e # %d\n", (double) ok * S / (double) tot, ++count); */
  /* fflush (stdout); */
  return (double) ok * S / (double) tot;
}

tabular_decomp_t*
generate_all_decomp (unsigned int mfb, unsigned long lim)
{
  tabular_decomp_t* res = tabular_decomp_create();  
  unsigned int l[256];
  double p0, p1;

  unsigned int imin = ceil (log2 ((double) (lim + 1)));
  for (unsigned int i = 1; (i < 256) && i <= mfb; i++)
    {
      p0 = ldexp (1.0, (int) i - 1);
      p0 = (p0 < (double) lim) ? (double) lim / log ((double) lim)
        : prime_pi (i - 1);
      p1 = ldexp (1.0, (int) i);
      p1 = (p1 <= (double) lim) ? (double) lim / log ((double) lim)
        : prime_pi (i);
      T[i] = p1 - p0;
      // if (T[i] != 0) printf ("T[%lu]=%.0f\n", i, T[i]);
    }
  unsigned int kmax = floor (log (ldexp (1.0, (int) mfb)) / log ((double) lim));
  for (unsigned int k = 2; k <= kmax; k++)
    {
      /* a product l[0]-bit * l[1]-bit * ... * l[k-1]-bit can have
         l[0] + l[1] + ... + l[k-1] - (k-1) bits */
      for (unsigned int j = 0; j < k; j++)
        {
          for (unsigned int i = 0; i < k - 1; i++)
            l[i] = imin;
          l[k-1] = mfb + j - (k - 1) * imin;
	  if (l[k-1] < imin)
	    continue;
          while (1)
            {
              double val = psi (l, k, mfb, lim);
	      decomp_t* tmp = decomp_create (k, l, val);
	      tabular_decomp_add_decomp (res, tmp); 
	      decomp_free (tmp);

              /* next partition of mfb + j */
              unsigned int t;
              for (t = k - 1; t > 0; t--)
                /* try to increase l[t-1] and decrease l[t] */
                if (l[t-1] + 1 <= l[t])
                  {
                    l[t-1] ++;
                    unsigned int s = 1;
                    for (unsigned int u = t; u < k; u++)
                      {                  
                        s -= l[u];
                        l[u] = l[u-1];
                        s += l[u];
                      }
                    /* the total has increased by s */
                    l[k-1] -= s;
                    if (l[k-2] <= l[k-1])
                      break;
                  }
              if (t == 0)
                break;
            }
        }    
    }
  return res;
}
