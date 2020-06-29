/* Make SNFS polynomials for cyclotomic numbers.

  Copyright 2005, 2006, 2007, 2008 Alexander Kruppa.

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.
*/

/* Version: 0.1.2.2

History: 0.1 initial release
         0.1.1 Added GGNFS format, estimate of cost and time,
               changed x to mpz_t type.
	 0.1.2 Makes better polynomials for numbers with composite 
               base and no useful algebraic factors
         0.1.2.1 Allow reading N from stdin
         0.1.2.2 Allow printing polynomial in msieve format. Make better
                 polynomials if exponent is 3 (mod 6). Fixed sign in
                 printed number (i.e. was 2^128-1 instead of 2^128+1).

   Todo: Detect Aurifeullian factors and make polynomials for those.
         Compute root properties that depend on x.
	 Test some more to see if any of this actually works.
*/

#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gmp.h>

#define VERSION "0.1.2.1"

#define uint unsigned int

/* Output formats */
#define FRANKE 0
#define CWI 1
#define GGNFS 2
#define MSIEVE 3

static void
mpz_add_si (mpz_t R, mpz_t S, long i)
{
  if (i < 0)
    mpz_sub_ui (R, S, (unsigned long) -i);
  else
    mpz_add_ui (R, S, (unsigned long) i);
}


/* Compute p[0]^v[0] * ... * p[k-1]^v[k-1] */

static void
mpz_mul_list (mpz_t e, int *p, int *v, uint k)
{
    int i;
    uint j;

    mpz_set_ui (e, 1);
    for (j = 0; j < k; j++)
	for (i = 0; i < v[j]; i++)
	    mpz_mul_ui (e, e, p[j]);
}

/* Find prime factors < maxp of R. Cofactor is put in R, 
   the i-th prime factor in p_i with exponent v_i. 
   Returns number of distinct prime factors found. */

static uint
trialdiv (mpz_t R, int *p, int *v, uint maxp, uint flen)
{
    uint i, k = 0;

#ifdef DEBUG
    printf ("trialdiv: ");
#endif

    for (i = 2; i < maxp && k < flen; i = (i + 1) | 1) /* 2 3 5 7 9... */
	if (mpz_divisible_ui_p (R, i))
	{
	    p[k] = i;
	    v[k] = 0;
	    do
	    {
		v[k]++;
		mpz_divexact_ui (R, R, i);
	    } while (mpz_divisible_ui_p (R, i));
#ifdef DEBUG
	    printf ("%d^%d ", i, v[k]);
#endif
	    k++;
	}

#ifdef DEBUG
    printf ("\n");
#endif

    return k;
}

static void 
usage ()
{
  printf ("phi \n" VERSION);
  printf ("Usage: phi [options] [n x] N\n");
  printf ("Makes SNFS polynomial for factoring the primitive part Phi_n(x) of x^n-1\n");
  printf ("Cofactor N must divide x^n-1 if n is odd, or x^(n/2)+1 if n is even\n");
  printf ("If you do not specify x and n, the program will try to auto-detect them.\n");
  printf ("Options:\n");
  printf ("-deg4    Make a degree 4 polynomial\n");
  printf ("-deg5    Make a degree 5 polynomial\n");
  printf ("-deg6    Make a degree 6 polynomial\n");
  printf ("         By default chooses the degree that minimises the SNFS difficulty\n");
  printf ("-franke  Use Franke file format for polynomial\n");
  printf ("-ggnfs   Use GGNFS file format for polynomial\n");
  printf ("-cwi     Use CWI file format for polynomial\n");
  printf ("-msieve  Use msieve file format for polynomial\n");
  printf ("-short   Omit unneeded lines for Franke file format\n");
  printf ("-        Read N from stdin\n");
}

static int
autodet (int *x, uint *n, mpz_t N)
{
    const int flen = 32;
    const int maxp = 2000;    /* Limit for prime divisors */
    const int maxbase = 1000; /* Limit for bases to try */
    int p[flen], v[flen];     /* p are primes, v exponents */
    int b, i, j, k;
    mpz_t R; /* Remaining cofactor */
    mpz_t e; /* exponent */

    mpz_init (R);
    mpz_sub_ui (R, N, 1);
    mpz_init (e);

    k = trialdiv (R, p, v, maxp, flen);

    mpz_mul_list (e, p, v, k);
#ifdef DEBUG
    gmp_printf ("\ne = %Zd\n", e);
#endif

    for (b = 2; b <= maxbase; b++)
    {
	/* Skip perfect powers */
	mpz_set_ui (R, b);
	if (mpz_perfect_power_p (R))
	    continue;
#ifdef DEBUG
/*	printf ("Trying base %d\n", b); */
#endif
	mpz_powm (R, R, e, N);
	if (mpz_cmp_ui (R, 1) == 0)
	{
#ifdef DEBUG
	    printf ("Discovered! base = %d\n", b);
#endif
	    break;
	}
    }

    if (b > maxbase)
    {
	/* No suitable base found */
	mpz_clear (R);
	mpz_clear (e);
	return 0;
    }

    /* Find exponent */
    for (j = 0; j < k; j++)
    {
	while (v[j] > 0)
	{
	    /* See if lowering this exponent by one still satisfies ord|e */
	    v[j]--;
	    mpz_mul_list (e, p, v, k);
#ifdef DEBUG
	    gmp_printf ("Trying if ord | e = %Zd\n", e);
#endif
	    mpz_set_ui (R, b);
	    mpz_powm (R, R, e, N);
	    if (mpz_cmp_ui (R, 1) != 0)
	    {
		/* No, cannot decrease this exponent any further. */
		v[j]++;	   /* Undo decrease */
#ifdef DEBUG
		printf ("No, leaving exponent of %d at %d\n", p[j], v[j]);
#endif
		break;     /* Try next prime */
	    }
#ifdef DEBUG
	    else
		printf ("Yes, decreasing exponent of %d to %d\n", p[j], v[j]);
#endif

	}
    }
    
    mpz_mul_list (e, p, v, k);
    if (mpz_fits_uint_p (e))
    {
	if (x != NULL)
	    *x = b;
	if (n != NULL)
	    *n = mpz_get_ui (e);
	i = 1;
    }
    else
	i = 0;

    mpz_clear (R);
    mpz_clear (e);

    return i;
}

/* Difficulty must be in base 10 */
static double 
snfscost (const double difficulty)
{
    const double c = 1.5262856567; /* (32/9)^(1/3) */
    const double logd = difficulty * 2.3025850930; /* *log(10) */

    return exp (c * pow (logd, 1./3) * pow (log (logd), 2./3));
}


int 
main (int argc, char **argv)
{
  int format = FRANKE; /* Output format */
  mpz_t x;
  int n = 0;    /* We factor phi_n(x) */
  int k;               /* k|n */
  int degree = 0;
  int sign, i, j;
  int halved = 0;      /* Did we use the degree halving trick? */
  int shortform = 0;   /* Short output for Franke format? */
  int read_stdin = 0;  /* Read number from stdin? */
  mpz_t N;             /* The cofactor of the number to factor */
  mpz_t f[7];          /* Coefficients on algebraic side, max degree = 6 */
  mpz_t g[2];          /* Coefficients on rational side, max degree = 1 */
  mpz_t M;             /* Common root */
  mpz_t t;             /* A temporary */
  double dx;           /* x as a double */
  double difficulty;
  double skewness;
  double alpha;
  const double timefactor = 1. / 2.1e15;
  
  /* Parse command line parameters */
  while (argc > 1 && argv[1][0] == '-')
    {
      if (strcmp (argv[1], "-cwi") == 0)
        {
          format = CWI;
          argc--;
          argv++;
        }
      else if (strcmp (argv[1], "-franke") == 0)
        {
          format = FRANKE;
          argc--;
          argv++;
        }
      else if (strcmp (argv[1], "-ggnfs") == 0)
        {
          format = GGNFS;
          argc--;
          argv++;
        }
      else if (strcmp (argv[1], "-msieve") == 0)
        {
          format = MSIEVE;
          argc--;
          argv++;
        }
      else if (strcmp (argv[1], "-short") == 0)
        {
          shortform = 1;
          argc--;
          argv++;
        }
      else if (strcmp (argv[1], "-deg4") == 0)
        {
          degree = 4;
          argc--;
          argv++;
        }
      else if (strcmp (argv[1], "-deg5") == 0)
        {
          degree = 5;
          argc--;
          argv++;
        }
      else if (strcmp (argv[1], "-deg6") == 0)
        {
          degree = 6;
          argc--;
          argv++;
        }
      else if (strcmp (argv[1], "-") == 0)
        {
          read_stdin = 1;
          argc--;
          argv++;
        }
      else
        {
          fprintf (stderr, "Error, unknown command line option: %s\n", 
                   argv[1]);
          usage();
          exit (EXIT_FAILURE);
        }
    }

  mpz_init (N);
  mpz_init (x);
  
  if (argc + read_stdin == 2)
    {
      int xi = 0;
      if (read_stdin)
        mpz_inp_str (N, stdin, 10);
      else
        mpz_set_str (N, argv[1], 10);
      n = 0;
      if (autodet (&xi, (uint *) &n, N))
      {
	  mpz_set_ui (x, xi);
#ifdef DEBUG
	  printf ("Auto-detected Phi_%d(%d)\n", n, xi);
      }
      else
      {
	  printf ("Auto-detectedion failed\n");
#endif
      }
    }
  else if (argc + read_stdin > 3)
    {
      n = atoi (argv[1]);
      mpz_set_str (x, argv[2], 10);
      if (read_stdin)
        mpz_inp_str (N, stdin, 10);
      else
        mpz_set_str (N, argv[3], 10);
    }

  if (n <= 0 || mpz_sgn (x) <= 0 || mpz_sgn (N) <= 0)
    {
      gmp_printf ("Error, invalid parameters: n = %d, x = %Zd\n", n, x);
      usage();
      exit (EXIT_FAILURE);
    }
  
  mpz_init (t);
  mpz_init (g[0]);
  mpz_init (g[1]);
  for (i = 0; i < 7; i++)
    mpz_init (f[i]);
  mpz_init (M);

  sign = -1;
  /* Simplify even n, that is x^(n/2)+1 numbers */
  if (n % 2 == 0)
    {
      sign = +1;
      n /= 2;
    }
  /* Thus the number we factor is x^n+sign */

  /* Check that N | x^n+sign */
  mpz_powm_ui (t, x, n, N);
  mpz_add_si (t, t, sign);
  mpz_mod (t, t, N);
  if (mpz_sgn (t) != 0)
    {
      gmp_fprintf (stderr, "Error, N does not divide %Zd^%d%c1\n", 
	       x, n, sign < 0 ? '-' : '+');
      exit (EXIT_FAILURE);
    }
  
  /* Find out what kind of number we have (algebraic factors etc.) and
     choose polynomial */
  
  dx = mpz_get_d (x);
  
  if ((degree == 0 || degree == 4) && n % 15 == 0)
    {
      /* x^(15*k)+-1, remaining size 8/15 = 0.5333... */
      degree = 4;
      k = n / 15;
      mpz_set_si (f[4], 1);
      mpz_set_si (f[3], sign);
      mpz_set_si (f[2], -4);
      mpz_set_si (f[1], -sign * 4);
      mpz_set_si (f[0], 1);
      mpz_pow_ui (M, x, k);
      halved = 1;
      skewness = 1.;
      difficulty = log10 (dx) * 8. * k;
      alpha = 1.684;
    }
  else if ((degree == 0 || degree == 6) && n % 21 == 0)
    {
      /* x^(21*k)+-1, remaining size 12/21 = 0.571428... */
      degree = 6;
      k = n / 21;
      mpz_set_si (f[6], 1);
      mpz_set_si (f[5], sign);
      mpz_set_si (f[4], -6);
      mpz_set_si (f[3], -sign * 6);
      mpz_set_si (f[2], 8);
      mpz_set_si (f[1], sign * 8);
      mpz_set_si (f[0], 1);
      mpz_pow_ui (M, x, k);
      halved = 1;
      skewness = 1.;
      difficulty = log10 (dx) * 12. * k;
      alpha = 2.116;
      }
  else if ((degree == 0 || degree == 6) && n % 9 == 0)
    {
      /* x^(9*k)+-1 /  x^(3*k)+-1 = x^(6k) -+ x^(3k) + 1
       remaining size 2/3 = 0.666... */
      degree = 6;
      k = n / 9;
      mpz_set_ui (f[6], 1);
      mpz_set_si (f[3], -sign);
      mpz_set_ui (f[0], 1);
      mpz_pow_ui (M, x, k);
      skewness = 1.;
      difficulty = log10 (dx) * 6. * k;
      alpha = 1.621;
    }
  else if ((degree == 0 || degree == 6) && n % 9 == 3)
    {
      /* x^(9*k+3)+-1, remaining size 2/3 = 0.666... */
      degree = 6;
      k = (n - 3) / 9;
      mpz_mul (f[6], x, x);
      mpz_mul_si (f[3], x, -sign);
      mpz_set_ui (f[0], 1);
      mpz_pow_ui (M, x, k);
      skewness = pow (dx, -1./3.);
      difficulty = log10 (dx) * (n / 3 * 2);
      alpha = 0.;
    }
  else if ((degree == 0 || degree == 6) && n % 9 == 6)
    {
      /* x^(9*k+3)+-1, remaining size 2/3 = 0.666... */
      degree = 6;
      k = (n + 3) / 9;
      mpz_set_ui (f[6], 1);
      mpz_mul_si (f[3], x, -sign);
      mpz_mul (f[0], x, x);
      mpz_pow_ui (M, x, k);
      skewness = pow (dx, 1./3.);
      difficulty = log10 (dx) * (n / 3 * 2 + 2);
      alpha = 0.;
    }
  else if ((degree == 0 || degree == 4) && n % 6 == 0)
    {
      /* x^(6*k)+1, remaining size 2/3 = 0.666... */
      degree = 4;
      k = n / 6;
      mpz_set_si (f[4], 1);
      mpz_set_si (f[2], -1);
      mpz_set_si (f[0], 1);
      mpz_pow_ui (M, x, k);
      skewness = 1.;
      difficulty = log10 (dx) * 4. * k;
      alpha = 1.292;
    }
  else if ((degree == 0 || degree == 4) && n % 6 == 3)
    {
      /* x^(6*k+3)+-1, remaining size 2/3 = 0.666... */
      degree = 4;
      k = (n - 3) / 6;
      mpz_mul (f[4], x, x);
      mpz_mul_si (f[2], x, -sign);
      mpz_set_si (f[0], 1);
      mpz_pow_ui (M, x, k);
      skewness = pow (dx, -0.5);
      difficulty = log10 (dx) * 4. * k;
      alpha = 0.;
    }
  else if ((degree == 0 || degree == 4) && n % 5 == 0)
    {
      /* x^(5*k)+-1, remaining size 4/5 = 0.8 */
      degree = 4;
      k = n / 5;
      mpz_set_si (f[4], 1);
      mpz_set_si (f[3], -sign);
      mpz_set_si (f[2], 1);
      mpz_set_si (f[1], -sign);
      mpz_set_si (f[0], 1);
      mpz_pow_ui (M, x, k);
      skewness = 1.;
      difficulty = log10 (dx) * 4. * k;
      alpha = 1.453;
    }
  else if ((degree == 0 || degree == 6) && n % 7 == 0)
    {
      /* x^(7*k)+-1, remaining size 6/7 = 0.857142... */
      degree = 6;
      k = n / 7;
      mpz_set_si (f[6], 1);
      mpz_set_si (f[5], -sign);
      mpz_set_si (f[4], 1);
      mpz_set_si (f[3], -sign);
      mpz_set_si (f[2], 1);
      mpz_set_si (f[1], -sign);
      mpz_set_si (f[0], 1);
      mpz_pow_ui (M, x, k);
      skewness = 1.;
      difficulty = log10 (dx) * 6. * k;
      alpha = 2.244;
    }
  else if ((degree == 0 || degree == 5) && n % 11 == 0)
    {
      /* x^(11*k)+-1, remaining size 10/11 = 0.90... */
      degree = 5;
      k = n / 11;
      mpz_set_si (f[5], 1);
      mpz_set_si (f[4], -sign*1);
      mpz_set_si (f[3], -4);
      mpz_set_si (f[2], sign*3);
      mpz_set_si (f[1], 3);
      mpz_set_si (f[0], -sign);
      mpz_pow_ui (M, x, k);
      halved = 1;
      skewness = 1.;
      difficulty = log10 (dx) * 10. * k;
      alpha = 2.224;
    }
  else if ((degree == 0 || degree == 6) && n % 13 == 0)
    {
      /* x^(13*k)+-1, remaining size 12/13 = 0.923076... */
      degree = 6;
      k = n / 13;
      mpz_set_si (f[6], 1);
      mpz_set_si (f[5], -sign);
      mpz_set_si (f[4], -5);
      mpz_set_si (f[3], sign*4);
      mpz_set_si (f[2], 6);
      mpz_set_si (f[1], -sign*3);
      mpz_set_si (f[0], -1);
      mpz_pow_ui (M, x, k);
      halved = 1;
      skewness = 1.;
      difficulty = log10 (dx) * 12. * k;
      alpha = 3.095;
    }
  else
    {
      /* Make polynomial without using any algebraic factor */
      const int d = (degree == 0) ? 5 : degree;
      const int maxp = 10000;
      const int flen = 20; /* Only examine 21 different factors 
			      (including cofactor). Should easily suffice */
      int p[flen], vp[flen], best = 0;
      mpz_t R;
      double rating, bestrating = 0.;

      mpz_init_set (R, x);
      k = trialdiv (R, p, vp, maxp, flen);
      /* Now x = p[0] ^ v[0] * ... * p[k - 1] ^ v[k - 1] * R */
      /* Here, k is the number of distinct prime factors */

      for (i = 0; i < k; i++)
	  vp[i] *= n;
      /* Now N = p[0] ^ v[0] * ... * p[k - 1] ^ v[k - 1] * R^n - 1 */

      /* Make all exponents multiples of $d$. 
	 When decreasing the exponent of $p$ by $u$, we need to multiply 
	   $f_d$ by $p^u$.
	 When increasing the exponent of $p$ by $u$, we need to multiply 
	   $f_0$ by $p^u$ and $M$ by $p$.
	 For each exponent $v$, we can choose between increasing by 
	   $(-v) \bmod 5$ or decreasing by $v \bmod 5$.
	 So we have a search space of $2^{k+1}$ possibilities.

	 For now, we simply use $f_d + f_0 + c$ as the function to minimise,
	   where c is what we multiplied M by.
	 This rating function could and should be improved, but it seems to 
	   usually produce decent polynomials.
      */

      for (j = 0; j < (1 << (k + ((mpz_cmp_ui (R, 1) > 1) ? 1 : 0))); j++)
      {
	  double fd = 1., f0 = 1., c = 1.;
	  for (i = 0; i < k; i++)
	  {
	      if (j & (1 << i))
	      {
		  /* Reduce exponent by v[i] % d */
		  fd *= pow ((double)p[i], (double)(vp[i] % d));
	      }
	      else if (vp[i] % d != 0)
	      {
		  /* Multiply by -v[i] % d */
		  f0 *= pow ((double)p[i], (double)(d - vp[i] % d));
		  c *= p[i];
	      }
	  }
	  if (j & (1 << k))
	      fd *= pow (mpz_get_d (R), (double)(n % d));
	  else
	  {
	      f0 *= pow (mpz_get_d (R), (double)((d - n % d) % d));
	      c *= mpz_get_d (R);
	  }
	  /* rating = c * c * pow (fd, 1. + 1./d) * pow (f0, 1. - 1./d); */
	  rating = fd + f0 + c;
#ifdef DEBUG
	  printf ("j = %d, %.0f * x^%d - %.0f, c = %.0f, rating = %f\n", 
		  j, fd, d, f0, c, rating);
#endif
	  if (j == 0 || rating < bestrating)
	  {
	      best = j;
	      bestrating = rating;
#ifdef DEBUG
	      printf ("New best rating!\n");
#endif
	  }

      }

      degree = d;
      difficulty = log10 (dx) * n;
      mpz_set_ui (M, 1);

      /* Now make the actual coefficients f_d and f_0 */

      mpz_set_ui (f[d], 1);
      mpz_set_ui (f[0], 1);
      for (i = 0; i < k; i++)
      {
	  if (best & (1 << i))
	  {
	      uint u = vp[i] % d;
	      mpz_ui_pow_ui (t, p[i], u);
	      mpz_mul (f[d], f[d], t);
	      mpz_ui_pow_ui (t, p[i], (vp[i] - u) / d);
	      mpz_mul (M, M, t);
	  }
	  else
	  {
	      uint u = (d - vp[i] % d) % d; /* u = (-vp[i]) % d */
	      mpz_ui_pow_ui (t, p[i], u);
	      mpz_mul (f[0], f[0], t);
	      difficulty += log10 ((double) p[i]) * u;
	      mpz_ui_pow_ui (t, p[i], (vp[i] + u) / d);
	      mpz_mul (M, M, t);
	  }
      }
      if (best & (1 << k))
      {
	  uint u = n % d;
	  mpz_pow_ui (t, R, u);
	  mpz_mul (f[d], f[d], t);
	  mpz_pow_ui (t, R, (n - u) / d);
	  mpz_mul (M, M, t);
      }
      else 
      {
	  uint u = (d - n % d) % d;
	  mpz_pow_ui (t, R, u);
	  mpz_mul (f[0], f[0], t);
	  mpz_pow_ui (t, R, (n + u) / d);
	  mpz_mul (M, M, t);
	  difficulty += log10 (mpz_get_d (R)) * u;
      }

      skewness = pow (mpz_get_d (f[0]) / mpz_get_d (f[d]), 1. / d);

      if (sign == -1)
        mpz_neg (f[0], f[0]);
      alpha = 0.; /* Alpha may depend on x, can't calculate yet */
    }

  if (halved)
    {
      mpz_neg (g[1], M);          /* Y1 = -x^k */
      mpz_mul (g[0], M, M);
      mpz_add_ui (g[0], g[0], 1); /* Y0 = x^(2k)+1 */
      mpz_invert (t, M, N);
      mpz_add (M, M, t);          /* M = x^k + x^(-k) */
    }
  else
    {
      mpz_set (g[0], M);
      mpz_set_si (g[1], -1);
    }

  /* Test polynomials (check that values vanish (mod N) at root M) */
  mpz_set_ui (t, 0);
  for (i = degree; i >= 0; i--)
    {
      mpz_mul (t, t, M);
      mpz_add (t, t, f[i]);
      mpz_mod (t, t, N);
    }
  if (mpz_sgn (t) != 0)
    {
      gmp_fprintf (stderr, "Error: M=%Zd is not a root of f(x) % N\n", M);
      fprintf (stderr, "f(x) = ");
      for (i = degree; i >= 0; i--)
	  gmp_fprintf (stderr, "%Zd*x^%d +", f[i], i);
      gmp_fprintf (stderr, "\n""Remainder is %Zd\n", t);
      gmp_fprintf (stderr, "%Zd^%d%c1\n", x, n, sign < 0 ? '-' : '+');
      fprintf (stderr, "Please report this bug.\n");
      exit (EXIT_FAILURE);
    }

  mpz_mul (t, g[1], M);
  mpz_add (t, t, g[0]);
  mpz_mod (t, t, N);
  if (mpz_sgn (t) != 0)
    {
      gmp_fprintf (stderr, "Error: M=%Zd is not a root of g(x) % N\n", M);
      gmp_fprintf (stderr, "Remainder is %Zd\n", t);
      fprintf (stderr, "Please report this bug.\n");
      exit (EXIT_FAILURE);
    }

  /* Output polynomials */
  if (format == FRANKE)
  {
      const double cost = snfscost (difficulty + alpha / log (10.) * degree);

      gmp_printf ("%Zd\n", N);
      gmp_printf ("# %Zd^%d%c1, difficulty: %.2f, skewness: %.2f, alpha: %.2f\n"
		  "# Cost: %g, est. time: %.2f GHz days (not accurate yet!)\n",
		  x, n, sign == 1 ? '+' : '-', difficulty, skewness, alpha, 
		  cost, cost * timefactor);

      for (i = degree; i >= 0; i--)
        if (!shortform || mpz_sgn (f[i]) != 0)
          gmp_printf ("X%d %Zd\n", i, f[i]);

      if (!shortform || halved) /* Franke's siever uses g(x)=-x+M by default */
        {
          gmp_printf ("Y1 %Zd\n", g[1]);
          gmp_printf ("Y0 %Zd\n", g[0]);
        }
      gmp_printf ("M %Zd\n", M);
    }
  else if (format == GGNFS)
    {
      const double cost = snfscost (difficulty);

      gmp_printf ("n: %Zd\n", N);
      gmp_printf ("# %Zd^%d%c1, difficulty: %.2f, skewness: %.2f, alpha: %.2f\n"
		  "# cost: %g, est. time: %.2f GHz days (not accurate yet!)\n",
		  x, n, sign == 1 ? '+' : '-', difficulty, skewness, alpha, 
		  cost, cost * timefactor);

      printf ("skew: %.3f\n", skewness);
      for (i = degree; i >= 0; i--)
        if (mpz_sgn (f[i]) != 0)
          gmp_printf ("c%d: %Zd\n", i, f[i]);

      gmp_printf ("Y1: %Zd\n", g[1]);
      gmp_printf ("Y0: %Zd\n", g[0]);
      gmp_printf ("m: %Zd\n", M);
      printf ("type: snfs\n");
    }
  else if (format == MSIEVE)
    {
      gmp_printf ("N %Zd\n", N);
      gmp_printf ("R0 %Zd\n", g[0]);
      gmp_printf ("R1 %Zd\n", g[1]);
      for (i = 0; i <= degree; i++)
        gmp_printf ("A%d %Zd\n", i, f[i]);
    }
  else
    {
      gmp_printf ("%Zd\n", N);
      gmp_printf ("%Zd\n\n2\n\n1\n", M);
      gmp_printf ("%Zd %Zd\n", g[0], g[1]);
      printf ("\n\n%d\n", degree);
      for (i = 0; i <= degree; i++)
        gmp_printf ("%Zd ", f[i]);
      printf ("\n");
    }

  return 0;
}
