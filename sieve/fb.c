/*****************************************************************
 *                Functions for the factor base                  *
 *****************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#define rdtscll(x)
#include "basicnt.h"
#include "fb.h"
#include "utils.h"


void 
fb_fprint_entry (FILE *fd, const factorbase_degn_t *fb)
{
  int i;
  fprintf (fd, "Prime " FBPRIME_FORMAT " with rounded log %d and roots ", 
	   fb->p, (int) fb->plog);
  for (i = 0; i < fb->nr_roots; i++)
    {
      fprintf (fd, FBROOT_FORMAT, fb->roots[i]);
      if (i + 1 < fb->nr_roots)
	fprintf (fd, ", ");
    }
  fprintf (fd, "\n");
}

void 
fb_fprint (FILE *fd, const factorbase_degn_t *fb)
{
  while (fb->p != FB_END)
    {
      fb_fprint_entry (fd, fb);
      fb = fb_next (fb);
    }
}

/* Add fb_add to (void *)fb + fbsize. If a realloc failed, returns NULL.
   fb_add->size need not be set by caller, this function does it */

static factorbase_degn_t *
fb_add_to (factorbase_degn_t *fb, size_t *fbsize, size_t *fballoc,
	   const size_t allocblocksize, factorbase_degn_t *fb_add)
{
  const size_t fb_addsize  = fb_entrysize (fb_add); 
  factorbase_degn_t *newfb = fb;

  ASSERT(fb_addsize <= allocblocksize); /* Otherwise we still might not have
					   enough mem after the realloc */

  /* Do we need more memory for fb? */
  if (*fballoc < *fbsize + fb_addsize)
    {
      *fballoc += allocblocksize;
      newfb = (factorbase_degn_t *) realloc (fb, *fballoc);
      if (newfb == NULL)
	{
	  fprintf (stderr, 
		   "Could not reallocate factor base to %lu bytes\n",
		   (unsigned long) *fballoc);
	  return NULL;
	}
    }
  memcpy (fb_skip(newfb, *fbsize), fb_add, fb_addsize);
  fb_skip(newfb, *fbsize)->size = fb_addsize;
  
  *fbsize += fb_addsize;

  return newfb;
}

/* Sort n primes in array *primes into ascending order */

void
fb_sortprimes (fbprime_t *primes, const unsigned int n)
{
  unsigned int k, l, m;
  fbprime_t t;

  for (l = n; l > 1; l = m)
    for (k = m = 0; k < l - 1; k++)
      if (primes[k] > primes[k + 1])
	{
	  t = primes[k];
	  primes[k] = primes[k + 1];
	  primes[k + 1] = t;
	  m = k + 1;
	}
#ifndef NDEBUG
  for (k = 1; k < n; k++)
    ASSERT(primes[k - 1] <= primes[k]);
#endif
}

#if 0
static double
fb_log_double (double n, double log_scale, double offset)
{
  return floor (log (n) * log_scale + offset + 0.5);
}
#endif

unsigned char
fb_log (double n, double log_scale, double offset)
{
  return (unsigned char) floor (log (n) * log_scale + offset + 0.5);
}

/* Assume q is a prime power p^k with k>=1, return p if k > 1, 0 otherwise. */
static uint32_t
is_prime_power (uint32_t q)
{
  unsigned int maxk, k;
  uint32_t p;

  for (maxk = 0, p = q; p > 1; p /= 2, maxk ++);
  for (k = maxk; k >= 2; k--)
    {
      p = (uint32_t) (pow ((double) q, 1.0 / (double) k) + 0.5);
      if (q % p == 0)
        return p;
    }
  return 0;
}


/* Make one factor base entry for a linear polynomial poly[1] * x + poly[0]
   and the prime (power) q. We assume that poly[0] and poly[1] are coprime.
   Non-projective roots a/b such that poly[1] * a + poly[0] * b == 0 (mod q) 
   with gcd(poly[1], q) = 1 are stored as a/b mod q.
   If do_projective != 0, also stores projective roots with gcd(q, f_1) > 1, 
   but stores the reciprocal root plus p, p + b/a mod p^k.
   If no root was stored (projective root with do_projective = 0) returns 0, 
   returns 1 otherwise. */

static int
fb_make_linear1 (factorbase_degn_t *fb_entry, const mpz_t *poly,
		 const fbprime_t p, const unsigned char logp, 
		 const int do_projective)
{
  modulusul_t m;
  residueul_t r0, r1;
  modintul_t gcd;
  int is_projective, rc;
  
  fb_entry->p = p;
  fb_entry->plog = logp;
  
  modul_initmod_ul (m, p);
  modul_init_noset0 (r0, m);
  modul_init_noset0 (r1, m);
  
  modul_set_ul_reduced (r0, mpz_fdiv_ui (poly[0], p), m);
  modul_set_ul_reduced (r1, mpz_fdiv_ui (poly[1], p), m);
  modul_gcd (gcd, r1, m);
  is_projective = (modul_intcmp_ul (gcd, 1UL) > 0);
  
  if (is_projective)
    {
      if (!do_projective)
	{
	  modul_clear (r0, m);
	  modul_clear (r1, m);
	  modul_clearmod (m);
	  return 0; /* If we don't do projective roots, just return */
	}
      modul_swap (r0, r1, m);
    }
  
  /* We want poly[1] * a + poly[0] * b == 0 <=> 
     a/b == - poly[0] / poly[1] */
  rc = modul_inv (r1, r1, m); /* r1 = 1 / poly[1] */
  ASSERT_ALWAYS (rc != 0);
  modul_mul (r1, r0, r1, m); /* r1 = poly[0] / poly[1] */
  modul_neg (r1, r1, m); /* r1 = - poly[0] / poly[1] */
  
  fb_entry->roots[0] = modul_get_ul (r1, m) + (is_projective ? p : 0);
  if (p % 2 != 0)
    fb_entry->invp = - modul_invmodlong (modul_getmod_ul (m));
  
  modul_clear (r0, m);
  modul_clear (r1, m);
  modul_clearmod (m);
  
  return 1;
}


/* Generate a factor base with primes <= bound and prime powers <= powbound 
   for a linear polynomial. If projective != 0, adds projective roots 
   (for primes that divide leading coefficient) */

factorbase_degn_t *
fb_make_linear (const mpz_t *poly, const fbprime_t bound, 
		const fbprime_t powbound, const double log_scale, 
		const int verbose, const int projective, FILE *output)
{
  fbprime_t p;
  factorbase_degn_t *fb = NULL, *fb_cur, *fb_new;
  size_t fbsize = 0, fballoc = 0, pow_len = 0;
  const size_t allocblocksize = 1 << 20;
  unsigned long logp;
  int had_proj_root = 0;
  fbprime_t *powers = NULL, min_pow = 0; /* List of prime powers that yet 
					    need to be included, and the 
					    minimum among them */

  rdtscll (tsc1);
  fb_cur = (factorbase_degn_t *) malloc (fb_entrysize_uc (1));
  ASSERT (fb_cur != NULL);

  fb_cur->nr_roots = 1;
  fb_cur->size = fb_entrysize_uc (1);

  if (verbose)
    gmp_fprintf (output, 
		 "# Making factor base for polynomial g(x) = %Zd * x + %Zd,\n"
		 "# including primes up to " FBPRIME_FORMAT 
		 " and prime powers up to " FBPRIME_FORMAT ".\n",
		 poly[1], poly[0], bound, powbound);
  
  p = 2;
  while (p <= bound)
    {
      fbprime_t q;
      /* Handle any prime powers that are smaller than p */
      if (pow_len > 0 && min_pow < p)
	{
	  size_t i;
	  /* Find this p^k in the list of prime powers and replace it 
	     by p^(k+1) if p^(k+1) does not exceed powbound, otherwise
	     remove p^k from the list */
	  for (i = 0; i < pow_len && powers[i] != min_pow; i++);
	  ASSERT_ALWAYS (i < pow_len);
	  q = is_prime_power (min_pow);
	  ASSERT (min_pow / q >= q && min_pow % (q*q) == 0);
	  logp = fb_log (min_pow, log_scale, 0.) - 
	         fb_log (min_pow / q, log_scale, 0.);
	  if (powers[i] <= powbound / q)
	    powers[i] *= q; /* Increase exponent */
	  else
	    {
	      for ( ; i < pow_len - 1; i++) /* Remove from list */
		powers[i] = powers[i + 1];
	      pow_len--;
	    }
	  q = min_pow;
	  /* Find new minimum among the prime powers */
	  if (pow_len > 0)
	    min_pow = powers[0];
	  for (i = 1; i < pow_len; i++)
	    if (powers[i] < min_pow)
	      min_pow =  powers[i];
	}
      else
	{
	  q = p;
	  logp = fb_log (q, log_scale, 0.);
	  /* Do we need to add this prime to the prime powers list? */
	  if (q <= powbound / q)
	    {
	      size_t i;
	      powers = realloc (powers, (++pow_len) * sizeof (fbprime_t));
	      ASSERT_ALWAYS(powers != NULL);
	      powers[pow_len - 1] = q*q;
	      /* Find new minimum among the prime powers */
	      min_pow = powers[0];
	      for (i = 1; i < pow_len; i++)
		if (powers[i] < min_pow)
		  min_pow =  powers[i];
	    }
	  p = getprime (p);
	}
      
      if (!fb_make_linear1 (fb_cur, poly, q, logp, projective))
	continue; /* If root is projective and we don't want those, 
		     skip to next prime */
      
      if (verbose && fb_cur->p > q)
	{
	  if (!had_proj_root)
	    {
	      fprintf (output, "# Primes with projective roots:");
	      had_proj_root = 1;
	    }
	  fprintf (output, " " FBPRIME_FORMAT , q);
	}
      
      fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
      if (fb_new == NULL)
	{
	  free (fb);
	  fb = NULL;
	  break;
	}
      /* fb_fprint_entry (stdout, fb_cur); */
      fb = fb_new;
      
      /* FIXME: handle prime powers */
    }

  getprime (0); /* free prime iterator */

  if (fb != NULL) /* If nothing went wrong so far, put the end-of-fb mark */
    {
      fb_cur->p = FB_END;
      fb_cur->invp = -1L;
      fb_cur->nr_roots = 0;
      fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
      if (fb_new == NULL)
	free (fb);
      fb = fb_new;
    }

  free (fb_cur);
  free (powers);

  if (had_proj_root)
    fprintf (output, "\n");

  rdtscll (tsc2);

  return fb;
}

#if 0
static int fb_comp(const void *a, const void *b) {
  factorbase_degn_t *pa, *pb;
  pa = (factorbase_degn_t *)a;
  pb = (factorbase_degn_t *)b;
  if (pa->p < pb->p) 
    return -1;
  if (pa->p > pb->p) 
    return 1;
  return 0;
}
#endif

#if 0
/* Generate a factor base with prime powers <= bound for a linear polynomial */
static factorbase_degn_t *
fb_make_linear_powers (mpz_t *poly, const fbprime_t bound, const double log_scale, 
		const int verbose)
{
  fbprime_t p;
  factorbase_degn_t *fb = NULL, *fb_cur, *fb_new;
  size_t fbsize = 0, fballoc = 0;
  const size_t allocblocksize = 1 << 20;
  
  rdtscll (tsc1);
  fb_cur = (factorbase_degn_t *) malloc (fb_entrysize_uc (1));
  ASSERT (fb_cur != NULL);

  fb_cur->nr_roots = 1;
  fb_cur->size = fb_entrysize_uc (1);

  if (verbose)
    {
      printf ("# Making factor base for polynomial g(x) = ");
      printf ("\n");
    }

  for (p = 2 ; p <= bound; p = getprime (p))
    {
      modulusul_t m;
      residueul_t r1, r2;
      uint64_t llq;
      fbprime_t q;
      unsigned char old_logpk;

      llq = (uint64_t)p;
      old_logpk = 0;
      while (llq <= (uint64_t)bound) {
	q = (fbprime_t)llq;
	if (q > p)
	  break;

	fb_cur->p = q;
	// TODO: adjust that!
	if (old_logpk == 0) {
	  old_logpk = fb_log (fb_cur->p, log_scale, 0.);
	  fb_cur->plog = old_logpk;
	} else {
	  double newlog;
	  newlog = fb_log_double (fb_cur->p, log_scale, 0.) - (double)old_logpk;
	  old_logpk = (unsigned char) newlog;
	  fb_cur->plog = old_logpk;
	}

	modul_initmod_ul (m, q);
	modul_init_noset0 (r1, m);
	modul_init_noset0 (r2, m);

	/* We want f_1 * a + f_0 * b == 0 (mod p^k) 
	   with a, b coprime and f_1, f_0 coprime

	   p^m = gcd (f_1, p^k)

	   p^k | F(a,b) <=>
	   ( f_0 * b == 0 (mod p^m) <=> b == 0 (mod p^m)  && 
	   (f_1/p^m) * a + f_0 * (b/p^m) == 0 (mod p^(k-m)) )

	   p^m = (p^k, f_1) => 
	   p^k | F(a,b) <=> p^m | b && f(a/b) == 0 (mod p^(k-m))

*/

	modul_set_ul_reduced (r2, mpz_fdiv_ui (poly[1], q), m);
	/* If p | g1 and !(p | g0), p|G(a,b) <=> p|b and that will be handeled
	   by the projective roots */
	if (modul_get_ul (r2, m) == 0)
	  continue;
	modul_set_ul_reduced (r1, mpz_fdiv_ui (poly[0], q), m);

	/* We want g_1 * a + g_0 * b == 0 <=> a/b == - g0 / g1 */
	if (!modul_inv (r2, r2, m)) /* r2 = 1 / g1 */
	  continue;

	modul_mul (r2, r1, r2, m); /* r2 = g0 / g1 */
	modul_neg (r2, r2, m); /* r2 = - g0 / g1 */

#ifndef NDEBUG
	{
	  residueul_t r3;
	  modul_init_noset0 (r3, m);
	  modul_set_ul_reduced (r3, mpz_fdiv_ui (poly[1], q), m);
	  modul_mul (r3, r3, r2, m); /* g1 * (- g0 / g1) */
	  modul_add (r3, r3, r1, m); /* g1 * (- g0 / g1) + g0 */
	  ASSERT (modul_get_ul (r3, m) == 0);
	  modul_clear (r3, m);
	}
#endif

	if (p % 2 != 0)
	  fb_cur->invp = - modul_invmodlong (m);
	fb_cur->roots[0] = modul_get_ul (r2, m);

	modul_clear (r1, m);
	modul_clear (r2, m);
	modul_clearmod (m);

	fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
	if (fb_new == NULL)
	{
	  free (fb);
	  fb = NULL;
	  break;
	}
	/* fb_fprint_entry (stdout, fb_cur); */
	fb = fb_new;

	/* get next power of p */
	llq *= (uint64_t)p;
      }
    }
  getprime (0); /* free prime iterator */

  /* Sort the factor base */
  qsort((void *)fb, fbsize / fb_entrysize_uc(1), fb_entrysize_uc(1), fb_comp);


  if (fb != NULL) /* If nothing went wrong so far, put the end-of-fb mark */
    {
      fb_cur->p = FB_END;
      fb_cur->invp = -1L;
      fb_cur->nr_roots = 0;
      fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
      if (fb_new == NULL)
	free (fb);
      fb = fb_new;
    }

  free (fb_cur);

  rdtscll (tsc2);

  return fb;
}
#endif

factorbase_degn_t *
fb_read_addproj (const char *filename, const double log_scale, 
                 const int verbose, const fbprime_t *proj_primes)
{
  factorbase_degn_t *fb = NULL, *fb_cur, *fb_new;
  FILE *fbfile;
  size_t fbsize = 0, fballoc = 0;
  // too small linesize led to a problem with rsa768;
  // it would probably be a good idea to get rid of fgets
  const size_t linesize = 1000;
  size_t linelen;
  char line[linesize];
  char *lineptr;
  int ok;
  unsigned int i;
  unsigned long linenr = 0;
  const size_t allocblocksize = 1<<20; /* Allocate in MB chunks */
  fbprime_t p, q; /* q is factor base entry q = p^k */
  fbprime_t maxprime = 0;
  unsigned long nr_primes = 0;

  fbfile = fopen (filename, "r");
  if (fbfile == NULL)
    {
      fprintf (stderr, "# Could not open file %s for reading\n", filename);
      return NULL;
    }

  fb_cur = (factorbase_degn_t *) malloc (sizeof (factorbase_degn_t) + 
				      MAXDEGREE * sizeof(fbroot_t));
  if (fb_cur == NULL)
    {
      fprintf (stderr, "# Could not allocate memory for factor base\n");
      fclose (fbfile);
      return NULL;
    }

  while (!feof(fbfile))
    {
      if (fgets (line, linesize, fbfile) == NULL)
	break;
      linenr++;
      linelen = strlen (line);
      if (linelen > 0 && line[linelen - 1] == '\n')
	linelen--; /* Remove newline */
      for (i = 0; i < linelen; i++) /* Skip comments */
        if (line[i] == '#')
          linelen = i;
      while (linelen > 0 && isspace (line[linelen - 1]))
        linelen--; /* Skip whitespace at end of line */
      if (linelen == 0) /* Skip empty/comment lines */
	continue;
      line[linelen] = '\0';

      /* Parse the line */
      lineptr = line;
      q = strtoul (lineptr, &lineptr, 10);
      ok = 0;
      if (q == 0)
	fprintf (stderr, "# fb_read: prime is not an integer on line %lu\n", 
		 linenr);
      else if (*lineptr != ':')
	fprintf (stderr, "# fb_read: prime is not followed by colon on line %lu",
		 linenr);
      else
	ok = 1;
      if (!ok)
	continue;

      if (q <= maxprime)
	{
	  fprintf (stderr, "# Error, primes in factor base file are "
		   "not in increasing order\n");
	  free (fb);
	  fb = NULL;
	  break;
	}

      /* we assume q is a prime or a prime power */
      p = is_prime_power (q);
      if (p != 0)
        {
          fbprime_t cof = q;
          while (cof % p == 0)
            cof /= p;
          if (cof != 1)
            {
              fprintf (stderr, "# Error, " FBPRIME_FORMAT " on line %lu in "
                       "factor base is not a prime or prime power\n",
                       q, linenr);
              break;                       
            }
        }
      else
	p = q; /* If q is prime, p = q */

      /* We have a valid prime or prime power for q. 
         Maybe there's one (or several) prime with only a projective root 
	 smaller than q that we need to add to the factor base first. */
      while (proj_primes != NULL && *proj_primes != 0 && *proj_primes < q)
        {
          fb_cur->p = *proj_primes++;
          fb_cur->plog = fb_log (fb_cur->p, log_scale, 0.);
          fb_cur->nr_roots = 1;
          fb_cur->roots[0] = fb_cur->p;
	  /* Compute invp */
	  if (fb_cur->p % 2 != 0)
	    {
	      modulusul_t m;
	      
	      modul_initmod_ul (m, fb_cur->p);
	      fb_cur->invp = - modul_invmodlong (modul_getmod_ul (m));
	      modul_clearmod (m);
	    }
	  
          fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
          if (fb_new == NULL)
            {
              free (fb);
              fb = NULL;
              break;
            }
          /* fb_fprint_entry (stdout, fb_cur); */
          fb = fb_new;
          maxprime = fb_cur->p;
          nr_primes++;
        }

      fb_cur->p = q;
      if (p == 0) /* It's a prime, do a normal log */
	fb_cur->plog = fb_log (fb_cur->p, log_scale, 0.);
      else
	{
	  /* Take into account the logs of the smaller powers of p that
	     have been added already. We should not just use 
	     fb_cur->plog = log(p) as that would allow rounding errors to 
	     accumulate. */
	  double oldlog;
	  oldlog = fb_log (fb_cur->p / p, log_scale, 0.);
	  fb_cur->plog = fb_log (fb_cur->p, log_scale, - oldlog);
	}

      lineptr++; /* Skip colon */
      ok = 1;
      fb_cur->nr_roots = 0;
      /* Read roots */
      while (ok && *lineptr != '\0')
	{
	  if (fb_cur->nr_roots == MAXDEGREE)
            {
              fprintf (stderr, "# Error, too many roots for " FBPRIME_FORMAT 
                      "in factor base line %lu\n", fb_cur->p, linenr);
              ok = 0;
              break;
            }
	  fbroot_t root = strtoul (lineptr, &lineptr, 10);
	  fb_cur->roots[fb_cur->nr_roots++] = root;
	  /* Check if root is properly reduced. Projective roots s/t (mod p^k)
	     with p|t are stored as p + (t/s mod p^k). */
	  if (root >= 2 * fb_cur->p || (root >= fb_cur->p && root % p != 0))
	    ok = 0;
	  if (*lineptr != '\0' && *lineptr != ',')
	    ok = 0;
	  if (*lineptr == ',')
	    lineptr++;
	}

      if (ok && proj_primes != NULL && fb_cur->p == proj_primes[0])
        {
          int i;
          /* Was the projective root already given in the factor base file? */
          for (i = 0; i < fb_cur->nr_roots; i++)
            if (fb_cur->roots[i] == fb_cur->p)
              break;
          if (i == fb_cur->nr_roots) /* No, need to add it */
            {
              if (fb_cur->nr_roots == MAXDEGREE)
                fprintf (stderr, "# Error, can't add projective root for " 
                         FBPRIME_FORMAT ", already has %d roots\n", 
                        fb_cur->p, fb_cur->nr_roots);
              else
                fb_cur->roots[fb_cur->nr_roots++] = fb_cur->p;
            }
          proj_primes++;
        }

      if (!ok)
	{
	  fprintf (stderr, "# Incorrect format in factor base file line %lu\n",
		   linenr);
	  continue;
	}

      /* Sort the roots into ascending order. We assume that fbprime_t and 
         fbroot_t are typecaste compatible. This is ugly. */
      ASSERT_ALWAYS (sizeof (fbprime_t) == sizeof (fbroot_t));
      fb_sortprimes ((fbprime_t *) fb_cur->roots, fb_cur->nr_roots);
      /* Eliminate duplicate roots. Pari's polrootsmod() can produce them */
      {
        unsigned int i = 0, j;
        for (j = 1; j < (unsigned int) fb_cur->nr_roots; j++)
          {
            if (fb_cur->roots[i] == fb_cur->roots[j])
              {
                fprintf (stderr, "Warning, factor base has repeated root " 
                         FBROOT_FORMAT " for prime " FBPRIME_FORMAT ". "
                         "Keeping only one.\n", fb_cur->roots[i], fb_cur->p);
              }
            else
              fb_cur->roots[++i] = fb_cur->roots[j];
          }
          fb_cur->nr_roots = i + 1;
      }
      
      /* Compute invp */
      if (fb_cur->p % 2 != 0)
	{
	  modulusul_t m;
	  
	  modul_initmod_ul (m, fb_cur->p);
	  fb_cur->invp = - modul_invmodlong (modul_getmod_ul (m));
	  modul_clearmod (m);
	}

      fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
      if (fb_new == NULL)
	{
	  free (fb);
	  fb = NULL;
	  break;
	}
      /* fb_fprint_entry (stdout, fb_cur); */
      fb = fb_new;
      maxprime = fb_cur->p;
      nr_primes++;
    }      

  if (fb != NULL) /* If nothing went wrong so far, put the end-of-fb mark */
    {
      fb_cur->p = FB_END;
      fb_cur->invp = -1L;
      fb_cur->nr_roots = 0;
      fb_new = fb_add_to (fb, &fbsize, &fballoc, allocblocksize, fb_cur);
      if (fb_new == NULL)
	free (fb);
      fb = fb_new;
    }

  if (fb != NULL && verbose)
    {
      printf ("# Factor base sucessfully read, %lu primes, largest was "
	      FBPRIME_FORMAT "\n", nr_primes, maxprime);
    }
  
  fclose (fbfile);
  free (fb_cur);
  return fb;
}

/* return the total size (in bytes) used by fb */
size_t
fb_size (const factorbase_degn_t *fb)
{
  size_t mem = 0;
  while (fb->p != FB_END)
    {
      mem += sizeof (fbprime_t) + sizeof (unsigned long)
        + 4 * sizeof (unsigned char)
        + sizeof (fbroot_t) * fb->nr_roots;
      fb = fb_next (fb);
    }
  return mem;
}

factorbase_degn_t *
fb_read (const char *filename, const double log_scale, const int verbose)
{
  return fb_read_addproj (filename, log_scale, verbose, NULL);
}

/* For all primes p that divide b, disable p and powers of p in fb */
/* Extracts primes (not prime powers) p <= plim with 
   p/nr_roots <= costlim. List ends with FB_END. Allocates memory */
fbprime_t *
fb_extract_bycost (const factorbase_degn_t *fb, const fbprime_t plim, 
                   const fbprime_t costlim)
{
  const factorbase_degn_t *fb_ptr;
  fbprime_t *primes;
  int i;
  
  fb_ptr = fb;
  i = 0;
  while (fb_ptr->p != FB_END && fb_ptr->p <= plim)  
    {
      if (fb_ptr->p <= costlim * fb_ptr->nr_roots &&
	  !is_prime_power(fb_ptr->p))
        i++;
      fb_ptr = fb_next (fb_ptr);
    }

  primes = (fbprime_t *) malloc ((i + 1) * sizeof (fbprime_t));
  ASSERT (primes != NULL);

  fb_ptr = fb;
  i = 0;
  while (fb_ptr->p != FB_END && fb_ptr->p <= plim)  
    {
      if (fb_ptr->p <= costlim * fb_ptr->nr_roots &&
	  !is_prime_power(fb_ptr->p))
        primes[i++] = fb_ptr->p;
      fb_ptr = fb_next (fb_ptr);
    }

  primes[i] = FB_END;
  
  return primes;
}
