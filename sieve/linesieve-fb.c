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
#include "linesieve-fb.h"
#include "utils.h"


static factorbase_degn_t *
fb_find_p (factorbase_degn_t *fb, const fbprime_t p)
{
  while (fb->p != FB_END && fb->p < p) /* Assumes fb is sorted in asc. order */
    {
      ASSERT (fb->p < fb_next (fb)->p || fb_next (fb)->p == FB_END);
      fb = fb_next (fb);
    }

  if (fb->p == p)
    return fb;

  return NULL;
}


/* Sort n primes in array *primes into ascending order */

void
fb_init_firstlog (factorbase_t fb)
{
  int i;
  factorbase_degn_t *fbp;

  for (i = 0; i < 256; i++)
    fb->firstlog[i] = NULL;

  for (fbp = fb->fullfb; fbp->p != FB_END; fbp = fb_next (fbp))
    {
      if (fb->firstlog[fbp->plog] == NULL)
	fb->firstlog[fbp->plog] = fbp;
    }
}


/* For all primes p that divide b, disable p and powers of p in fb */

void
fb_disable_roots (factorbase_degn_t *fb, const unsigned long b, const int verbose)
{
  unsigned long t;
  
  /* Remove the roots of primes that divide b, and of the powers of those 
     primes, from factor base */
  rdtscll (tsc1);
  t = b;
  while (t > 1)
    {
      factorbase_degn_t *fb_del;
      unsigned long p, ppow;
      p = iscomposite (t);
      if (p == 0)
	p = t;
      do t /= p; while (t % p == 0);
      ppow = p;
      while ((fb_del = fb_find_p (fb, ppow)) != NULL)
	{
	  if (verbose)
	    {
	      printf ("# Temporarily removing from factor base ");
	      fb_fprint_entry (stdout, fb_del);
	    }
	  fb_del->nr_roots = 0;
	  ppow *= p;
	}
    }
  rdtscll (tsc2);
}

/* For all primes p that divide b, re-enable p and powers of p in fb */

void
fb_restore_roots (factorbase_degn_t *fb, const unsigned long b, const int verbose)
{
  unsigned long t;

  /* Put the roots back in */
  rdtscll (tsc1);
  t = b;
  while (t > 1)
    {
      factorbase_degn_t *fb_restore;
      unsigned long p, ppow;
      p = iscomposite (t);
      if (p == 0)
	p = t;
      t /= p;
      ppow = p;
      while ((fb_restore = fb_find_p (fb, ppow)) != NULL)
	{
	  fb_restore->nr_roots = 
	    (fb_restore->size - sizeof (factorbase_degn_t)) 
	    / sizeof (fbroot_t);
	  if (verbose)
	    {
	      printf ("# Restored to factor base ");
	      fb_fprint_entry (stdout, fb_restore);
	    }
	  ppow *= p;
	}
    }
  rdtscll (tsc2);
}

#if 0
static fbprime_t
fb_maxprime (factorbase_degn_t *fb)
{
    fbprime_t maxp = 0;
    factorbase_degn_t *fbptr = fb;
    while (fbptr->p != FB_END)
    {
	ASSERT (maxp <= fbptr->p);
	maxp = fbptr->p;
	fbptr = fb_next (fbptr);
    }
    
    return maxp;
}
#endif

/* Extracts primes < bound from fb->fblarge and store them in
   factorbase_small_t format. fb->fblarge gets updated to point
   at smallest prime not extracted. */

void
fb_extract_small (factorbase_t fb, const fbprime_t bound, 
		  const int lvl, const int verbose)
{
  factorbase_degn_t *fbptr;
  unsigned int i, j, size;

  ASSERT_ALWAYS (lvl >= 0 && lvl <= SIEVE_BLOCKING);

  /* Count how many roots there are among the primes <= bound, i.e. how
     many entries the small fb will need */
  size = 0;
  for (fbptr = fb->fblarge; fbptr->p <= bound && fbptr->p != FB_END; 
       fbptr = fb_next (fbptr))
    size += fbptr->nr_roots;

  /* Add one entry for the stop marker */
  size++;

  /* Allocate memory for the small prime factor base and the small prime
     sieve-ready factor base (which may have fewer primes, if some are not 
     coprime to b) */
  fb->fbsmall[lvl] = (factorbase_small_t *) 
                       malloc (size * sizeof (factorbase_small_t));
  fb->fbinit[lvl] = (factorbase_small_inited_t *) 
                        malloc (size * sizeof (factorbase_small_t));
  ASSERT_ALWAYS (fb->fbsmall[lvl] != NULL);
  ASSERT_ALWAYS (fb->fbinit[lvl] != NULL);
  fb->fbsmallsize[lvl] = size;
  fb->fbsmallbound[lvl] = bound;
  
  /* Copy info into the small prime factor base */
  j = 0;
  for (fbptr = fb->fblarge; fbptr->p < bound && fbptr->p != FB_END; 
       fbptr = fb_next (fbptr))
    {
      for (i = 0; i < fbptr->nr_roots; i++, j++)
	{
	  fb->fbsmall[lvl][j].p = fbptr->p;
	  fb->fbsmall[lvl][j].root_and_log = fbptr->roots[i] + 
	    (fbptr->plog << 24);
	}
    }

  /* Put end marker */
  fb->fbsmall[lvl][j].p = FB_END;
  fb->fbsmall[lvl][j].root_and_log = 0;
  j++;

  /* Update the pointer to the large prime factor base */
  
  fb->fblarge = fbptr;

  ASSERT (j == size);

  if (verbose)
    printf ("# There are %u entries in the level %d small factor base\n", 
            j, lvl + 1);
}

/* Test if the data in factorbase is correct. If side == 0, checks the roots
   for the algebraic polynomial; otherwise for the rational one */

int
fb_check (factorbase_t fb, cado_poly poly, int side)
{
  factorbase_degn_t *fbptr;
  fbprime_t p, lastp;
  unsigned long q;
  unsigned char i;

  lastp = 1;
  for (fbptr = fb->fullfb; fbptr->p != FB_END; fbptr = fb_next (fbptr))
    {
      p = fbptr->p;

      if (p <= lastp)
	{
	  fprintf (stderr, "Prime " FBPRIME_FORMAT " in factorbase is not "
		   "greater than previous prime " FBPRIME_FORMAT " \n", 
		   p, lastp);
	  return 1;
	}

      if ((q = iscomposite (p)) != 0)
	{
	  fbprime_t cofac = p;
	  while (cofac > 1 && cofac % (fbprime_t) q == 0)
	    cofac /= (fbprime_t) q;
	  if (cofac != 1)
	    {
	      fprintf (stderr, "Prime " FBPRIME_FORMAT " in factorbase is not "
		       "a prime or prime power\n", p);
	      return 1;
	    }
	}

      if (p % 2 != 0 && fbptr->invp * (unsigned long)p != ~(0UL))
	{
	  fprintf (stderr, "For prime " FBPRIME_FORMAT " in factorbase, "
		   "invp %lu is not -1/p (mod 2^wordsize)\n", p, fbptr->invp);
	  return 1;
	}

      if (fb_entrysize (fbptr) != fbptr->size)
	{
	  fprintf (stderr, "For prime " FBPRIME_FORMAT " in factorbase, size "
		   "entry is %d instead of %d\n", 
		   p, (int) fbptr->size, (int) fb_entrysize (fbptr));
	  return 1;
	}

      for (i = 0; i < fbptr->nr_roots; i++)
	{
	  modulusul_t m;
	  residueul_t val, r, c;
	  unsigned long res;

	  modul_initmod_ul (m, p);
	  modul_init_noset0 (r, m);
	  modul_init_noset0 (c, m);
	  modul_init (val, m);
	  modul_set_ul (r, fbptr->roots[i], m);

	  if (side == 0)
	    {
	      int j;
	      modul_set_ul_reduced (val, mpz_fdiv_ui (poly->f[poly->degree], p), 
				  m);
	      for (j = 1; j <= poly->degree; j++)
		{
		  modul_mul (val, val, r, m);
		  modul_set_ul_reduced 
		    (c, mpz_fdiv_ui (poly->f[poly->degree - j], p), m);
		  modul_add (val, val, c, m);
		}
	    }
	  else
	    {
	      modul_set_ul (val, mpz_fdiv_ui (poly->g[1], p), m);
	      modul_mul (val, val, r, m);
	      modul_set_ul_reduced (c, mpz_fdiv_ui (poly->g[0], p), m);
	      modul_add (val, val, c, m);
	    }

	  res = modul_get_ul (val, m);

	  modul_clear (c, m);
	  modul_clear (r, m);
	  modul_clear (r, m);
	  modul_clearmod (m);

	  if (res != 0)
	    {
	      fprintf (stderr, "False root in factor base: %c(" FBROOT_FORMAT 
		       ") = %lu (mod " FBPRIME_FORMAT ")\n", 
		       (side == 0 ? 'f' : 'g'), fbptr->roots[i], res, p);
	      return 1;
	    }
	}

      /* FIXME: add the tests for the small factor bases */
    }
  return 0;
}

void
fb_clear (factorbase_t fb)
{
  int i;
  free (fb->fullfb);
  fb->fullfb = NULL;
  fb->fblarge = NULL;
  for (i = 0; i < SIEVE_BLOCKING; i++)
    {
      free (fb->fbsmall[i]);
      fb->fbsmall[i] = NULL;
      free (fb->fbinit[i]);
      fb->fbinit[i] = NULL;
    }
}

