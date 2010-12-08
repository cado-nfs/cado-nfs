#include "cado.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for strcmp() */
#include <math.h> /* for sqrt and floor */
#include "utils.h"


typedef struct {
  unsigned long pk;
  unsigned long nroots;
  unsigned long r[10];  // memory is not a problem, here, so take some margin.
} power_ideal_t;


typedef struct {
  power_ideal_t *ideals;
  int len;
  int alloc;
} ideal_heap_t;


/*
 * TODO: here, an insertion sort is implemented. A heap would be much
 * better (quasi-linear vs quadratic).
 * Also: an array of pointer would facilitate swapping elements!
 *
 * Invariant: the ideals are storted in _decreasing_ order, so that
 * removing the minimum costs nothing.
 */


void
include_power_ideal(ideal_heap_t *id_list, const power_ideal_t *ideal) {
  // reallocate if necessary
  if (id_list->alloc == id_list->len) {
    id_list->alloc += 100;
    id_list->ideals = (power_ideal_t *)realloc(id_list->ideals, (id_list->alloc)*sizeof(power_ideal_t));
    ASSERT_ALWAYS(id_list->ideals != NULL);
  }
  // find place
  int i = 0;
  unsigned long pk = ideal->pk;
  while ((i < id_list->len) && (id_list->ideals[i].pk > pk)) 
    i++;
  // if this is the end, then easy.
  if (i == id_list->len) {
    id_list->ideals[i] = *ideal;
    id_list->len++;
    return;
  }
  // otherwise, shift the rest of the array
  power_ideal_t tmp, tmp2;
  tmp = id_list->ideals[i];
  id_list->ideals[i] = *ideal;
  while (i < id_list->len) {
    tmp2 = id_list->ideals[i+1];
    id_list->ideals[i+1] = tmp;
    tmp = tmp2;
    i++;
  }
  id_list->len++;
}

// return 0 if empty heap.
unsigned long 
get_min_power_ideal(const ideal_heap_t *id_list) {
  if (id_list->len == 0)
    return 0;
  return (id_list->ideals[id_list->len-1].pk);
}

void
prune_min_power_ideal(power_ideal_t *ideal, ideal_heap_t *id_list) {
  ASSERT(id_list->len != 0);
  *ideal = id_list->ideals[id_list->len-1];
  id_list->len--;
}


// TODO:
// Stolen from sieve.c. Should be shared, at some point.
void
mp_poly_eval (mpz_t r, mpz_t *poly, int deg, long a)
{
  int i;

  mpz_set (r, poly[deg]);
  for (i = deg - 1; i >= 0; i--)
    {
      mpz_mul_si (r, r, a);
      mpz_add (r, r, poly[i]);
    }
}

// Evaluate the derivative of poly at a.
void
mp_poly_eval_diff (mpz_t r, mpz_t *poly, int deg, long a)
{
  int i;

  mpz_mul_ui (r, poly[deg], (unsigned long)deg);
  for (i = deg - 1; i >= 1; i--)
    {
      mpz_mul_si (r, r, a);
      mpz_addmul_ui (r, poly[i], (unsigned long)i);
    }
}

// return 1 if ok, 0 if failed (multiple root).
//   r <- r - f(r)/f'(r) mod p^k
// Assuming that r was a root mod p^(k-1) this gives a root mod p^k.
int
lift_root(const cado_poly cpoly, const unsigned long pk, unsigned long * r) {
  unsigned long res;
  mpz_t aux, aux2, mp_p;

  mpz_init(aux);
  mpz_init(aux2);
  mpz_init_set_ui(mp_p, pk);
  mp_poly_eval_diff(aux, cpoly->f, cpoly->degree, *r);
  mpz_mod_ui(aux, aux, pk);
  if (!mpz_invert(aux, aux, mp_p)) {
    return 0;
  }
  mp_poly_eval(aux2, cpoly->f, cpoly->degree, *r);
  mpz_mod(aux2, aux2, mp_p);
  mpz_mul(aux2, aux2, aux);
  mpz_neg(aux2, aux2);
  mpz_add_ui(aux2, aux2, *r);
  mpz_mod(aux2, aux2, mp_p);

  res = mpz_get_ui(aux2);
  mpz_clear(aux);
  mpz_clear(aux2);
  mpz_clear(mp_p);

  *r = res;
  return 1;
}

/* This assumes that factor base primes always fit within a
 * long. If they ever grow beyond, then getprime() mmust be changed, and
 * the rootfinding routines must be expanded. */
void
makefb_with_powers(FILE *fp, cado_poly cpoly)
{
  unsigned long p;
  int d = cpoly->degree;
  unsigned long *roots;
  int nroots, i;
  unsigned long power_lim;

  power_lim = (unsigned long)(floor(sqrt((double)(cpoly->alim))));

  check_polynomials (cpoly);

  fprintf (fp, "# Roots for polynomial ");
  fprint_polynomial (fp, cpoly->f, d);

  fprintf (fp, "# DEGREE: %d\n", d);

  roots = (unsigned long*) malloc (d * sizeof (long));

  ideal_heap_t id_list;
  id_list.len = 0;
  id_list.alloc = 0;
  id_list.ideals = NULL;

  for (p = 2; p <= cpoly->alim; p = getprime (p))
    {
      /* Check if there is some power ideal to print before prime p
       * Lemma: there is a prime between a prime and its square.
       */
      unsigned long q = get_min_power_ideal(&id_list);
      while (q != 0 && q < p) {
	power_ideal_t id;
	prune_min_power_ideal(&id, &id_list);
          
	fprintf (fp, "%lu: %lu", id.pk, id.r[0]);
	for (i = 1; i < (int) id.nroots; i++)
	  fprintf (fp, ",%lu", id.r[i]);
	fprintf (fp, "\n");

	q = get_min_power_ideal(&id_list);
      }
      nroots = poly_roots_ulong(roots, cpoly->f, d, p);

      if (nroots != 0)
        {
          fprintf (fp, "%lu: %lu", p, roots[0]);
          for (i = 1; i < nroots; i++)
            fprintf (fp, ",%lu", roots[i]);
          fprintf (fp, "\n");
        }

      /* Create prime powers, if relevant */
      if (nroots != 0 && p <= power_lim) {
	uint64_t pk = (uint64_t)(p*p);
	power_ideal_t id;
	do {
	  id.pk = pk;
	  id.nroots = 0;
	  for (i = 0; i < nroots; ++i) {
	    if (lift_root(cpoly, pk, &(roots[i]))) {
	      id.r[id.nroots] = roots[i];
	      id.nroots++;
	    }
	  }
	  if (id.nroots == 0)
	    break;
	  include_power_ideal(&id_list, &id);
  	  pk *= p;
	} while (pk <= cpoly->alim);
      }
    }

  getprime (0); /* free the memory used by getprime() */
  free (roots);
}





void
makefb (FILE *fp, cado_poly cpoly)
{
  unsigned long p;
  int d = cpoly->degree;
  unsigned long *roots;
  int nroots, i;

  check_polynomials (cpoly);

  fprintf (fp, "# Roots for polynomial ");
  fprint_polynomial (fp, cpoly->f, d);

  fprintf (fp, "# DEGREE: %d\n", d);

  roots = (unsigned long*) malloc (d * sizeof (unsigned long));

  for (p = 2; p <= cpoly->alim; p = getprime (p))
    {
        nroots = poly_roots_ulong(roots, cpoly->f, d, p);
        // TODO: poly_roots_ulong returns 0 if f mod p is 0.
        // This corresponds to complicated roots that should maybe go to
        // the factor base (hum... maybe not!)
      if (nroots != 0)
        {
          fprintf (fp, "%lu: %lld", p, (long long int) roots[0]);
          for (i = 1; i < nroots; i++)
            fprintf (fp, ",%lld", (long long int) roots[i]);
          fprintf (fp, "\n");
        }
    }

  getprime (0); /* free the memory used by getprime() */
  free (roots);
}

static void usage()
{
    fprintf (stderr, "Usage: makefb [-v] [-powers] [-poly file]\n");
    exit (1);
}

int
main (int argc, char *argv[])
{
  int verbose = 0;
  int use_powers = 0;
  param_list pl;
  cado_poly cpoly;
  FILE * f;

  param_list_init(pl);
  cado_poly_init(cpoly);

  param_list_configure_knob(pl, "-v", &verbose);
  param_list_configure_knob(pl, "-powers", &use_powers);

  argv++, argc--;
  for( ; argc ; ) {
      if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }

      /* Could also be a file */
      if ((f = fopen(argv[0], "r")) != NULL) {
          param_list_read_stream(pl, f);
          fclose(f);
          argv++,argc--;
          continue;
      }

      fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
      usage();
  }

  const char * filename;
  if ((filename = param_list_lookup_string(pl, "poly")) != NULL) {
      param_list_read_file(pl, filename);
  }

  if (!cado_poly_set_plist (cpoly, pl))
    {
      fprintf (stderr, "Error reading polynomial file\n");
      exit (EXIT_FAILURE);
    }
  param_list_clear(pl);

  if (use_powers)
    makefb_with_powers (stdout, cpoly);
  else
    makefb (stdout, cpoly);

  cado_poly_clear (cpoly);

  return 0;
}
