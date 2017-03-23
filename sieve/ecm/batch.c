#include "cado.h"
#include <stdio.h>
#include <math.h>
#ifdef  HAVE_OPENMP
#include <omp.h>
#endif
#include "batch.h"
#include "utils.h"


static void
ulong_list_init (ulong_list L)
{
  L->l = NULL;
  L->alloc = L->size = 0;
}

/* at any point, all values up to L->alloc should be mpz_init()'ed,
   even beyond L->size */
static void
ulong_list_add (ulong_list L, unsigned long n)
{
  if (L->size == L->alloc)
    {
      L->alloc = 3 * (L->alloc / 2) + 2;
      L->l = realloc (L->l, L->alloc * sizeof (mpz_t));
    }
  ASSERT(L->size < L->alloc);
  L->l[L->size] = n;
  L->size ++;
}

static void
ulong_list_clear (ulong_list L)
{
  free (L->l);
  L->alloc = L->size = 0;
}

static void
mpz_product_tree_init (mpz_product_tree t)
{
  t->l = NULL;
  t->n = NULL;
  t->size = 0;
}

/* add a new entry n to product tree t */
static void
mpz_product_tree_add_ui (mpz_product_tree t, unsigned long n)
{
  if (t->size == 0) /* tree was empty */
    {
      t->l = malloc (sizeof (mpz_t));
      mpz_init_set_ui (t->l[0], n);
      t->n = malloc (sizeof (unsigned long));
      t->n[0] = 1;
      t->size = 1;
    }
  else
    {
      /* first accumulate in l[0] */
      if (t->n[0] == 0)
        mpz_set_ui (t->l[0], n);
      else
        mpz_mul_ui (t->l[0], t->l[0], n);
      t->n[0] ++;
      for (unsigned int i = 0; t->n[i] == (2UL<<i); i++)
        {
          if (i+1 == t->size) /* realloc */
            {
              t->l = realloc (t->l, (t->size + 1) * sizeof (mpz_t));
              mpz_init_set_ui (t->l[t->size], 1);
              t->n = realloc (t->n, (t->size + 1) * sizeof (unsigned long));
              t->n[t->size] = 0;
              t->size++;
            }
          if (t->n[i+1] == 0)
            {
              mpz_swap (t->l[i+1], t->l[i]);
              mpz_set_ui (t->l[i], 1);
            }
          else /* accumulate */
            mpz_mul (t->l[i+1], t->l[i+1], t->l[i]);
          t->n[i+1] += t->n[i];
          t->n[i] = 0;
        }
    }
}

/* accumulate all products in t->l[0] */
static void
mpz_product_tree_accumulate (mpz_product_tree t)
{
  for (unsigned int i = 1; i < t->size; i++)
    {
      mpz_mul (t->l[0], t->l[0], t->l[i]);
      t->n[0] += t->n[i];
      t->n[i] = 0;
    }
}

static void
mpz_product_tree_clear (mpz_product_tree t)
{
  for (unsigned int i = 0; i < t->size; i++)
    mpz_clear (t->l[i]);
  free (t->l);
  free (t->n);
  t->size = 0;
}

void
cofac_list_init (cofac_list l)
{
  l->a = NULL;
  l->b = NULL;
  l->R = NULL;
  l->A = NULL;
  l->R0 = NULL;
  l->A0 = NULL;
  l->sq = NULL;
  l->side = NULL;
  l->perm = NULL;
  l->alloc = 0;
  l->size = 0;
}

static unsigned long
prime_list (ulong_list L, prime_info pi, unsigned long pmin,
            unsigned long pmax)
{
  unsigned long p;

  for (p = pmin; p < pmax; p = getprime_mt (pi))
    ulong_list_add (L, p);
  return p;
}

static unsigned long
prime_list_poly (ulong_list L, prime_info pi, unsigned long pmin,
                 unsigned long pmax, mpz_poly f)
{
  unsigned long p;

  if (f->deg == 1)
    return prime_list (L, pi, pmin, pmax);

  for (p = pmin; p < pmax; p = getprime_mt (pi))
    if (mpz_divisible_ui_p (f->coeff[f->deg], p) ||
        mpz_poly_roots_ulong (NULL, f, p) > 0)
      ulong_list_add (L, p);
  return p;
}

/* add in the product tree L all primes pmin <= p < pmax.
   Assume pmin is the current prime in 'pi'
   (pmin=2 when 'pi' was just initialized).
   Return the current prime in 'pi' at the end, i.e.,
   the smallest prime >= pmax. */
static unsigned long
prime_tree (mpz_product_tree L, prime_info pi, unsigned long pmin,
            unsigned long pmax)
{
  unsigned long p;

  for (p = pmin; p < pmax; p = getprime_mt (pi))
    mpz_product_tree_add_ui (L, p);
  return p;
}

/* same as prime_tree, but only adds primes p for which f has at least one
   root modulo p, or the leading coefficient of f vanishes modulo p */
static unsigned long
prime_tree_poly (mpz_product_tree L, prime_info pi, unsigned long pmin,
                 unsigned long pmax, mpz_poly f)
{
  unsigned long p, *q;
  int i, j, nthreads = 1;

  if (f->deg == 1)
    return prime_tree (L, pi, pmin, pmax);

#ifdef HAVE_OPENMP
#pragma omp parallel
  nthreads = omp_get_num_threads ();
#endif
  nthreads *= 10; /* to amortize the varying cost of mpz_poly_roots_ulong */
  q = malloc (nthreads * sizeof (unsigned long));

  for (p = pmin; p < pmax;)
    {
      /* sequential part: getprime_mt is fast */
      for (i = 0; i < nthreads && p < pmax; p = getprime_mt (pi), i++)
        q[i] = p;

      /* parallel part: mpz_poly_roots_ulong is the bottleneck */
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (j = 0; j < i; j++)
        if (mpz_divisible_ui_p (f->coeff[f->deg], q[j]) == 0 &&
            mpz_poly_roots_ulong (NULL, f, q[j]) == 0)
          q[j] = 0;

      /* sequential part: mpz_product_tree_add_ui is fast */
      for (j = 0; j < i; j++)
        if (q[j])
          mpz_product_tree_add_ui (L, q[j]);
    }

  free (q);

  return p;
}

void
cofac_list_realloc (cofac_list l, size_t newsize)
{
  unsigned long i;

  /* if we shrink the list, clear the mpz_t's */
  for (i = newsize; i < l->size; i++)
    {
      mpz_clear (l->R[i]);
      mpz_clear (l->A[i]);
      mpz_clear (l->R0[i]);
      mpz_clear (l->A0[i]);
      mpz_clear (l->sq[i]);
    }
  l->a = realloc (l->a, newsize * sizeof (int64_t));
  l->b = realloc (l->b, newsize * sizeof (uint64_t));
  l->R = realloc (l->R, newsize * sizeof (mpz_t));
  l->A = realloc (l->A, newsize * sizeof (mpz_t));
  l->R0 = realloc (l->R0, newsize * sizeof (mpz_t));
  l->A0 = realloc (l->A0, newsize * sizeof (mpz_t));
  l->sq = realloc (l->sq, newsize * sizeof (mpz_t));
  l->side = realloc (l->side, newsize * sizeof (int));
  l->perm = realloc (l->perm, newsize * sizeof (uint32_t));
  l->alloc = newsize;
  if (newsize < l->size)
    l->size = newsize;
}

void
cofac_list_add (cofac_list l, long a, unsigned long b, mpz_t R, mpz_t A,
                int side, mpz_t sq)
{
  if (l->size == l->alloc)
    cofac_list_realloc (l, 2 * l->alloc + 1);
  l->a[l->size] = a;
  l->b[l->size] = b;
  mpz_init_set (l->R[l->size], R);
  mpz_init_set (l->A[l->size], A);
  mpz_init_set (l->R0[l->size], R);
  mpz_init_set (l->A0[l->size], A);
  mpz_init_set (l->sq[l->size], sq);
  l->side[l->size] = side;
  l->perm[l->size] = l->size;
  (l->size)++;
}

void
cofac_list_clear (cofac_list l)
{
  unsigned long i;
  for (i = 0; i < l->size; i++)
    {
      mpz_clear (l->R[i]);
      mpz_clear (l->A[i]);
      mpz_clear (l->R0[i]);
      mpz_clear (l->A0[i]);
      mpz_clear (l->sq[i]);
    }
  free (l->a);
  free (l->b);
  free (l->R);
  free (l->A);
  free (l->R0);
  free (l->A0);
  free (l->sq);
  free (l->side);
  free (l->perm);
}

/* put in P the product of primes p for which the given polynomial has factors
   modulo p */
static void
prime_product_poly (mpz_t P, prime_info pi, unsigned long p_max,
                    unsigned long p_last, mpz_poly f)
{
  mpz_product_tree L;

  mpz_product_tree_init (L);
  p_last = prime_tree_poly (L, pi, p_last, p_max, f);

  mpz_product_tree_accumulate (L);

  mpz_set (P, L->l[0]);

  mpz_product_tree_clear (L);
}

static unsigned long
tree_height (unsigned long n)
{
  unsigned long h = 0;

  while (n > 1)
    {
      h ++;
      n = (n + 1) / 2;
    }
  return h;
}

/* return the product tree formed from R[0..n-1].
   Put in w[i] the number of elements of level i:
   w[0] = n, w[1] = ceil(n/2), ... */
static mpz_t**
product_tree (mpz_t *R, uint32_t *perm, unsigned long n, unsigned long *w)
{
  unsigned long h = tree_height (n), i, j;
  mpz_t **T;

  ASSERT_ALWAYS(n >= 1);

  T = (mpz_t**) malloc ((h + 1) * sizeof (mpz_t*));

  /* initialize tree */
  w[0] = n;
  for (i = 0; i <= h; i++)
    {
      w[i] = 1 + ((n - 1) >> i);
      T[i] = (mpz_t*) malloc (w[i] * sizeof (mpz_t));
      for (j = 0; j < w[i]; j++)
        mpz_init (T[i][j]);
    }

  /* initialize T[0] to R */
  for (j = 0; j < n; j++)
    mpz_set (T[0][j], R[perm[j]]);

  /* compute product tree */
  for (i = 1; i <= h; i++)
    {
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (j = 0; j < w[i-1] / 2; j++)
        mpz_mul (T[i][j], T[i-1][2*j], T[i-1][2*j+1]);
      if (w[i-1] & 1)
        mpz_set (T[i][w[i]-1], T[i-1][w[i-1]-1]);
    }

  return T;
}

static void
remainder_tree_aux (mpz_t **T, unsigned long **nbits, unsigned long i,
                    unsigned long j, unsigned long guard)
{
  /* T[i][j]/2^(nbits[i][j] + guard) ~ P/T[i][j] */
  mpz_mul (T[i-1][2*j], T[i][j], T[i-1][2*j]);
  /* same for the right part */
  mpz_mul (T[i-1][2*j+1], T[i][j], T[i-1][2*j+1]);
  /* swap */
  mpz_swap (T[i-1][2*j], T[i-1][2*j+1]);

  /* get the fractional part, i.e., the low nbits[i][j] + guard bits */
  mpz_tdiv_r_2exp (T[i-1][2*j], T[i-1][2*j], nbits[i][j] + guard);
  /* now keep only nbits[i-1][2*j] + guard significant bits */
  mpz_div_2exp (T[i-1][2*j], T[i-1][2*j], nbits[i][j] - nbits[i-1][2*j]);

  mpz_tdiv_r_2exp (T[i-1][2*j+1], T[i-1][2*j+1], nbits[i][j] + guard);
  mpz_div_2exp (T[i-1][2*j+1], T[i-1][2*j+1], nbits[i][j] - nbits[i-1][2*j+1]);
}

/* Compute the remainder tree using the "scaled" variant
   (http://cr.yp.to/arith/scaledmod-20040820.pdf).
   At the root, we compute a floating-point
   approximation of P/T[h][0] with m+guard bits, where m = nbits(T[h][0]). */
static void
remainder_tree (mpz_t **T, unsigned long n, unsigned long *w, mpz_t P,
                mpz_t *R, uint32_t *perm)
{
  unsigned long h = tree_height (n), i, j, guard;
  unsigned long **nbits;
  mpz_t Q;

  guard = h;
  nbits = malloc ((h + 1) * sizeof (unsigned long*));
  for (i = 0; i <= h; i++)
    {
      nbits[i] = malloc (w[i] * sizeof (unsigned long));
      for (j = 0; j < w[i]; j++)
        nbits[i][j] = mpz_sizeinbase (T[i][j], 2);
    }

  mpz_init (Q);
  mpz_mod (Q, P, T[h][0]); /* first reduce modulo T[h][0] in case P is huge */
  mpz_mul_2exp (Q, Q, nbits[h][0] + guard);
  mpz_tdiv_q (T[h][0], Q, T[h][0]);
  /* P/T[h][0] ~ Q/2^(m+guard) */
  for (i = h; i > 0; i--)
    {
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (j = 0; j < w[i-1] / 2; j++)
        remainder_tree_aux (T, nbits, i, j, guard);
      if (w[i-1] & 1)
        mpz_swap (T[i-1][w[i-1]-1], T[i][w[i]-1]);
    }

  /* from T[0][j] ~ P/R[j]*2^(nbits[0][j] + guard) mod 2^(nbits[0][j] + guard),
     get T[0][j]*R[j]/2^(nbits[0][j] + guard) ~ P mod R[j] */
  for (j = 0; j < n; j++)
    {
      mpz_mul (T[0][j], T[0][j], R[perm[j]]);
      /* T[0][j] ~ P*2^(nbits[0][j] + guard) mod R[j]*2^(nbits[0][j]+guard) */
      mpz_div_2exp (T[0][j], T[0][j], nbits[0][j]);
      /* T[0][j] ~ P*2^guard mod R[j]*2^guard */
      mpz_add_ui (T[0][j], T[0][j], (1UL << guard) - 1UL);
      mpz_div_2exp (T[0][j], T[0][j], guard);
    }

  mpz_clear (Q);
  for (i = 0; i <= h; i++)
    free (nbits[i]);
  free (nbits);
}

/* Clear the product tree. */
static void
clear_product_tree (mpz_t **T, unsigned long n, unsigned long *w)
{
  unsigned long i, j, h = tree_height (n);

  for (i = 0; i <= h; i++)
    {
      for (j = 0; j < w[i]; j++)
        mpz_clear (T[i][j]);
      free (T[i]);
    }
  free (T);
}

#define MAX_DEPTH 32

/* Input:
   R[0], ..., R[n-1] are cofactors
   P is the product of primes
   Output:
   Each R[j] has been divided by its P-smooth part.
*/
static void
smoothness_test (mpz_t *R, uint32_t *perm, unsigned long n, mpz_t P, FILE *out)
{
  unsigned long j, w[MAX_DEPTH];
  mpz_t **T;
  double st, wct;

  if (n == 0)
    return;

  st = seconds ();
  wct = wct_seconds ();
  T = product_tree (R, perm, n, w);
  unsigned long h = tree_height (n);
  fprintf (out, "# batch: took %.2fs (wct %.2fs) to compute product tree of %zu bits\n",
           seconds () - st, wct_seconds () - wct, mpz_sizeinbase (T[h][0], 2));

  /* compute remainder tree */
  st = seconds ();
  wct = wct_seconds ();
  remainder_tree (T, n, w, P, R, perm);
  fprintf (out, "# batch: took %.2fs (wct %.2fs) to compute remainder tree\n",
           seconds () - st, wct_seconds () - wct);

  /* now T[0][j] = P mod R[j] for 0 <= j < n */
  for (j = 0; j < n; j++)
    {
      mpz_gcd (T[0][j], T[0][j], R[perm[j]]);
      mpz_divexact (R[perm[j]], R[perm[j]], T[0][j]);
    }

  clear_product_tree (T, n, w);
}

/* invariant:
   relations 0 to *nb_smooth-1 are smooth
   relations *nb_smooth to *nb_unknown-1 are unknown
   relations >= *nb_unknown are non-smooth (useless).

   on input, we know that n has no factor<=B. 
   We're willing to accept prime cofactors <= L,
   and to try harder to factor composites which are <= M
*/
static void
update_status (mpz_t *R, uint32_t *perm,
               unsigned char *b_status_r, unsigned char *b_status_a,
               mpz_srcptr B,
               mpz_srcptr L,
               mpz_srcptr M,
               unsigned long *nb_smooth, unsigned long *nb_unknown)
{
  unsigned long i, j;

  mpz_t BB;
  mpz_init(BB);
  mpz_mul(BB,B,B);

  for (j = *nb_smooth; j < *nb_unknown; j++)
    {
      i = perm[j];
      if (b_status_r[i] == STATUS_UNKNOWN)
      {
        /* relation i is smooth iff R[i]=1 ; another option is in case
         * the remaining cofactor is below the mfb we've been given. */
        if (mpz_cmp_ui (R[i], 1) == 0
                || mpz_cmp(R[i], L) <= 0
                || (
                    mpz_cmp(R[i], BB) >= 0
                    && mpz_cmp(R[i], M) <= 0
                    && !mpz_probab_prime_p(R[i], 1)
                ))
          {
            b_status_r[i] = STATUS_SMOOTH;
            if (b_status_a[i] == STATUS_SMOOTH)
              {
                /* relation j is smooth, swap it with relation *nb_smooth */
                perm[j] = perm[*nb_smooth];
                perm[*nb_smooth] = i;
                (*nb_smooth)++;
              }
          }
        else /* not smooth */
          {
            /* relation j is useless, swap it with relation *nb_unknown - 1 */
            (*nb_unknown)--;
            perm[j] = perm[*nb_unknown];
            perm[*nb_unknown] = i;
            j--;
          }
      }
    }
  mpz_clear(BB);
}

/* return the number n of smooth relations in l,
   which should be at the end in locations perm[0], perm[1], ..., perm[n-1] */
unsigned long
find_smooth (cofac_list l, mpz_t batchP[2], mpz_t B[2], mpz_t L[2], mpz_t M[2], FILE *out,
             int nthreads MAYBE_UNUSED)
{
  unsigned long nb_rel_read = l->size;
  unsigned long nb_smooth;
  unsigned long nb_unknown;
  unsigned char *b_status_r;
  unsigned char *b_status_a;
  double start = seconds ();

#ifdef HAVE_OPENMP
  omp_set_num_threads (nthreads);
#endif

  b_status_r = (unsigned char *) malloc (nb_rel_read * sizeof(unsigned char));
  b_status_a = (unsigned char *) malloc (nb_rel_read * sizeof(unsigned char));
  memset (b_status_r, STATUS_UNKNOWN, nb_rel_read);
  memset (b_status_a, STATUS_UNKNOWN, nb_rel_read);

  nb_smooth = 0;
  nb_unknown = nb_rel_read;

  /* invariant: the smooth relations are in 0..nb_smooth-1,
     the unknown ones in nb_smooth..nb_unknown-1,
     the remaining ones are not smooth */

  /* it seems faster to start from the algebraic side */
  for (int z = 1; z >= 0; z--)
    {
      if (z == 0)
        smoothness_test (l->R, l->perm + nb_smooth, nb_unknown - nb_smooth,
                         batchP[0], out);
      else
        smoothness_test (l->A, l->perm + nb_smooth, nb_unknown - nb_smooth,
                         batchP[1],  out);

      /* we only need to update relations in [nb_smooth, nb_unknown-1] */
      if (z == 0)
        update_status (l->R, l->perm, b_status_r, b_status_a, B[0], L[0], M[0],
                       &nb_smooth, &nb_unknown);
      else
        update_status (l->A, l->perm, b_status_a, b_status_r, B[1], L[1], M[1],
                       &nb_smooth, &nb_unknown);
    }

  free (b_status_r);
  free (b_status_a);

  fprintf (out, "# batch: took %.2fs to detect %lu smooth relations out of %lu\n", seconds () - start, nb_smooth, nb_rel_read);

  return nb_smooth;
}

/* print all factors of n into 'str', where 'str0' is the string start.
   B is the small prime bound: n contain no prime factor <= B
   lpb is the large prime bound: all prime factor of n are < 2^lpb */
static char*
print_smooth_aux (mpz_t *factors, mpz_t n, facul_method_t *methods,
                  struct modset_t *fm, struct modset_t *cfm,
                  int lpb, double B, char *str0, char *str, size_t str0size)
{
  unsigned long i;
  int j, res_fac;
  double BB = B * B, BBB = B * B * B;

  /* any factor < B^2 is necessarily prime */
  for (i = 0; methods[i].method != 0 && mpz_cmp_d (n, BB) >= 0; i++)
    {
      res_fac = facul_doit_onefm_mpz (factors, n, methods[i], fm, cfm,
                                      lpb, BB, BBB);

      /* Could happen if we allowed a cofactor bound after batch
       * cofactorization */
      if (res_fac == FACUL_NOT_SMOOTH) {
          gmp_fprintf(stderr, "# C=%Zd not smooth\n", n);
          return str;
      }


      ASSERT_ALWAYS(res_fac != FACUL_NOT_SMOOTH);

      /* factors[0..res_fac-1] are prime factors of n, 0 <= res_fac <= 2 */

      for (j = 0; j < res_fac; j++)
        {
          if (str > str0)
            str += snprintf (str, str0size - (str - str0), ",");
          str += gmp_snprintf (str, str0size - (str - str0), "%Zx", factors[j]);
          mpz_divexact (n, n, factors[j]);
        }

      if (fm->arith != CHOOSE_NONE) /* fm is a non-trivial composite factor */
        {
          mpz_t t;
          struct modset_t cfm2;
          mpz_init (t);
          modset_get_z (t, fm);
          /* t should be composite, i.e., t >= BB */
          ASSERT(mpz_cmp_d (t, BB) >= 0);
          mpz_divexact (n, n, t);
          str = print_smooth_aux (factors, t, methods + i + 1, fm, &cfm2,
                                  lpb, B, str0, str, str0size);
          modset_clear (fm);
          fm->arith = CHOOSE_NONE;
          mpz_clear (t);
        }

      if (cfm->arith != CHOOSE_NONE)
        {
          mpz_t t;
          struct modset_t fm2;
          mpz_init (t);
          modset_get_z (t, cfm);
          /* t should be composite, i.e., t >= BB */
          ASSERT(mpz_cmp_d (t, BB) >= 0);
          mpz_divexact (n, n, t);
          str = print_smooth_aux (factors, t, methods + i + 1, &fm2, cfm,
                                  lpb, B, str0, str, str0size);
          modset_clear (cfm);
          cfm->arith = CHOOSE_NONE;
          mpz_clear (t);
        }
    }

  if (mpz_cmp_ui (n, 1) > 0)
    {
      ASSERT (mpz_cmp_d (n, BB) < 0);
      if (str > str0)
        str += snprintf (str, str0size - (str - str0), ",");
      str += gmp_snprintf (str, str0size - (str - str0), "%Zx", n);
    }

  return str;
}

/* return the end of the written string */
static char*
trial_divide (mpz_t n, unsigned long *sp, unsigned long spsize, char *str, size_t strsize)
{
  unsigned long i;
  char *str0 = str;

  for (i = 0; i < spsize; i++)
    {
      while (mpz_divisible_ui_p (n, sp[i]))
        {
          if (str > str0)
            str += snprintf (str, strsize - (str - str0), ",");
          str += snprintf (str, strsize - (str - str0), "%lx", sp[i]);
          mpz_divexact_ui (n, n, sp[i]);
        }
    }

  return str;
}

/* Print the prime factors of the input 'n' separated by spaces.
   The list sp[] contains small primes (less than B).
   B is the small prime bound: any factor < B^2 is necessarily prime.
   'cofac' is the initial cofactor (without special-q).
   factors[2] is a scratch space to store factors.
   Return the pointer to the next character to write.
*/
static char*
print_smooth (mpz_t *factors, mpz_t n, facul_method_t *methods,
              struct modset_t *fm, struct modset_t *cfm,
              unsigned int lpb, double B, unsigned long *sp,
              unsigned long spsize, mpz_t cofac, mpz_ptr sq,
              char *str, size_t strsize)
{
  char *str0 = str;

  ASSERT(mpz_divisible_p (n, cofac));
  mpz_divexact (n, n, cofac);

  if (sq != NULL)
    {
      ASSERT(mpz_divisible_p (n, sq));
      mpz_divexact (n, n, sq);
    }

  /* remove small primes */
  str = trial_divide (n, sp, spsize, str, strsize - (str - str0));

  /* factor the cofactor */
  str = print_smooth_aux (factors, cofac, methods, fm, cfm, lpb, B, str0, str, strsize - (str - str0));

  /* factor rest of factor base primes */
  str = print_smooth_aux (factors, n, methods, fm, cfm, lpb, B, str0, str, strsize - (str - str0));

  if (sq != NULL)
    {
      if (str > str0)
        str += snprintf (str, strsize - (str - str0), ",");
      str += gmp_snprintf (str, strsize - (str - str0), "%Zx", sq);
    }

  return str;
}

/* strip integers in l[0..n-1] which do not divide P */
unsigned long
strip (unsigned long *l, unsigned long n, mpz_t P)
{
  unsigned long i, j;

  for (i = j = 0; i < n; i++)
    if (mpz_divisible_ui_p (P, l[i]))
      l[j++] = l[i];
  return j;
}

/* sqside = 1 if the special-q is on side 1 (algebraic) */
static void
factor_one (cofac_list L, cado_poly pol, unsigned long *lim, int *lpb,
            FILE *out, facul_method_t *methods, unsigned long *sp[],
            unsigned long spsize[], unsigned long i)
{
  mpz_t norm;
  mpz_t factors[2];
  uint32_t *perm = L->perm;
  struct modset_t fm, cfm;
  char s0[1024];        /* output relation */
  char *s = s0;

  mpz_init (norm);
  mpz_init (factors[0]);
  mpz_init (factors[1]);

  /* compute norms F(a,b) and G(a,b) */

  mpz_poly_homogeneous_eval_sisi (norm, pol->pols[0],
                                  L->a[perm[i]], -L->b[perm[i]]);

  int bn = L->b[perm[i]] < 0;

  s += snprintf (s, sizeof(s0)-(s-s0), "X %" PRIx64 ",%s%" PRIx64 ":",
          L->a[perm[i]],
          bn ? "-" : "",
          bn ? -L->b[perm[i]] : L->b[perm[i]]);

  s = print_smooth (factors, norm, methods, &fm, &cfm, lpb[0], (double) lim[0],
                    sp[0], spsize[0], L->R0[perm[i]],
                    (L->side[perm[i]] == 0) ? L->sq[perm[i]] : NULL, s, sizeof(s0)-(s-s0));
  s += snprintf (s, sizeof(s0)-(s-s0), ":");

  mpz_poly_homogeneous_eval_sisi (norm, pol->pols[1],
                                  L->a[perm[i]], -L->b[perm[i]]);
  s = print_smooth (factors, norm, methods, &fm, &cfm, lpb[1], (double) lim[1],
                    sp[1], spsize[1], L->A0[perm[i]],
                    (L->side[perm[i]] == 1) ? L->sq[perm[i]] : NULL, s, sizeof(s0)-(s-s0));

  /* avoid two threads writing a relation simultaneously */
#ifdef  HAVE_OPENMP
#pragma omp critical
#endif
  {
    fprintf (out, "%s\n", s0);
    fflush (out);
  }

  mpz_clear (norm);
  mpz_clear (factors[0]);
  mpz_clear (factors[1]);
}

/* Given a list L of bi-smooth cofactors, print the corresponding relations
   on "out".
   n is the number of bi-smooth cofactors in L.
*/
void
factor (cofac_list L, unsigned long n, cado_poly pol, int lpb[],
        FILE *out, int nthreads MAYBE_UNUSED)
{
  unsigned long i, *sp[2], spsize[2], B[2];
  int nb_methods;
  facul_method_t *methods;
  ulong_list SP0, SP1;
  prime_info pi;
  double start = seconds (), wct_start = wct_seconds ();

  ulong_list_init (SP0);
  prime_info_init (pi);
  B[0] = (unsigned long) ceil (pow (2.0, (double) lpb[0] / 2.0));
  prime_list_poly (SP0, pi, 2, B[0], pol->pols[0]);
  spsize[0] = SP0->size;
  sp[0] = SP0->l;
  prime_info_clear (pi);

  ulong_list_init (SP1);
  prime_info_init (pi);
  B[1] = (unsigned long) ceil (pow (2.0, (double) lpb[1] / 2.0));
  prime_list_poly (SP1, pi, 2, B[1], pol->pols[1]);
  spsize[1] = SP1->size;
  sp[1] = SP1->l;
  prime_info_clear (pi);

  nb_methods = 30;
  if (nb_methods >= NB_MAX_METHODS)
    nb_methods = NB_MAX_METHODS - 1;
  methods = facul_make_default_strategy (nb_methods - 3, 0);

#ifdef HAVE_OPENMP
  omp_set_num_threads (nthreads);
#pragma omp parallel for schedule(static)
#endif
  for (i = 0; i < n; i++)
    factor_one (L, pol, B, lpb, out, methods, sp, spsize, i);

  ulong_list_clear (SP0);
  ulong_list_clear (SP1);

  facul_clear_methods (methods);

  fprintf (out, "# batch: took %.2fs (wct %.2fs) to factor %lu smooth relations\n",
           seconds () - start, wct_seconds () - wct_start, n);
}

static void
create_batch_product (mpz_t P, unsigned long B, unsigned long L, mpz_poly pol)
{
  prime_info pi;
  unsigned long p;

  ASSERT_ALWAYS (L > B);

  prime_info_init (pi);

  for (p = 2; p < B; p = getprime_mt (pi));

  prime_product_poly (P, pi, L, p, pol);

  prime_info_clear (pi);
}

/* We have 3 cases:
   1) if f == NULL: P is computed but not stored
   2) if f != NULL but file is non-existing: P is computed and saved in f
   3) if f != NULL and file is existing: P is read from file
*/
void
create_batch_file (const char *f, mpz_t P, unsigned long B, unsigned long L,
                   mpz_poly pol, FILE *out, int nthreads MAYBE_UNUSED)
{
  FILE *fp;
  double s = seconds (), wct = wct_seconds ();
  size_t ret;

  // the product of primes up to B takes \log2(B)-\log\log 2 / \log 2
  // bits. The added constant is 0.5287.
  if (log2(B/GMP_LIMB_BITS) + 0.5287 >= 31) {
      fprintf(stderr, "Gnu MP cannot deal with primes product that large (max 37 bits)\n");
      abort();
  } else if (log2(B) + 0.5287 >= 34) {
      fprintf(stderr, "Gnu MP's mpz_inp_raw and mpz_out_raw functions are limited to integers of at most 34 bits\n");
      abort();
  }

#ifdef HAVE_OPENMP
  omp_set_num_threads (nthreads);
#endif

  if (f == NULL) /* case 1 */
    {
      fprintf (out, "# batch: creating large prime product");
      fflush (out);
      create_batch_product (P, B, L, pol);
      goto end;
    }

  fp = fopen (f, "r");
  if (fp != NULL) /* case 3 */
    {
      fprintf (out, "# batch: reading large prime product");
      fflush (out);
      ret = mpz_inp_raw (P, fp);
      if (ret == 0)
        {
          fprintf (stderr, "Error while reading batch product from %s\n", f);
          exit (1);
        }
      goto end;
    }

  /* case 2 */
  fprintf (out, "# batch: creating large prime product");
  fflush (out);
  create_batch_product (P, B, L, pol);

  fp = fopen (f, "w");
  ASSERT_ALWAYS(fp != NULL);

  ret = mpz_out_raw (fp, P);
  if (ret == 0)
    {
      fprintf (stderr, "Error while writing batch product to %s\n", f);
      exit (1);
    }

  fclose (fp);

 end:
  gmp_fprintf (out, " of %zu bits took %.2fs (wct %.2fs)\n",
               mpz_sizeinbase (P, 2), seconds () - s, wct_seconds () - wct);
  fflush (out);
}

