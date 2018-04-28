#include "cado.h"
#include <stdio.h>
#include <math.h>
#ifdef  HAVE_OPENMP
#include <omp.h>
#endif
#include <vector>
#include <list>
#include <sstream>
#include "batch.h"
#include "utils.h"

#define STATUS_SMOOTH  0
#define STATUS_UNKNOWN 1

using namespace std;

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
      t->l = (mpz_t *) malloc (sizeof (mpz_t));
      mpz_init_set_ui (t->l[0], n);
      t->n = (unsigned long *) malloc (sizeof (unsigned long));
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
              t->l = (mpz_t *) realloc (t->l, (t->size + 1) * sizeof (mpz_t));
              mpz_init_set_ui (t->l[t->size], 1);
              t->n = (unsigned long *) realloc (t->n, (t->size + 1) * sizeof (unsigned long));
              t->n[t->size] = 0;
              t->size++;
            }
          if (t->n[i+1] == 0)
            mpz_swap (t->l[i+1], t->l[i]);
          else /* accumulate */
            mpz_mul (t->l[i+1], t->l[i+1], t->l[i]);
          t->n[i+1] += t->n[i];
          /* primes from l[i] are now in l[i+1], thus reset l[i] to 1: */
          mpz_set_ui (t->l[i], 1);
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

static void
prime_list (std::vector<unsigned long> & L, prime_info pi,
            unsigned long pmax)
{
  for (unsigned long p = 2 ; p <= pmax ; p = getprime_mt(pi))
      L.push_back(p);
}

static void
prime_list_poly (std::vector<unsigned long> & L, prime_info pi,
                 unsigned long pmax, mpz_poly f)
{
  if (f->deg == 1)
    return prime_list (L, pi, pmax);

  for (unsigned long p = 2 ; p <= pmax ; p = getprime_mt(pi))
    if (mpz_divisible_ui_p (f->coeff[f->deg], p) ||
        mpz_poly_roots_ulong (NULL, f, p) > 0)
        L.push_back(p);
}

/* add in the product tree L all primes pmin <= p < pmax.
   Assume pmin is the current prime in 'pi'
   (pmin=2 when 'pi' was just initialized).
   Return the current prime in 'pi' at the end, i.e.,
   the smallest prime >= pmax. */
static void
prime_tree (mpz_product_tree L, unsigned long pmax)
{
  prime_info pi;
  prime_info_init (pi);
  for (unsigned long p = 2; p < pmax; p = getprime_mt (pi))
    mpz_product_tree_add_ui (L, p);
  prime_info_clear (pi);
}

/* same as prime_tree, but only adds primes p for which f has at least one
   root modulo p, or the leading coefficient of f vanishes modulo p */
static void
prime_tree_poly (mpz_product_tree L, unsigned long pmax, mpz_poly_srcptr f)
{
  if (f->deg == 1) {
    prime_tree (L, pmax);
    return;
  }

  int nthreads = 1;
#ifdef HAVE_OPENMP
#pragma omp parallel
  nthreads = omp_get_num_threads ();
#endif
  nthreads *= 10; /* to amortize the varying cost of mpz_poly_roots_ulong */

  unsigned long * q = (unsigned long *) malloc (nthreads * sizeof (unsigned long));

  prime_info pi;
  prime_info_init (pi);
  for (unsigned long p = 2; p < pmax;)
    {
        int i;
      /* sequential part: getprime_mt is fast */
      for (i = 0; i < nthreads && p < pmax; p = getprime_mt (pi), i++)
        q[i] = p;

      /* parallel part: mpz_poly_roots_ulong is the bottleneck */
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (int j = 0; j < i; j++)
        if (mpz_divisible_ui_p (f->coeff[f->deg], q[j]) == 0 &&
            mpz_poly_roots_ulong (NULL, f, q[j]) == 0)
          q[j] = 0;

      /* sequential part: mpz_product_tree_add_ui is fast */
      for (int j = 0; j < i; j++)
        if (q[j])
          mpz_product_tree_add_ui (L, q[j]);
    }

  prime_info_clear (pi);
  free (q);
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
  l->a = (int64_t*) realloc (l->a, newsize * sizeof (int64_t));
  l->b = (uint64_t*) realloc (l->b, newsize * sizeof (uint64_t));
  l->R = (mpz_t*) realloc (l->R, newsize * sizeof (mpz_t));
  l->A = (mpz_t*) realloc (l->A, newsize * sizeof (mpz_t));
  l->R0 = (mpz_t*) realloc (l->R0, newsize * sizeof (mpz_t));
  l->A0 = (mpz_t*) realloc (l->A0, newsize * sizeof (mpz_t));
  l->sq = (mpz_t*) realloc (l->sq, newsize * sizeof (mpz_t));
  l->side = (int*) realloc (l->side, newsize * sizeof (int));
  l->perm = (uint32_t*) realloc (l->perm, newsize * sizeof (uint32_t));
  l->alloc = newsize;
  if (newsize < l->size)
    l->size = newsize;
}

void
cofac_list_add (cofac_list l, long a, unsigned long b,
        mpz_srcptr R, mpz_srcptr A,
                int side, mpz_srcptr sq)
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
prime_product_poly (mpz_t P, unsigned long p_max, mpz_poly_srcptr f)
{
  mpz_product_tree L;

  mpz_product_tree_init (L);
  prime_tree_poly (L, p_max, f);

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
  nbits = (unsigned long**) malloc ((h + 1) * sizeof (unsigned long*));
  for (i = 0; i <= h; i++)
    {
      nbits[i] = (unsigned long *) malloc (w[i] * sizeof (unsigned long));
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

  if (mpz_cmp_ui(P, 1) == 0) {
      /* XXX do we have something to do with perm[] ? */
      return;
  }

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
        /* Divide out R by gcd(P, R) as much as we can. The first gcd may
         * have some cost, the subsequent ones are supposedly cheap
         * enough */
        for(;;) {
            mpz_gcd (T[0][j], T[0][j], R[perm[j]]);
            if (mpz_cmp_ui(T[0][j], 1) == 0)
                break;
            mpz_divexact (R[perm[j]], R[perm[j]], T[0][j]);
        }
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
  mpz_t BB;
  mpz_init(BB);
  mpz_mul(BB,B,B);

  for (unsigned long j = *nb_smooth; j < *nb_unknown; j++)
    {
      unsigned long i = perm[j];
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

static void
trial_divide (std::vector<cxx_mpz>& factors, cxx_mpz & n, std::vector<unsigned long> const& SP)
{
    for (auto p : SP) {
        while (mpz_divisible_ui_p (n, p)) {
            factors.push_back(p);
            mpz_divexact_ui (n, n, p);
        }
    }
}

/* Puts in factors[] the prime factors of n. Additional info provided:
 *  - cofac must be a divisor of n (possibly composite)
 *  - sq, if not null, must be a prime divisor of n.
 *
 * The list SP contains small primes (less than B).
 * B is the small prime bound: any factor < B^2 is necessarily prime.
 */
static bool
factor_simple_minded (std::vector<cxx_mpz> &factors,
              cxx_mpz & n,
              facul_method_t *methods,
              unsigned int lpb, double B,
              std::vector<unsigned long> const& SP,
              cxx_mpz& cofac, mpz_ptr sq)
{
    double BB = B * B, BBB = B * B * B;
    ASSERT(mpz_divisible_p (n, cofac));
    mpz_divexact (n, n, cofac);

    if (sq) {
        ASSERT(mpz_divisible_p (n, sq));
        mpz_divexact (n, n, sq);
    }

    /* remove small primes */
    trial_divide (factors, n, SP);

    /* the cofactor that we found with the product tree will have primes
     * between lim and 2^batchlpb. We typically have batchlpb < lpb, and
     * lim > sqrt(2^lpb), so that we know that [cofac] has no prime factor
     * below B=sqrt(2^lpb). It is not guaranteed, though: we may have
     * elected to get rid of sieving on that side, in which case [cofac]
     * may very well contain small primes.
     */
    trial_divide (factors, cofac, SP);

    std::list<std::pair<cxx_mpz, facul_method_t *>> composites;
    if (mpz_cmp_ui(n, 1) > 0) 
        composites.push_back(std::make_pair(std::move(n), methods));
    if (mpz_cmp_ui(cofac, 1) > 0)
        composites.push_back(std::make_pair(std::move(cofac), methods));

    struct modset_t fm[2];

    for (; !composites.empty() ; ) {
        cxx_mpz & n0 = composites.front().first;
        facul_method_t * pm = composites.front().second;
        if (mpz_cmp_d (n0, BB) < 0) {
            if (mpz_cmp_ui(n0, 1) > 0)
                factors.push_back(move(n0));
            composites.pop_front();
            continue;
        }

        if (!pm->method) {
            mpz_set(cofac, n0);
            return false;
        }

        /* We're not going to try this method mode than once */
        fm[0].arith = modset_t::CHOOSE_NONE;
        fm[1].arith = modset_t::CHOOSE_NONE;

        std::vector<cxx_mpz> temp;
        int nf = facul_doit_onefm_mpz (temp, n0, *pm, &fm[0], &fm[1], lpb, BB, BBB);
        pm++;

        /* Could happen if we allowed a cofactor bound after batch
         * cofactorization */
        if (nf == FACUL_NOT_SMOOTH) {
            mpz_set(cofac, n0);
            return false;
        }

        /* In this case, no prime factor was stored, no composite was
         * stored: the input number has not been changed, we move on to
         * the next method */
        if (nf == 0) {
            composites.front().second = pm;
            continue;
        }

        /* factors[] contains the primes, and fm[] the composite
         * cofactors found */
        composites.pop_front();

        ASSERT_ALWAYS(temp.size() == (size_t) nf);
        /* temp[0..nf-1] are prime factors of n, 0 <= nf <= 2 */
        for (int j = 0; j < nf; j++) {
            factors.push_back(std::move(temp[j]));
        }

        /* we may also have composites */
        for (int j = 0; j < 2; j++) {
            if (fm[j].arith == modset_t::CHOOSE_NONE) continue;

            /* fm is a non-trivial composite factor */
            cxx_mpz t;
            modset_get_z (t, &fm[j]);

            /* t should be composite, i.e., t >= BB */
            ASSERT(mpz_cmp_d (t, BB) >= 0);

            composites.push_back(std::make_pair(std::move(t), pm));
        }
    }

    if (sq) {
        cxx_mpz tz;
        mpz_set(tz, sq);
        factors.push_back(std::move(tz));
    }

    return true;
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
static bool
factor_one (cofac_list L, cado_poly pol, unsigned long *lim, int * batchlpb, int *lpb,
            FILE *out, facul_method_t *methods,
            std::vector<unsigned long> (&SP)[2],
            unsigned long i)
{
    uint32_t *perm = L->perm;
    std::ostringstream os;

    int64_t a = L->a[perm[i]];
    uint64_t b = L->b[perm[i]];

    os << a << "," << b;

    bool smooth = true;
    std::vector<cxx_mpz> factors[2];
    cxx_mpz norm, cofac;
    int side;
    for(side = 0 ; smooth && side < 2 ; side++) {
        mpz_set(cofac,(side ? L->A0 : L->R0)[perm[i]]);
        mpz_poly_homogeneous_eval_siui (norm, pol->pols[side], a, b);
        smooth = smooth && factor_simple_minded (factors[side], norm, methods,
                lpb[side], (double) lim[side], SP[side],
                cofac,
                (L->side[perm[i]] == side) ? L->sq[perm[i]] : NULL);
    }

    /* avoid two threads writing a relation simultaneously */
#ifdef  HAVE_OPENMP
#pragma omp critical
#endif
    {
        if (!smooth) {
            --side;
            /* when we've knowingly decided to _do_ some cofactoring
             * after the product-tree on that side, then it's normal to
             * have non-smooth values after all.
             */
            if (batchlpb[side] == lpb[side]) {
                gmp_fprintf(out,
                        "# failed on %s on side %d; non-smooth cofactor %Zd\n",
                        os.str().c_str(), side, (mpz_srcptr) cofac);
            }
        } else {
            for(int side = 0 ; smooth && side < 2 ; side++) {
                os << ":";
                auto it = factors[side].begin();
                os << std::hex << *it++;
                for( ; it < factors[side].end() ; ++it) {
                    os << "," << std::hex << *it;
                }
            }
            fprintf (out, "%s\n", os.str().c_str());
        }
        fflush (out);
    }

    return smooth;
}

/* Given a list L of bi-smooth cofactors, print the corresponding relations
   on "out".
   n is the number of bi-smooth cofactors in L.
*/
unsigned long
factor (cofac_list L, unsigned long n, cado_poly pol, int batchlpb[], int lpb[],
        FILE *out, int nthreads MAYBE_UNUSED)
{
  unsigned long i, B[2];
  int nb_methods;
  facul_method_t *methods;
  std::vector<unsigned long> SP[2];
  prime_info pi;
  double start = seconds (), wct_start = wct_seconds ();

  for(int side = 0 ; side < 2 ; side++) {
      prime_info_init (pi);
      B[side] = (unsigned long) ceil (pow (2.0, (double) lpb[side] / 2.0));
      prime_list_poly (SP[side], pi, B[side], pol->pols[side]);
      prime_info_clear (pi);
  }

  nb_methods = 30;
  if (nb_methods >= NB_MAX_METHODS)
    nb_methods = NB_MAX_METHODS - 1;
  methods = facul_make_default_strategy (nb_methods - 3, 0);

  unsigned long rep = 0;
#ifdef HAVE_OPENMP
  omp_set_num_threads (nthreads);
#pragma omp parallel for schedule(static)
#endif
  for (i = 0; i < n; i++)
    rep += factor_one (L, pol, B, batchlpb, lpb, out, methods, SP, i);

  facul_clear_methods (methods);

  fprintf (out,
          "# batch: took %.2fs (wct %.2fs) to factor %lu smooth relations (%lu final cofac misses)\n",
           seconds () - start, wct_seconds () - wct_start, rep, n-rep);

  return rep;
}

static void
create_batch_product (mpz_t P, unsigned long L, mpz_poly_srcptr pol)
{
  prime_product_poly (P, L, pol);
}

/* output P in the batch file fp, with a header consisting of 3 lines:
   1) the factor base bound B
   2) the large prime bound L
   3) the polynomial, in the form "f0 f1 ... fd"
   Then the integer P is written using mpz_out_raw.
   The header can be read by a human with head -3 batch_file. */
static void
output_batch (FILE *fp, unsigned long B, unsigned long L,
              mpz_poly pol, mpz_t P, const char *f)
{
  int ret;

  ret = fprintf (fp, "%lu\n", B);
  ASSERT_ALWAYS (ret > 0);
  ret = fprintf (fp, "%lu\n", L);
  ASSERT_ALWAYS (ret > 0);
  mpz_poly_fprintf_coeffs (fp, pol, ' ');
  ret = mpz_out_raw (fp, P);
  if (ret == 0)
    {
      fprintf (stderr, "Error while writing batch product to %s\n", f);
      exit (1);
    }
}

/* read a batch file from fp, and check the header is consistent with
   B, L and pol. See #21459. */
static void
input_batch (FILE *fp, unsigned long B, unsigned long L, mpz_poly pol,
             mpz_t P, const char *f)
{
  unsigned long Bread, Lread;
  mpz_poly pol_read;
  int ret;
  char msg[1024];
  msg[0]='\0';

#define CHECK_Z(condition, error_message) do {				\
  if (!(condition)) {							\
      snprintf(msg, sizeof(msg), error_message);			\
      goto parse_error;							\
  }									\
} while (0)
#define CHECK_2(condition, error_message, arg1, arg2) do {		\
  if (!(condition)) {							\
      snprintf(msg, sizeof(msg), error_message, arg1, arg2);		\
      goto parse_error;							\
  }									\
} while (0)
  ret = fscanf (fp, "%lu\n", &Bread);
  CHECK_Z(ret == 1, "Cannot read B\n");
  CHECK_2(Bread == B, "Inconsistent B: expected %lu, file has %lu\n", B, Bread);
  ret = fscanf (fp, "%lu\n", &Lread);
  CHECK_Z(ret == 1, "Cannot read L\n");
  CHECK_2(Lread == L, "Inconsistent L: expected %lu, file has %lu\n", L, Lread);
  mpz_poly_init (pol_read, pol->deg);
  mpz_poly_fscanf_coeffs (fp, pol_read, ' ');
  if (mpz_poly_cmp (pol_read, pol) != 0)
    {
      fprintf (stderr, "Error while reading batch product from %s:\n", f);
      fprintf (stderr, "Inconsistent polynomial in batch file\n");
      fprintf (stderr, "expected ");
      mpz_poly_fprintf (stderr, pol);
      fprintf (stderr, "file has ");
      mpz_poly_fprintf (stderr, pol_read);
      exit (EXIT_FAILURE);
    }
  mpz_poly_clear (pol_read);
  /* now that the header is consistent, we read the integer P */
  ret = mpz_inp_raw (P, fp);
  CHECK_Z(ret > 0, "Could not read large integer\n");
  return;
#undef CHECK_2
#undef CHECK_Z
parse_error:
  fprintf (stderr, "Error while reading batch product from %s:\n%s", f, msg);
  exit(EXIT_FAILURE);
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

  if (L <= B) {
      /* We may be content with having a product tree on one side only.
       * In which case we'll quietly return. The product is over an empty
       * set of primes, so it's best defined as being 1. */
      mpz_set_ui(P, 1);
      return;
  }

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
      create_batch_product (P, L, pol);
      goto end;
    }

  fp = fopen (f, "r");
  if (fp != NULL) /* case 3 */
    {
      fprintf (out, "# batch: reading large prime product");
      fflush (out);
      input_batch (fp, B, L, pol, P, f);
      goto end;
    }

  /* case 2 */
  fprintf (out, "# batch: creating large prime product");
  fflush (out);
  create_batch_product (P, L, pol);

  fp = fopen (f, "w");
  ASSERT_ALWAYS(fp != NULL);

  output_batch (fp, B, L, pol, P, f);

  fclose (fp);

 end:
  gmp_fprintf (out, " of %zu bits took %.2fs (wct %.2fs)\n",
               mpz_sizeinbase (P, 2), seconds () - s, wct_seconds () - wct);
  fflush (out);
}

