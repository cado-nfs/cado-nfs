/* Batch cofactorization.

   References:
   [1] Daniel J. Bernstein, How to find small factors of integers,
       http://cr.yp.to/papers.html#sf, 2002.
   We implement Algorithm 7.1 from [1] here, except that once the smooth
   cofactors have been detected, instead of using another remainder tree
   approach to factor them, we factor them naively. */

#include "cado.h" // IWYU pragma: keep

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>

#include <iterator>
#include <list>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <gmp.h>
#include "fmt/format.h"

#include "batch.hpp"
#include "cado_poly.h"
#include "facul.hpp"
#include "facul_doit.hpp"
#include "facul_method.hpp"
#include "facul_strategies.hpp"
#include "getprime.h"
#include "gmp_aux.h"
#include "special-q.hpp"
#include "macros.h"
#include "modset.hpp"
#include "mpz_poly.h"
#include "omp_proxy.h"
#include "relation.hpp"
#include "rootfinder.h"
#include "timing.h"
#include "utils_cxx.hpp"

/* This function is useful in the openmp context. This segment goes
 * parallel, and all threads except the calling thread subtract their
 * registered RUSAGE_THREAD counter (a.k.a. seconds_thread()) to the
 * provided timer. This must obviously be paired with
 * add_openmp_subtimings().
 *
 * The following construct is ok:
 *
 * double e0 = extra_time;
 *
 * call a function that has a double& extra_time argument
 * 
 *   in there, do subtract_openmp_subtimings()
 *   do something openmp
 *   add_openmp_subtimings()
 * 
 * print extra_time - e0.
 *
 * It is not ok to call subtract_openmp_subtimings twice. Implicitly, a
 * function that has a double& extra_time argument embeds some sort of
 * openmp computation that is enclosed in a subtract/add pair. Therefore,
 * it is also *not* ok to call such a function within a subtract/add
 * pair.
 */
static void subtract_openmp_subtimings(double & extra_time MAYBE_UNUSED)
{
#ifdef HAVE_OPENMP
#pragma omp parallel
  {
      if (omp_get_thread_num() != omp_get_ancestor_thread_num(omp_get_level()-1)) {
#pragma omp critical
          {
              extra_time -= seconds_thread();
          }
      }
  }
#endif
}
static void add_openmp_subtimings(double & extra_time MAYBE_UNUSED)
{
#ifdef HAVE_OPENMP
#pragma omp parallel
  {
      if (omp_get_thread_num() != omp_get_ancestor_thread_num(omp_get_level()-1)) {
#pragma omp critical
          {
              extra_time += seconds_thread();
          }
      }
  }
#endif
}

/* structure to compute on-line a product tree, avoiding to first compute a
   list of mpz_t (which might take too much memory) */
typedef struct {
  mpz_t *l;     /* the value stored is l[0] * l[1] * ... * l[size-1],
                   where l[0] is the product of n[0] elements, l[1] is
                   the product of n[1] elements, ..., with n[0]=0 or 1,
                   n[1]=0 or 2, ..., n[k]=0 or 2^k */
  unsigned long *n;
  size_t size;
} mpz_product_tree_t;
typedef mpz_product_tree_t mpz_product_tree[1];


static void
mpz_product_tree_init (mpz_product_tree t)
{
  t->l = nullptr;
  t->n = nullptr;
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
              checked_realloc(t->l, t->size + 1);
              mpz_init_set_ui (t->l[t->size], 1);
              checked_realloc(t->n, t->size + 1);
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

static void
prime_list (std::vector<unsigned long> & L, prime_info pi,
            unsigned long pmax)
{
  for (unsigned long p = 2 ; p <= pmax ; p = getprime_mt(pi))
      L.push_back(p);
}

static void
prime_list_poly (std::vector<unsigned long> & L, prime_info pi,
                 unsigned long pmax, cxx_mpz_poly const & f)
{
  if (f->deg == 1)
    return prime_list (L, pi, pmax);

  if (pmax < 2)
      return;

  gmp_randstate_t rstate;
  gmp_randinit_default(rstate);

  for (unsigned long p = 2 ; p <= pmax ; p = getprime_mt(pi))
    if (mpz_divisible_ui_p (mpz_poly_lc(f), p) ||
        mpz_poly_roots_ulong (nullptr, f, p, rstate) > 0)
        L.push_back(p);

  gmp_randclear(rstate);
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
 * root modulo p, or the leading coefficient of f vanishes modulo p
 *
 * 
 * Adds to extra_time the cpu time (RUSAGE_THREAD, seconds_thread())
 * spent in openmp helper threads, NOT counting the time spent in the
 * main thread.
 */
static void
prime_tree_poly (mpz_product_tree L, unsigned long pmax, cxx_mpz_poly const& f, double & extra_time)
{
  if (f->deg == 1) {
    prime_tree (L, pmax);
    return;
  }

  int nthreads = 1;

  /* TODO: there are (far) better ways to deal with the parallelization
   * task.  My favorite would be to use subdivide_primes_interval and
   * prime_info_seek (see 12186c0ab when it gets merged)
   */

#ifdef HAVE_OPENMP
#pragma omp parallel
  nthreads = omp_get_num_threads ();
#endif
  nthreads *= 10; /* to amortize the varying cost of mpz_poly_roots_ulong */

  unsigned long * q = (unsigned long *) malloc (nthreads * sizeof (unsigned long));

  subtract_openmp_subtimings(extra_time);

  prime_info pi;
  prime_info_init (pi);

  std::vector<cxx_gmp_randstate> rstate_per_thread(omp_get_max_threads());

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
        if (mpz_divisible_ui_p (mpz_poly_lc(f), q[j]) == 0 &&
            mpz_poly_roots_ulong (nullptr, f, q[j], rstate_per_thread[omp_get_thread_num()]) == 0)
          q[j] = 0;

      /* sequential part: mpz_product_tree_add_ui is fast */
      for (int j = 0; j < i; j++)
        if (q[j])
          mpz_product_tree_add_ui (L, q[j]);
    }

  prime_info_clear (pi);

  add_openmp_subtimings(extra_time);

  free (q);
}

/* put in P the product of primes p for which the given polynomial has factors
   modulo p */
static void
prime_product_poly (mpz_t P, unsigned long p_max, cxx_mpz_poly const & f, double & extra_time)
{
  mpz_product_tree L;

  mpz_product_tree_init (L);
  prime_tree_poly (L, p_max, f, extra_time);

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
 * Put in w[i] the number of elements of level i:
 * w[0] = n, w[1] = ceil(n/2), ...
 * 
 * Adds to extra_time the cpu time (RUSAGE_THREAD, seconds_thread())
 * spent in openmp helper threads, NOT counting the time spent in the
 * main thread.
 */
static mpz_t**
product_tree (std::vector<cxx_mpz> const & R, size_t *w, double & extra_time)
{
  size_t const n = R.size();
  int const h = tree_height (n);
  ASSERT_ALWAYS(n >= 1);

  mpz_t ** T = (mpz_t**) malloc ((h + 1) * sizeof (mpz_t*));

  /* initialize tree */
  w[0] = n;
  for (int i = 0; i <= h; i++)
    {
      w[i] = 1 + ((n - 1) >> i);
      T[i] = (mpz_t*) malloc (w[i] * sizeof (mpz_t));
      for (size_t j = 0; j < w[i]; j++)
        mpz_init (T[i][j]);
    }

  /* initialize T[0] to R */
  for (size_t j = 0 ; j < R.size() ; ++j)
    mpz_set (T[0][j], R[j]);

  subtract_openmp_subtimings(extra_time);

  /* compute product tree */
  for (int i = 1; i <= h; i++)
    {
#ifdef HAVE_OPENMP
#pragma omp parallel for schedule(static)
#endif
      for (size_t j = 0; j < w[i-1] / 2; j++)
        mpz_mul (T[i][j], T[i-1][2*j], T[i-1][2*j+1]);
      if (w[i-1] & 1)
        mpz_set (T[i][w[i]-1], T[i-1][w[i-1]-1]);
    }

  add_openmp_subtimings(extra_time);

  return T;
}

/* Auxiliary routine: a node T[i][j] with left son T[i-1][2*j] and right
   son T[i-1][2*j+1] corresponds to a node say t with left son u and right
   son v in the original product tree, where t = u*v.
   Now the current value of T[i][j], say t', is an approximation of
   (P*2^k)/t mod 2^k, where k = nbits(t) + g (g is the guard).
   This routine puts in T[i-1][2*j] and T[i-1][2*j+1] an approximation
   u' of (P*2^ku)/u mod 2^ku and v' of (P*2^kv)/v mod 2^kb respectively, where
   ku = nbits(u) + g, and kv = nbits(v) + g.
   Assume at input we have |t' - (P*2^k)/t| < x mod 2^k, where "mod 2^k"
   means that all quantities are reduced modulo 2^k with remainder in
   [-2^k/2, 2^k/2]. Then we first compute v*t':
   |v*t' - (P*2^k)/u| < v*x mod 2^k   [since t = u*v]
   then we divide by 2^(k-ku):
   |v*t'/2^(k-ku) - (P*2^ku)/u| < v*x/2^(k-ku) mod 2^ku
   thus since v = t/u < 2^k/2^(ku-1) and |v*t'/2^(k-ku) - u'| < 1,
   |u' - (P*2^ku)/u| < 2*x+1 mod 2^ku
   This proves that the error goes from x to 2x+1 at each step
   of the tree, thus since it is at most 1 at the root of the tree, it is
   at most 2^(h+1)-1 at the leaves. */
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
 * (http://cr.yp.to/arith/scaledmod-20040820.pdf).
 * At the root, we compute a floating-point
 * approximation of P/T[h][0] with m+guard bits, where m = nbits(T[h][0]).
 *
 * Adds to extra_time the cpu time (RUSAGE_THREAD, seconds_thread())
 * spent in openmp helper threads, NOT counting the time spent in the
 * main thread.
 */
static void
remainder_tree (mpz_t **T, size_t *w, mpz_srcptr P,
        std::vector<cxx_mpz> & R,
        double & extra_time)
{
  size_t const n = R.size();
  unsigned long h = tree_height (n), i, j, guard;
  unsigned long **nbits;
  mpz_t Q;

  guard = h + 2; /* see error analysis above */
  nbits = (unsigned long**) malloc ((h + 1) * sizeof (unsigned long*));
  for (i = 0; i <= h; i++)
    {
      nbits[i] = (unsigned long *) malloc (w[i] * sizeof (unsigned long));
      for (j = 0; j < w[i]; j++)
        nbits[i][j] = mpz_sizeinbase (T[i][j], 2);
    }

  subtract_openmp_subtimings(extra_time);

  mpz_init (Q);
  mpz_mod (Q, P, T[h][0]); /* first reduce modulo T[h][0] in case P is huge */
  mpz_mul_2exp (Q, Q, nbits[h][0] + guard);
  mpz_tdiv_q (T[h][0], Q, T[h][0]);
  /* |T' - 2^k*P/T| < 1 mod 2^k, where T is the original value of T[h][0],
     T' is the new value of T[h][0], k = nbits(T) + guard, and "mod 2^k"
     means that all values are taken modulo 2^k, with remainder in
     [-2^k/2,2^k/2]. */
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

  /* now for all leaves, if R = R[j], and R' = T[0][j],
     we have |R' - 2^k*P/R| < 2^h mod 2^k, where k = nbits(R) + g,
     thus R' = 2^k*P/R + a*2^k + b, with a integer, and |b| < 2^(h+1),
     thus R*R' = 2^k*P + a*R*2^k + b*R
     thus P mod R = R*R'/2^k - b*R/2^k, with |b*R/2^k| < 2^(h+1-g).
     Now it suffices to have g >= h+2 so that the 2^(h+1-g) term
     is less than 1/2, and rounding R*R'/2^k to the nearest integer
     gives P mod R. */

  /* from T[0][j] ~ P/R[j]*2^(nbits[0][j] + guard) mod 2^(nbits[0][j] + guard),
     get T[0][j]*R[j]/2^(nbits[0][j] + guard) ~ P mod R[j] */
  ASSERT_ALWAYS(R.size() == w[0]);
  for (size_t j = 0 ; j < R.size() ; ++j)
    {
      mpz_mul (T[0][j], T[0][j], R[j]);
      /* T[0][j] ~ P*2^(nbits[0][j] + guard) mod R[j]*2^(nbits[0][j]+guard) */
      mpz_div_2exp (T[0][j], T[0][j], nbits[0][j]);
      /* T[0][j] ~ P*2^guard mod R[j]*2^guard */
      mpz_add_ui (T[0][j], T[0][j], (1UL << guard) - 1UL);
      mpz_div_2exp (T[0][j], T[0][j], guard);
    }

  mpz_clear (Q);

  add_openmp_subtimings(extra_time);

  for (i = 0; i <= h; i++)
    free (nbits[i]);
  free (nbits);
}

/* Clear the product tree. */
static void
clear_product_tree (mpz_t **T, unsigned long n, size_t *w)
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
 * R[0], ..., R[n-1] are cofactors
 * P is the product of primes
 * Output:
 * Each R[j] has been divided by its P-smooth part.
 *
 * Adds to extra_time the cpu time (RUSAGE_THREAD, seconds_thread())
 * spent in openmp helper threads, NOT counting the time spent in the
 * main thread.
 */
static void
smoothness_test (std::vector<cxx_mpz> & R, mpz_srcptr P, FILE *out, double& extra_time)
{
  size_t w[MAX_DEPTH];
  mpz_t **T;
  double s, st, e0, wct;
  size_t const n = R.size();

  if (n == 0)
    return;

  if (mpz_cmp_ui(P, 1) == 0) {
      /* XXX do we have something to do with perm[] ? */
      return;
  }

  s = seconds ();
  st = seconds_thread ();
  e0 = extra_time;
  wct = wct_seconds ();

  T = product_tree (R, w, extra_time);
  size_t const h = tree_height (n);
  fprintf (out, "# batch: took %.2fs (%.2fs + %.2fs ; wct %.2fs) to compute product tree of %zu bits\n",
           seconds() - s,
           seconds_thread () - st,
           extra_time - e0,
           wct_seconds () - wct, mpz_sizeinbase (T[h][0], 2));

  /* compute remainder tree */
  s = seconds ();
  st = seconds_thread ();
  e0 = extra_time;
  wct = wct_seconds ();

  remainder_tree (T, w, P, R, extra_time);

  fprintf (out, "# batch: took %.2fs (%.2fs + %.2fs ; wct %.2fs) to compute remainder tree\n",
           seconds() - s,
           seconds_thread () - st,
           extra_time - e0,
           wct_seconds () - wct);


  /* now T[0][j] = P mod R[j] for 0 <= j < n */
  ASSERT_ALWAYS(R.size() == w[0]);
  for (size_t j = 0 ; j < R.size() ; ++j) {
        /* Divide out R by gcd(P, R) as much as we can. The first gcd may
         * have some cost, the subsequent ones are supposedly cheap
         * enough */
        for(;;) {
            mpz_gcd (T[0][j], T[0][j], R[j]);
            if (mpz_cmp_ui(T[0][j], 1) == 0)
                break;
            mpz_divexact (R[j], R[j], T[0][j]);
        }
    }

  clear_product_tree (T, n, w);
}

/* return the number n of smooth relations in l (same as l.size()) */
size_t
find_smooth (std::list<cofac_candidate> & l,
        std::vector<cxx_mpz> const & batchP,
        std::vector<unsigned int> const & batchlpb,
        std::vector<unsigned int> const & lpb,
        std::vector<unsigned int> const & batchmfb,
        FILE *out,
        int nthreads MAYBE_UNUSED, double & extra_time)
{
    int const nsides = lpb.size();
    ASSERT_ALWAYS(batchP.size() == (size_t) nsides);
    ASSERT_ALWAYS(batchlpb.size() == (size_t) nsides);
    ASSERT_ALWAYS(batchmfb.size() == (size_t) nsides);
#ifdef HAVE_OPENMP
    omp_set_num_threads (nthreads);
#endif

    /* FIXME FIXME FIXME. We have some static side-decision logic here,
     * and that doesn't make sense if we want to be really side-agnostic
     * in the code. It doesn't make sense if the side of the largest 
     */
    int first_smoothness_test_side;

    if (nsides == 2) {
        /* it seems faster to start from the algebraic side */
        first_smoothness_test_side = 1;
    } else {
        first_smoothness_test_side = 0;
    }

    for (int xside = 0 ; xside < nsides ; ++xside)
    {
        double s, st, e0, wct;

        s = seconds ();
        st = seconds_thread ();
        e0 = extra_time;
        wct = wct_seconds ();

        size_t const input_candidates = l.size();

        int const side = xside ^ first_smoothness_test_side;

        cxx_mpz B, BB, L, M;

        mpz_ui_pow_ui(B, 2, batchlpb[side]);
        mpz_ui_pow_ui(L, 2, lpb[side]);
        mpz_ui_pow_ui(M, 2, batchmfb[side]);
        mpz_mul(BB,B,B);

        /* The post-tree cofactor is useless. Therefore, we:
         *  - copy the cofactors to an auxiliary structure,
         *  - tamper with these auxiliary integers
         *  - keep or discard the original cofactors, depending on
         *  whether we've found that something is smooth.
         */

        std::vector<cxx_mpz> temp;
        temp.reserve(l.size());
        for(auto const & x : l) {
            cxx_mpz const & c(x.cofactor[side]);
            /* If a cofactor is marked as zero, it means "not smooth".
             * This might be the case if the mfb check didn't succeed
             * after sieving for one of the sides > 0 (hence only in an
             * MNFS setting)
             */
            if (c == 0) continue;
            temp.push_back(c);
        }

        smoothness_test (temp, batchP[side], out, extra_time);

        size_t smooth = 0;
        auto jt = begin(temp);
        for(auto it = begin(l) ; it != end(l) ; /* it++ within loop */) {
            /* If we read a zero cofactor in the input set, then it did
             * not enter the product tree. Therefore, we must skip it.
             */
            for( ; it != end(l) && it->cofactor[side] == 0 ; ++it);
            if (it == end(l))
              break;
            /* check if the cofactor on the side we've just tested is
             * smooth. If it isn't, we put it at the end of the array,
             * and we free it.
             *
             * on input, we know that cofactors have no factor<=B. 
             * We're willing to accept prime cofactors <= L,
             * and to try harder to factor composites which are <= M
             *
             * relation i is smooth iff R[i]=1 ; another option is in case
             * the remaining cofactor is below the mfb we've been given.
             */
            cxx_mpz const& c(*jt++);

            if (mpz_cmp_ui (c, 1) == 0
                    || mpz_cmp(c, L) <= 0
                    || (mpz_cmp(c, BB) >= 0 && mpz_cmp(c, M) <= 0
                        && !mpz_probab_prime_p(c, 1)
                       ))
            {
                /* cofactor is smooth, we keep it where it is.  */
                ++it;
                smooth++;
            } else {
                /* cofactor is not smooth. Get rid of it if we're dealing
                 * with the first of the sides to consider, or if there
                 * are more than two sides. but otherwise
                 * just set it to 0, and use that as a marker for
                 * "non-smooth". We do so because in an MNFS context,
                 * further sides might be smooth, and match with the
                 * first side.
                 */
                if (xside == 0 || nsides == 2) {
                    auto nit = it++;
                    l.erase(nit);
                } else {
                    it++->cofactor[side] = 0;
                }
            }
        }

        fprintf (out, "# batch (side %d): took %.2fs (%.2fs + %.2fs ; wct %.2fs) to detect %zu smooth relations out of %zu\n",
                side,
                seconds() - s,
                seconds_thread () - st,
                extra_time - e0,
                wct_seconds () - wct,
                smooth, input_candidates);
    }

    return l.size();
}

size_t
find_smooth (std::list<std::pair<special_q, std::list<cofac_candidate>>> & L,
        std::vector<cxx_mpz> const & batchP,
        std::vector<unsigned int> const & batchlpb,
        std::vector<unsigned int> const & lpb,
        std::vector<unsigned int> const & batchmfb,
        FILE *out,
        int nthreads MAYBE_UNUSED, double & extra_time)
{
    size_t n = 0;
    std::list<std::pair<special_q, std::list<cofac_candidate>>> R;
    for( ; !L.empty() ; ) {
        auto M = std::move(L.front());
        L.pop_front();
        const size_t m = find_smooth(M.second, batchP, batchlpb, lpb, batchmfb, out, nthreads, extra_time);
        n += m;
        if (m)
            R.emplace_back(std::move(M));
    }
    std::swap(L, R);
    return n;
}


static void
trial_divide (std::vector<cxx_mpz>& factors, cxx_mpz & n, std::vector<unsigned long> const& SP)
{
    for (auto p : SP) {
        while (mpz_divisible_ui_p (n, p)) {
            factors.emplace_back(p);
            mpz_divexact_ui (n, n, p);
        }
    }
}

/* Puts in factors[] the prime factors of n. Additional info provided:
 *  - cofac must be a divisor of n (possibly composite)
 *  - sq_factors, is the (possibly empty if on the wrong side) list of
 *    prime factors of the special-q 
 *
 * The list SP contains small primes (less than B).
 * B is the small prime bound: any factor < B^2 is necessarily prime.
 */
static bool
factor_simple_minded (std::vector<cxx_mpz> &factors,
              cxx_mpz & n,
              std::vector<facul_method> const & methods,
              unsigned int lpb, double B,
              std::vector<unsigned long> const& SP,
              cxx_mpz& cofac,
              std::vector<uint64_t> const& sq_factors)
{
    const double BB = B * B, BBB = B * B * B;

    if (mpz_cmp_ui(cofac, 0) == 0) {
        for (auto sqf : sq_factors) {
            cxx_mpz Sqf;
            mpz_set_uint64(Sqf, sqf);
            ASSERT(mpz_divisible_p (n, Sqf));
            mpz_divexact (n, n, Sqf);
        }
        trial_divide (factors, n, SP);
    } else {
        ASSERT(mpz_divisible_p (n, cofac));
        mpz_divexact (n, n, cofac);

        for (auto sqf : sq_factors) {
            cxx_mpz Sqf;
            mpz_set_uint64(Sqf, sqf);
            ASSERT(mpz_divisible_p (n, Sqf));
            mpz_divexact (n, n, Sqf);
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
    }

    std::list<std::pair<cxx_mpz, std::vector<facul_method>::const_iterator>> composites;
    if (mpz_cmp_ui(n, 1) > 0) 
        composites.emplace_back(std::move(n), methods.begin());
    if (mpz_cmp_ui(cofac, 1) > 0)
        composites.emplace_back(std::move(cofac), methods.begin());


    /* This calls facul_doit_onefm repeatedly,  until there's no work
     * left.
     *
     *
     * XXX I have the impression that the plain "facul" method does the
     * same, in fact.
     */
    for (; !composites.empty() ; ) {
        cxx_mpz & n0 = composites.front().first;
        auto pm = composites.front().second;
        if (mpz_cmp_d (n0, BB) < 0) {
            if (mpz_cmp_ui(n0, 1) > 0)
                factors.push_back(std::move(n0));
            composites.pop_front();
            continue;
        }

        if (pm == methods.end()) {
            cofac = std::move(n0);
            return false;
        }

        /* If fm[j] is still NULL after the call to
	   facul_doit_onefm_mpz, it means fm[j] has not been set. */

        std::vector<cxx_mpz> temp;
        std::vector<std::unique_ptr<FaculModulusBase>> comp;

        facul_status const nf = facul_doit_onefm (temp, n0, *pm, comp,
                lpb, BB, BBB);
        pm++;

        /* Could happen if we allowed a cofactor bound after batch
         * cofactorization */
        if (nf == FACUL_NOT_SMOOTH) {
            cofac = std::move(n0);
            return false;
        }

        /* In this case, no prime factor was stored, no composite was
         * stored: the input number has not been changed, we move on to
         * the next method */
        if (nf == 0 && comp.empty()) {
            composites.front().second = pm;
            continue;
        }

        /* factors[] contains the primes, and fm[] the composite
         * cofactors found */
        composites.pop_front();

        /* temp[0..nf-1] are prime factors of n, 0 <= nf <= 2 */
        for (auto & c : temp)
            factors.push_back(std::move(c));

        /* we may also have composites */
        for (auto & c : comp) {
            /* fm is a non-trivial composite factor */
            auto t = c->get_z();

            /* t should be composite, i.e., t >= BB */
            ASSERT(mpz_cmp_d (t, BB) >= 0);

            if (mpz_perfect_square_p (t))
            {
                /* Yes, this does happen, at least in the F9_batch test
                 */
                mpz_sqrt (t, t);
                if (mpz_probab_prime_p(t, 1)) {
                    factors.emplace_back(t);
                    factors.emplace_back(std::move(t));
                } else {
                    composites.emplace_back(t, pm);
                    composites.emplace_back(std::move(t), pm);
                }
            }
            else
              composites.emplace_back(std::move(t), pm);
        }
    }

    for (auto sqf : sq_factors) {
        cxx_mpz tz;
        mpz_set_uint64(tz, sqf);
        factors.push_back(std::move(tz));
    }

    return true;
}

#if 0
/* unused */
/* strip integers in l[0..n-1] which do not divide P */
static unsigned long
strip (unsigned long *l, unsigned long n, mpz_t P)
{
  unsigned long i, j;

  for (i = j = 0; i < n; i++)
    if (mpz_divisible_ui_p (P, l[i]))
      l[j++] = l[i];
  return j;
}
#endif

/* sqside = 1 if the special-q is on side 1 (algebraic) */
static bool
factor_one (
        std::list<relation> & smooth,
        cofac_candidate const & C,
        cxx_cado_poly const & cpoly,
        special_q const & doing,
        std::vector<unsigned long> const & lim,
        std::vector<unsigned int> const & batchlpb,
        std::vector<unsigned int> const & lpb,
        FILE *out,
        std::vector<facul_method> const & methods,
        std::vector<std::vector<unsigned long>> const & SP,
        int recomp_norm)
{
    int64_t const a = C.a;
    uint64_t const b = C.b;

    int const nsides = cpoly->nb_polys;

    std::vector<std::vector<cxx_mpz>> factors;

    factors.assign(nsides, {});

    cxx_mpz norm, cofac;
    for(int side = 0 ; side < nsides ; side++) {
        mpz_set(cofac, C.cofactor[side]);
        if (recomp_norm) {
            mpz_poly_homogeneous_eval_siui (norm, cpoly->pols[side], a, b);
            mpz_abs(norm, norm);
        } else {
            if (doing.side == side) {
                mpz_mul(norm, cofac, doing.p);
            } else {
                mpz_set(norm, cofac);
            }
        }
        std::vector<uint64_t> const empty;
        bool const smooth = factor_simple_minded (factors[side], norm, methods,
                lpb[side], (double) lim[side], SP[side],
                cofac,
                (doing.side == side) ? doing.prime_factors : empty);
        if (!smooth) {
            /* when we've knowingly decided to _do_ some cofactoring
             * after the product-tree on that side, then it's normal to
             * have non-smooth values after all.
             */
            if (batchlpb[side] == lpb[side]) {
#ifdef HAVE_OPENMP
#pragma omp critical
#endif
                {
                    std::ostringstream os;
                    os << a << "," << b;
                    gmp_fprintf(out,
                            "# failed on %s on side %d;"
                            " non-smooth cofactor %Zd\n",
                            os.str().c_str(), side, (mpz_srcptr) cofac);
                }
            }
            return false;
        }
    }

    relation rel(a,b);
    for (size_t side = 0; side < std::min(rel.sides.size(), factors.size()); side++) { // FIXME workaround for HARDCODED 2 in relation class
        for (auto const& z : factors[side])
            rel.add(side, z, 0);
    }
    rel.compress();
    smooth.push_back(rel);

    return true;
}

/* Given a list L of bi-smooth cofactors, print the corresponding relations
 * on "out".
 * n is the number of bi-smooth cofactors in L.
 * 
 * Adds to extra_time the cpu time (RUSAGE_THREAD, seconds_thread())
 * spent in openmp helper threads, NOT counting the time spent in the
 * main thread.
 */
std::list<relation>
factor (std::list<cofac_candidate> const & L,
        cxx_cado_poly const & cpoly,
        special_q const & doing,
        std::vector<unsigned int> const & batchlpb,
        std::vector<unsigned int> const & lpb,
        int ncurves,
        FILE *out, int nthreads MAYBE_UNUSED, double& extra_time,
        int recomp_norm)
{
  int const nsides = cpoly->nb_polys;
  ASSERT_ALWAYS(batchlpb.size() == (size_t) nsides);
  ASSERT_ALWAYS(lpb.size() == (size_t) nsides);
  std::vector<unsigned long> B(nsides);
  std::vector<std::vector<unsigned long>> SP(nsides);
  prime_info pi;
  double s, st, e0, wct;

  s = seconds ();
  st = seconds_thread ();
  e0 = extra_time;
  wct = wct_seconds ();

  for(int side = 0 ; side < nsides ; side++) {
      prime_info_init (pi);
      B[side] = (ceil (pow (2.0, (double) lpb[side] / 2.0)));
      if (!recomp_norm) {
          // If not recomp_norm, then the poly file might be fake...
          // This list of primes is rather small anyway.
          prime_list (SP[side], pi, B[side]);
      } else {
          prime_list_poly (SP[side], pi, B[side], cpoly->pols[side]);
      }
      prime_info_clear (pi);
  }

  /* this creates the factoring methods (and B1-adapted addition chains)
   * from the parameters that are given in the default strategy
   */
  std::vector<facul_method> methods;
  for(auto const & mp : facul_strategy_oneside::default_strategy (ncurves))
      methods.emplace_back(mp);

  std::list<relation> smooth;
  std::list<cofac_candidate>::const_iterator it;
  
#ifdef HAVE_OPENMP
  omp_set_num_threads (nthreads);
#endif

  subtract_openmp_subtimings(extra_time);

#ifdef HAVE_OPENMP
#pragma omp parallel private(it)
  {
      std::list<relation> smooth_local;
#else
      std::list<relation> & smooth_local(smooth);
#endif
      for (it = begin(L); it != end(L); ++it) {
#ifdef HAVE_OPENMP
#pragma omp single nowait
#endif
          factor_one (smooth_local, *it, cpoly, doing,
                  B, batchlpb, lpb, out, methods,
                  SP, recomp_norm);
      }
#ifdef HAVE_OPENMP
#pragma omp critical
      smooth.splice(smooth.end(), smooth_local);
  }
#endif

  add_openmp_subtimings(extra_time);

  fprintf (out,
          "# batch: took %.2fs (%.2f + %.2f ; wct %.2fs) to factor %zu smooth relations (%zd final cofac misses)\n",
          seconds() - s,
          seconds_thread () - st,
          extra_time - e0,
          wct_seconds () - wct,
          smooth.size(), L.size()-smooth.size());

  return smooth;
}

static void
create_batch_product (mpz_t P, unsigned long L, cxx_mpz_poly const & cpoly, double & extra_time)
{
  prime_product_poly (P, L, cpoly, extra_time);
}

/* output P in the batch file fp, with a header consisting of 3 lines:
   1) the factor base bound B
   2) the large prime bound L
   3) the polynomial, in the form "f0 f1 ... fd"
   Then the integer P is written using mpz_out_raw.
   The header can be read by a human with head -3 batch_file.

   XXX This must be kept in sync with sieve/inspect-batch-file.pl
*/
static void
output_batch (FILE *fp, unsigned long B, unsigned long L,
              cxx_mpz_poly const & cpoly, cxx_mpz const & P, const char *f)
{
  int ret;

  ret = fprintf (fp, "%lu\n", B);
  ASSERT_ALWAYS (ret > 0);
  ret = fprintf (fp, "%lu\n", L);
  ASSERT_ALWAYS (ret > 0);
  mpz_poly_fprintf_coeffs (fp, cpoly, " ");
  ret = mpz_out_raw (fp, P);
  if (ret == 0)
    {
      fprintf (stderr, "Error while writing batch product to %s\n", f);
      exit (1);
    }
}

/* read a batch file from fp, and check the header is consistent with
   B, L and cpoly. See #21459. */
static void
input_batch (FILE *fp, unsigned long B, unsigned long L, cxx_mpz_poly const & cpoly,
             cxx_mpz & P, const char *f)
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
  mpz_poly_init (pol_read, cpoly->deg);
  mpz_poly_fscanf_coeffs (fp, pol_read, " ");
  if (mpz_poly_cmp (pol_read, cpoly) != 0)
    {
      fprintf (stderr, "Error while reading batch product from %s:\n", f);
      fprintf (stderr, "Inconsistent polynomial in batch file\n");
      fprintf (stderr, "expected ");
      mpz_poly_fprintf (stderr, cpoly);
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
  throw std::runtime_error(fmt::format(
              "Error while reading batch product from {}:\n{}\n",
              f, msg));
}

/* We have 3 cases:
   1) if f == NULL: P is computed but not stored
   2) if f != NULL but file is non-existing: P is computed and saved in f
   3) if f != NULL and file is existing: P is read from file
*/
void
create_batch_file (std::string const & fs, cxx_mpz & P, unsigned long B, unsigned long L,
                   cxx_mpz_poly const & cpoly, FILE *out, int nthreads MAYBE_UNUSED, double & extra_time)
{
  double e0, s, st, wct;
  const char * f = fs.empty() ? nullptr : fs.c_str();

  if (L <= B) {
      /* We may be content with having a product tree on one side only.
       * In which case we'll quietly return. The product is over an empty
       * set of primes, so it's best defined as being 1. */
      mpz_set_ui(P, 1);
      return;
  }

  // the product of primes up to B takes \log2(B)-\log\log 2 / \log 2
  // bits. The added constant is 0.5287.
  if (log2(B) - log2(GMP_LIMB_BITS) + 0.5287 >= 31) {
    /* mpz_t "size" field is an int, thus can hold up to 2^31-1:
       on a 32-bit processor, this means up to 2^36 bits,
       on a 64-bit processor, this means up to 2^37 bits */
    fprintf(stderr, "Gnu MP cannot deal with prime products that large (maximum 2^%d bits)\n", 36 + (GMP_LIMB_BITS > 32));
      abort();
  } else if (log2(B) + 0.5287 >= 34) {
      fprintf(stderr, "Gnu MP's mpz_inp_raw and mpz_out_raw functions are limited to integers of at most 34 bits\n");
      abort();
  }

#ifdef HAVE_OPENMP
  omp_set_num_threads (nthreads);
#endif

  e0 = extra_time;
  s = seconds ();
  st = seconds_thread ();
  wct = wct_seconds ();

  if (f == nullptr) { /* case 1 */
      fprintf (out, "# batch: creating large prime product");
      fflush (out);
      create_batch_product (P, L, cpoly, extra_time);
  } else {
      auto fp = fopen_helper(f, "r", true);
      if (fp) { /* case 3 */
          fprintf (out, "# batch: reading large prime product");
          fflush (out);
          input_batch (fp.get(), B, L, cpoly, P, f);
      } else {
          /* case 2 */
          fprintf (out, "# batch: creating large prime product");
          fflush (out);
          create_batch_product (P, L, cpoly, extra_time);
          output_batch (fopen_helper (f, "w").get(), B, L, cpoly, P, f);
      }
  }


  gmp_fprintf (out, " of %zu bits took %.2fs (%.2fs + %.2fs ; wct %.2fs)\n",
               mpz_sizeinbase (P, 2),
               seconds() - s,
               seconds_thread () - st,
               extra_time - e0,
               wct_seconds () - wct);

  fflush (out);
}

