#include "cado.h" // IWYU pragma: keep
#include <cstdlib>
#include <cstdio>
#include <climits>
#include <algorithm>
#include <vector>
#include <type_traits>   // for is_same
#include <gmp.h>
#include "cxx_mpz.hpp"   // for cxx_mpz
#include "gmp_aux.h"
#include "macros.h"
#include "mod_ul.h"      // for residueul_t, modul_clear, modul_clearmod
#include "modul_poly.h"
#include "rootfinder.h"
#include "portability.h" // IWYU pragma: keep
#include "gmp_aux.h"

/* Entry point for rootfind routines, for prime p.
   Assume r is an array of deg(F) entries, which are mpz_init'ed. */
int mpz_poly_roots (mpz_t *r, mpz_poly_srcptr F, mpz_srcptr p, gmp_randstate_ptr rstate)
{
    int d = F->deg;

    if (mpz_cmp_ui(p, ULONG_MAX) <= 0) {
        /* There's a chance of using one of our layers. */
        unsigned long pp = mpz_get_ui(p);
        unsigned long * rr;
        int i;
        int n;
	
        if (r == NULL)
            return mpz_poly_roots_ulong (NULL, F, pp, rstate);

        rr = (unsigned long *) malloc(d * sizeof(unsigned long));
        n = mpz_poly_roots_ulong (rr, F, pp, rstate);

        for(i = 0 ; i < n ; i++) {
            /* The assumption is that p fits within an unsigned long
             * anyway. So the roots do as well.
             */
            mpz_set_ui(r[i], rr[i]);
        }
        free(rr);
        return n;
    } else {
      int n;
      n = mpz_poly_roots_mpz (r, F, p, rstate);
      return n;
    }
}


/* put in r[0], ..., r[n-1] the roots of F modulo p, where p is prime,
   and the return value n is the number of roots (without multiplicities) */
int
mpz_poly_roots_ulong (unsigned long *r, mpz_poly_srcptr F, unsigned long p, gmp_randstate_ptr rstate)
{
    int n;

    residueul_t * rr;
    modulusul_t pp;
    modul_initmod_ul(pp, p);
    int i;
    int d = F->deg;

    if (r == NULL)
      return modul_poly_roots(NULL, F, pp, rstate);

    rr = (residueul_t *) malloc(d * sizeof(residueul_t));
    for(i = 0 ; i < d ; i++) {
      modul_init_noset0(rr[i], pp);
    }
    n = modul_poly_roots(rr, F, pp, rstate);
    for(int i = 0 ; i < n ; i++) {
      /* The assumption is that p fits within an unsigned long
       * anyway. So the roots do as well.
       */
      r[i] = modul_get_ul(rr[i], pp);
    }
    for(i = 0 ; i < d ; i++) {
      modul_clear(rr[i], pp);
    }
    free(rr);
    modul_clearmod(pp);

    return n;
}



/* Note that parallelizing this makes no sense. The payload is too small.
 */
int
mpz_poly_roots_uint64 (uint64_t * r, mpz_poly_srcptr F, uint64_t p, gmp_randstate_ptr rstate)
{
    /* This is glue around poly_roots_ulong, nothing more. When uint64
       is larger than ulong, we call mpz_poly_roots_mpz as a fallback */
    unsigned long *rr;
    int i, n;

#if ULONG_BITS < 64
    if (p > (uint64_t) ULONG_MAX)
      {
        mpz_t pp;
        mpz_init (pp);
        mpz_set_uint64 (pp, p);
        if (r == NULL)
          n = mpz_poly_roots_mpz (NULL, F, pp, rstate);
        else
          {
            mpz_t *rr;
            rr = (mpz_t*) malloc ((F->deg + 1) * sizeof (mpz_t));
            for (i = 0; i <= F->deg; i++)
              mpz_init (rr[i]);
            n = mpz_poly_roots_mpz (NULL, F, pp, rstate);
            for (i = 0; i <= F->deg; i++)
              {
                if (i < n)
                  r[i] = mpz_get_uint64 (rr[i]);
                mpz_clear (rr[i]);
              }
            free (rr);
          }
        mpz_clear (pp);
        return n;
      }
#endif

    if (r == NULL)
      return mpz_poly_roots_ulong (NULL, F, p, rstate);

    if (sizeof (unsigned long) != sizeof (uint64_t))
      rr = (unsigned long *) malloc(F->deg * sizeof(unsigned long));
    else
      rr = (unsigned long *) r;
    n = mpz_poly_roots_ulong (rr, F, p, rstate);
    if (sizeof (unsigned long) != sizeof (uint64_t))
      {
        for(i = 0 ; i < n ; i++)
          r[i] = rr[i];
        free (rr);
      }
    return n;
}




typedef int (*sortfunc_t) (const void *, const void *);

static int mpz_poly_coeff_cmp(const mpz_t *a, const mpz_t *b) {
  return mpz_cmp(*a, *b) < 0 ? -1 : 1;
}

/* Assuming f is a (squarefree) product of linear factors mod p, splits it
   and put the corresponding roots mod p in r[]. Return number of roots
   which should be degree of f. Assumes p is odd, and deg(f) >= 1. */
static int
mpz_poly_cantor_zassenhaus (mpz_t *r, mpz_poly_srcptr f, mpz_srcptr p,
                            int depth, gmp_randstate_ptr rstate)
{
  mpz_t a, aux;
  mpz_poly q, h;
  int d = f->deg, dq, n, m;

  mpz_init (a);
  mpz_init (aux);

  /* linear polynomial */
  if (d == 1) {
    mpz_neg (aux, f->coeff[1]);
    mpz_invert (a, aux, p);
    mpz_mul (r[0], a, f->coeff[0]);
    mpz_fdiv_r (r[0], r[0], p);
    n = 1;
    goto clear_a;
  }

  ASSERT_ALWAYS(mpz_odd_p(p));

  /* if f has degree d, then q,h may have up to degree 2d-2 in the
     powering algorithm */
  mpz_poly_init (q, 2 * d - 2);
  mpz_poly_init (h, 2 * d - 2);

  /* random polynomial by a */
  mpz_urandomm(a, rstate, p);
  for (;;)
  {
    /* q=x+a */
    mpz_set_ui (aux, 1);
    mpz_poly_setcoeff (q, 1, aux);
    mpz_poly_setcoeff (q, 0, a);
    q->deg = 1;

    /* h=(x+a)^((p-1)/2) mod (f, p) */
    mpz_sub_ui (aux, p, 1);
    mpz_divexact_ui (aux, aux, 2);
    mpz_poly_pow_mod_f_mod_mpz (h, q, f, aux, p);
    mpz_poly_sub_ui (h, h, 1);

    /* q = gcd(f,h) */
    mpz_poly_gcd_mpz (q, f, h, p);
    dq = q->deg;
    ASSERT (dq >= 0);

    /* recursion-split */
    if (0 < dq && dq < d)
    {
      n = mpz_poly_cantor_zassenhaus (r, q, p, depth + 1, rstate);
      ASSERT (n == dq);

      mpz_poly_divexact (h, f, q, p);
      m = mpz_poly_cantor_zassenhaus (r + n, h, p, depth + 1, rstate);
      ASSERT (m == h->deg);
      n += m;
      break;
    }

    mpz_add_ui (a, a, 1); /* no need to reduce mod p, since it will be done
                             in mpz_poly_pow_mod_f_mod_mpz */
  }

  mpz_poly_clear (q);
  mpz_poly_clear (h);

clear_a:
  mpz_clear (a);
  mpz_clear (aux);
  return n;
}


/* Solve f(x)=0 (mod p), where p is a prime. Return the number of roots.
   Assume d (the degree of f) is at least 1.
 */
int
mpz_poly_roots_mpz (mpz_t *r, mpz_poly_srcptr f, mpz_srcptr p, gmp_randstate_ptr rstate)
{
  int nr = 0;
  mpz_poly fp, g, h;
  int d = f->deg;

  ASSERT(d >= 1);

  mpz_poly_init (fp, d);
  mpz_poly_init (g, 2*d-1);
  mpz_poly_init (h, 2*d-1);

  /* If f has small coefficients (like in Joux-Lercier polynomial selection)
     don't make f monic, since it might make the coefficients of fp blow up.
     In that case we only reduce coefficients in [-p+1, p-1], to keep
     negative coefficients small in absolute value. */
  if (mpz_poly_sizeinbase (f, 2) < mpz_sizeinbase (p, 2))
    mpz_poly_mod_mpz_lazy (fp, f, p);
  else
    mpz_poly_makemonic_mod_mpz (fp, f, p);
  if (fp->deg <= 0)
    goto clear_and_exit;
  /* h=x^p-x (mod mpz_poly_fp) */
  mpz_poly_setcoeff_ui (g, 1, 1);
  mpz_poly_pow_mod_f_mod_mpz (h, g, fp, p, p);
  /* FIXME: instead of computing x^p-x, we could compute x^(p-1) - 1 while
     saving the value of h = x^((p-1)/2). If several roots, gcd(h-1, f)
     might help to split them. */
  mpz_poly_sub(h, h, g);
  /* g = gcd (mpz_poly_fp, h) */
  mpz_poly_gcd_mpz (fp, fp, h, p);
  /* fp contains gcd(x^p-x, f) */
  nr = fp->deg;
  ASSERT (nr >= 0);

  /* If r is NULL, we only return the number of roots. */
  if (r != NULL && nr > 0)
  {
    int n MAYBE_UNUSED = mpz_poly_cantor_zassenhaus (r, fp, p, 0, rstate);
    ASSERT (n == nr);
  }

 clear_and_exit:
  mpz_poly_clear(fp);
  mpz_poly_clear(g);
  mpz_poly_clear(h);

  /* Sort the roots */
  if (r && nr)
    qsort(r, nr, sizeof(mpz_t), (sortfunc_t) &mpz_poly_coeff_cmp);

  return nr;
}



template<typename T>
struct poly_roots_impl_details
{
};

template<>
struct poly_roots_impl_details<unsigned long> {
    static inline void mul(unsigned long& a, unsigned long const & b, unsigned long const &c) { a = b * c; }
    static inline unsigned long mul(unsigned long const & b, unsigned long const &c) { return b * c; }
    static inline void add(unsigned long& a, unsigned long const & b, unsigned long const &c) { a = b + c; }
    static inline unsigned long add(unsigned long const & b, unsigned long const &c) { return b + c; }
    static inline void mod(unsigned long& a, unsigned long const & b, unsigned long const &c) { a = b % c; }
    static inline unsigned long mod(unsigned long const & b, unsigned long const &c) { return b % c; }
    static void invmod(unsigned long& y, unsigned long const & x, unsigned long const & p) {
        cxx_mpz xx,pp,yy;
        mpz_set_ui(xx,x);
        mpz_set_ui(pp,p);
        mpz_invert (yy,xx,pp);
        y = mpz_get_ui(yy);
    }
    static void invmod(unsigned long& y, cxx_mpz const & xx, unsigned long const & p) {
        cxx_mpz pp,yy;
        mpz_set_ui(pp,p);
        mpz_invert (yy,xx,pp);
        y = mpz_get_ui(yy);
    }
    static bool probab_prime_p(unsigned long const & x) { 
        cxx_mpz xx;
        mpz_set_ui(xx, x);
        return mpz_probab_prime_p(xx, 1);
    }
    static inline bool fits(cxx_mpz const & x) { return mpz_fits_ulong_p(x); }
    static inline unsigned long get_from_cxx_mpz(cxx_mpz const & x) { return mpz_get_ui(x); }
};

#ifndef UINT64_T_IS_EXACTLY_UNSIGNED_LONG
template<>
struct poly_roots_impl_details<uint64_t> {
    static inline void mul(uint64_t& a, uint64_t const & b, uint64_t const &c) { a = b * c; }
    static inline uint64_t mul(uint64_t const & b, uint64_t const &c) { return b * c; }
    static inline void add(uint64_t& a, uint64_t const & b, uint64_t const &c) { a = b + c; }
    static inline uint64_t add(uint64_t const & b, uint64_t const &c) { return b + c; }
    static inline void mod(uint64_t& a, uint64_t const & b, uint64_t const &c) { a = b % c; }
    static inline uint64_t mod(uint64_t const & b, uint64_t const &c) { return b % c; }
    static void invmod(uint64_t& y, uint64_t const & x, uint64_t const & p) {
        cxx_mpz xx,pp,yy;
        mpz_set_uint64(xx,x);
        mpz_set_uint64(pp,p);
        mpz_invert (yy,xx,pp);
        y = mpz_get_uint64(yy);
    }
    static void invmod(uint64_t& y, cxx_mpz const & xx, uint64_t const & p) {
        cxx_mpz pp,yy;
        mpz_set_uint64(pp,p);
        mpz_invert (yy,xx,pp);
        y = mpz_get_uint64(yy);
    }
    static bool probab_prime_p(uint64_t const & x) { 
        cxx_mpz xx;
        mpz_set_uint64(xx, x);
        return mpz_probab_prime_p(xx, 1);
    }
    static inline bool fits(cxx_mpz const & x) { return mpz_fits_uint64_p(x); }
    static inline uint64_t get_from_cxx_mpz(cxx_mpz const & x) { return mpz_get_uint64(x); }
};
#endif

template<>
struct poly_roots_impl_details<cxx_mpz> {
    static inline cxx_mpz mul(cxx_mpz const & b, cxx_mpz const &c) {
        cxx_mpz r;
        mul(r,b,c);
        return r;
    }
    static inline void mul(cxx_mpz& a, cxx_mpz const & b, cxx_mpz const &c) {
        mpz_mul(a, b, c);
    }
    static inline void mul(cxx_mpz& a, cxx_mpz const & b, unsigned long const &c) {
        mpz_mul_ui(a, b, c);
    }
#ifndef UINT64_T_IS_EXACTLY_UNSIGNED_LONG
    static inline void mul(cxx_mpz& a, cxx_mpz const & b, uint64_t const &c) {
        mpz_mul_uint64(a, b, c);
    }
#endif
    static inline cxx_mpz add(cxx_mpz const & b, cxx_mpz const &c) {
        cxx_mpz r;
        add(r,b,c);
        return r;
    }
    static inline void add(cxx_mpz& a, cxx_mpz const & b, cxx_mpz const &c) {
        mpz_add(a, b, c);
    }
    static inline cxx_mpz mod(cxx_mpz const & b, cxx_mpz const &c) {
        cxx_mpz r;
        mod(r,b,c);
        return r;
    }
    static inline void mod(cxx_mpz& a, cxx_mpz const & b, cxx_mpz const &c) {
        mpz_mod(a, b, c);
    }
    static void invmod(cxx_mpz& y, cxx_mpz const & x, cxx_mpz const & p) {
        mpz_invert (y,x,p);
    }
    static bool probab_prime_p(cxx_mpz const & x) { 
        return mpz_probab_prime_p(x, 10);
    }
    static inline bool fits(cxx_mpz const &) { return true; }
    static inline cxx_mpz get_from_cxx_mpz(cxx_mpz const & x) { return x; }
};

template<typename T, typename F>
struct poly_roots_impl : public poly_roots_impl_details<T>
{
    typedef poly_roots_impl_details<T> super;
    using super::mul;
    using super::add;
    using super::mod;
    using super::invmod;
    using super::fits;
    using super::probab_prime_p;
    using super::get_from_cxx_mpz;
    std::vector<T> operator()(cxx_mpz_poly const & f, T const & q, std::vector<F> const & qfac, gmp_randstate_ptr rstate);
};

/* Compute the roots of f modulo q, where q is a square-free number whose
 * factorization is given. */
template<typename T, typename F>
std::vector<T> poly_roots_impl<T,F>::operator()(cxx_mpz_poly const & f, T const & q, std::vector<F> const & qfac, gmp_randstate_ptr rstate)
{
    std::vector<T> res;

    if (probab_prime_p(q))
        return mpz_poly_roots (f, q, rstate);

    if (qfac.size() <= 1) {
        /* weird. Then we simply use the prime version, right ? */
        return mpz_poly_roots(f, q, rstate);
    }

    if (!std::is_same<T, cxx_mpz>::value) {
        /* We need that for the CRT to be ok, since we'll accumulate sums
         * in type T (this branch should be optimized away when
         * T==cxx_mpz, with or without the if above, in fact). */
        cxx_mpz test;
        mpz_mul_ui(test, cxx_mpz(q), qfac.size());
        if (!fits(test)) {
            /* compute with type cxx_mpz instead */
            std::vector<cxx_mpz> tr = mpz_poly_roots(f, cxx_mpz(q), qfac, rstate);
            /* convert back to type T */
            std::vector<T> res;
            res.reserve(tr.size());
            for(auto const & x : tr)
                res.push_back(get_from_cxx_mpz(x));
            return res;
        }
    }

    /* precompute the CRT: Qi*inverse mod qi of Qi, with Qi = q/qi.  */
    std::vector<T> Qi(qfac.size(), 1UL);
    std::vector<T> QQ(qfac.begin(), qfac.end());
    for(size_t s = 1 ; s < qfac.size() ; s <<= 1) {
        /* invariant:
         *
         * QQ[s*i] is the product of qfac[s*i...s*i+s-1]
         * Qi[s*i+k] is QQ[i] divided by qfac[s*i+k],
         *
         */
        for(size_t b = 0 ; b < qfac.size() ; b+=s) {
            if (b + s >= qfac.size()) break;
            for(size_t k = 0 ; k < s && b + k < qfac.size() ; k++)
                mul(Qi[b + k], Qi[b + k], QQ[b+s]);
            for(size_t k = 0 ; k < s && b + k + s < qfac.size() ; k++)
                mul(Qi[b + s + k], Qi[b + s + k], QQ[b]);
            mul(QQ[b], QQ[b], QQ[b+s]);
        }
    }
    std::vector<F> Ci(qfac.size(),1UL);
    if (qfac.size() >= 1) { /* otherwise a stupidly expensive no-op. */
        for(size_t i = 0 ; i < qfac.size() ; ++i) {
            /* Qi[i] has type T, which might be bigger than F.
             * poly_roots_impl_details<F>::invmod provides for this case
             */
            poly_roots_impl_details<F>::invmod(Ci[i], Qi[i], qfac[i]);
        }
    }

    /* compute the roots individually mod all prime factors. */
    std::vector<std::vector<F>> roots;
    for(size_t i = 0 ; i < qfac.size() ; ++i) {
        roots.push_back(mpz_poly_roots(f, qfac[i], rstate));
        if (roots.back().empty()) return res;
    }

    std::vector<T> sums;
    sums.push_back(0UL);
    for(size_t i = 0 ; i < qfac.size() ; ++i) {
        std::vector<T> newsums;
        for(auto r : roots[i]) {
            F rC;
            T rCQ;
            poly_roots_impl_details<F>::mul(rC, r, Ci[i]);
            poly_roots_impl_details<F>::mod(rC, rC, qfac[i]);
            mul(rCQ, Qi[i], rC);  /* rCQ is r mod qi and 0 mod all others */
            for(auto s : sums)
                newsums.push_back(add(s, rCQ));
        }
        std::swap(sums, newsums);
    }
    for(auto & x : sums)
        mod(x, x, q);
    std::sort(begin(sums), end(sums));
    res.insert(end(res), begin(sums), end(sums));
    return res;
}

template<typename T, typename F>
std::vector<T> mpz_poly_roots(cxx_mpz_poly const & f, T const & q, std::vector<F> const & qfac, gmp_randstate_ptr rstate)
{
    return poly_roots_impl<T,F>()(f, q, qfac, rstate);
}

template std::vector<cxx_mpz> mpz_poly_roots<cxx_mpz, cxx_mpz>(cxx_mpz_poly const & f, cxx_mpz const & q, std::vector<cxx_mpz> const & qfac, gmp_randstate_ptr rstate);
template std::vector<cxx_mpz> mpz_poly_roots<cxx_mpz, unsigned long>(cxx_mpz_poly const & f, cxx_mpz const & q, std::vector<unsigned long> const & qfac, gmp_randstate_ptr rstate);
template std::vector<unsigned long> mpz_poly_roots<unsigned long, unsigned long>(cxx_mpz_poly const & f, unsigned long const & q, std::vector<unsigned long> const & qfac, gmp_randstate_ptr rstate);
#ifndef UINT64_T_IS_EXACTLY_UNSIGNED_LONG
template std::vector<cxx_mpz> mpz_poly_roots<cxx_mpz, uint64_t>(cxx_mpz_poly const & f, cxx_mpz const & q, std::vector<uint64_t> const & qfac, gmp_randstate_ptr rstate);
template std::vector<uint64_t> mpz_poly_roots<uint64_t, uint64_t>(cxx_mpz_poly const & f, uint64_t const & q, std::vector<uint64_t> const & qfac, gmp_randstate_ptr rstate);
#endif

/* TODO: we'd like to expose the parallelized version via this interface
 * as well. However, the overloading approach show its limitations here.
 * The pinf handles would pollute the code all over the place...
 */
template<>
std::vector<cxx_mpz> mpz_poly_roots<cxx_mpz>(cxx_mpz_poly const & f, cxx_mpz const & q, gmp_randstate_ptr rstate)
{
    std::vector<cxx_mpz> tmp(f.degree());
    int n = mpz_poly_roots_mpz((mpz_t*) tmp.data(), f, q, rstate);
    tmp.erase(tmp.begin() + n, tmp.end());
    return tmp;
}
template<>
std::vector<unsigned long> mpz_poly_roots<unsigned long>(cxx_mpz_poly const & f, unsigned long const & q, gmp_randstate_ptr rstate)
{
    std::vector<unsigned long> tmp(f.degree());
    int n = mpz_poly_roots_ulong(tmp.data(), f, q, rstate);
    tmp.erase(tmp.begin() + n, tmp.end());
    return tmp;
}
#ifndef UINT64_T_IS_EXACTLY_UNSIGNED_LONG
template<>
std::vector<uint64_t> mpz_poly_roots<uint64_t>(cxx_mpz_poly const & f, uint64_t const & q, gmp_randstate_ptr rstate)
{
    std::vector<uint64_t> tmp(f.degree());
    int n = mpz_poly_roots_uint64(tmp.data(), f, q, rstate);
    tmp.erase(tmp.begin() + n, tmp.end());
    return tmp;
}
#endif

#if 0
int roots_for_composite_q(mpz_t* roots, mpz_poly_srcptr f,
        const mpz_t q, const unsigned long * fac_q)
{
    ASSERT_ALWAYS(fac_q[0] != 0);
    // only one prime left ?
    if (fac_q[1] == 0) {
        ASSERT(mpz_cmp_ui(q, fac_q[0]) == 0);
        return mpz_poly_roots(roots, f, q);
    }
    // First, a recursive call with q' = q / fac_q[0]
    mpz_t qp;
    mpz_init(qp);
    ASSERT(mpz_divisible_ui_p(q, fac_q[0]));
    mpz_divexact_ui(qp, q, fac_q[0]);
    int nr = roots_for_composite_q(roots, f, qp, fac_q+1);
    if (nr == 0) { // no roots modulo q'; we have finished.
        mpz_clear(qp);
        return 0;
    }

    // Second, compute the roots modulo fac_q[0]
    mpz_t fac_q0;
    mpz_init_set_ui(fac_q0, fac_q[0]);
    mpz_t roots2[MAX_DEGREE];
    for (int i = 0; i < MAX_DEGREE; ++i)
        mpz_init(roots2[i]);
    int nr2 = mpz_poly_roots(roots2, f, fac_q0);

    // Combine by CRT
    if (nr2 > 0) {
        mpz_t new_root, aux;
        mpz_init(new_root);
        mpz_init(aux);
        // pre-compute the coefficients of the CRT
        mpz_t c, c2;
        mpz_init(c);
        mpz_init(c2);
        int ret = mpz_invert(c, fac_q0, qp);
        ASSERT_ALWAYS(ret > 0);
        mpz_mul(c, c, fac_q0);
        ret = mpz_invert(c2, qp, fac_q0);
        ASSERT_ALWAYS(ret > 0);
        mpz_mul(c2, c2, qp);

        // reverse order to avoid erasing the input in roots[]
        for (int i = nr2-1; i >= 0; --i) {
            for (int j = 0; j < nr; ++j) {
                mpz_mul(new_root, roots[j], c);
                mpz_mul(aux, roots2[i], c2);
                mpz_add(new_root, new_root, aux);
                mpz_mod(new_root, new_root, q);
                mpz_set(roots[i*nr+j], new_root);
            }
        }
        mpz_clear(new_root);
        mpz_clear(aux);
    }

    for (int i = 0; i < MAX_DEGREE; ++i)
        mpz_clear(roots2[i]);
    mpz_clear(fac_q0);
    mpz_clear(qp);

    return nr*nr2; // can be 0.
}
#endif

/* Entry point for rootfind routines, for an integer n0 not necessarily prime.
   Since we cannot know in advance an easy bound on the number of
   roots, we allocate them in the function: if r = rp[0] at exit,
   the roots are r[0], r[1], ..., r[k-1] and the return value is k.
   Note: the elements r[j] must be mpz_clear'ed by the caller, and the
   array r also.

   TODO: this code is tested, but not used anywhere. I'm not sure that we want to keep it. Anyway it is wrong (roots of x^2-1 mod 4, for example). See test_rootfinder -v 4 1 0 -1

*/
unsigned long
mpz_poly_roots_gen (mpz_t **rp, mpz_poly_srcptr F, mpz_srcptr n, gmp_randstate_ptr rstate)
{
    ASSERT_ALWAYS (mpz_sgn (n) > 0);
    ASSERT_ALWAYS (F->deg >= 0);

    if (mpz_probab_prime_p (n, 1))
    {
        cxx_mpz p;
        mpz_set(p, n);
        std::vector<cxx_mpz> roots_p = mpz_poly_roots(F, p, rstate);
        *rp = (mpz_t *) malloc(roots_p.size() * sizeof(mpz_t));
        for(unsigned int i = 0 ; i < roots_p.size() ; i++)
            mpz_init_set((*rp)[i], roots_p[i]);
        return roots_p.size();
    }

    std::vector<cxx_mpz> results {{0}};

    /* now n is composite */

    /* invariants: Q = product of all factors we've dealt with so far, nn
     * = those that are still to be processed */
    cxx_mpz Q = 1, nn;
    mpz_set(nn, n);

    for (cxx_mpz p = 2 ; nn > 1 ; mpz_nextprime (p, p)) {
        if (mpz_probab_prime_p (nn, 1))
            mpz_set (p, nn);
        if (!mpz_divisible_p (nn, p))
            continue;

        std::vector<cxx_mpz> roots_p = mpz_poly_roots(F, p, rstate);
        
        /* Invariant: let nn0 = nn as it is here.
         * We'll maintain nn * q = nn0.
         */

        mpz_divexact (nn, nn, p);
        cxx_mpz q = p;

        for ( ;  mpz_divisible_p (nn, p) ; ) {
            mpz_mul (q, q, p);
            mpz_divexact (nn, nn, p);
            
            /* Lift each of the roots. We're doing it in a completely
             * stupid way for now. Of course we should compute the
             * derivative and so on, and special-case the ramified
             * situations.
             */
            std::vector<cxx_mpz> new_roots_p;

            for(auto const & r : roots_p) {
                cxx_mpz u = r;
                for(cxx_mpz x = 0 ; x < p ; x += 1, u += p) {
                    cxx_mpz v;
                    mpz_poly_eval (v, F, u);
                    if (mpz_divisible_p (v, q))
                        new_roots_p.push_back(u);
                }
            }
            roots_p = std::move(new_roots_p);
        }

        if (roots_p.empty()) {
            results.clear();
            break;
        }

        /* do a CRT between results[*] mod Q and roots_p[*] mod q */
        cxx_mpz xQ, xq;
        mpz_invert (xQ, Q, q); /* x = 1/Q mod q */
        mpz_invert (xq, q, Q); /* x = 1/q mod Q */

        cxx_mpz nQ = Q*q;

        std::vector<cxx_mpz> new_results;
        for(auto const & rQ : results) {
            cxx_mpz w0 = ((rQ * xq) % Q) * q;
            for(auto const & rq : roots_p) {
                cxx_mpz w = w0 + ((rq * xQ) % q) * Q;
                if (w >= nQ)
                    w -= nQ;
                new_results.push_back(w);
            }
        }
        mpz_swap(Q, nQ);
        results = std::move(new_results);
    }

    if (results.empty()) {
        *rp = NULL;
        return 0;
    }
    
    *rp = (mpz_t *) malloc(results.size() * sizeof(mpz_t));
    for(unsigned int i = 0 ; i < results.size() ; i++)
        mpz_init_set((*rp)[i], results[i]);
    return results.size();
}
