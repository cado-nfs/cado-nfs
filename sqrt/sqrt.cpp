/* square root, it can be used in two ways:

 * either do everything in one call:

   sqrt -poly cxxx.poly -prefix cxxx.dep.gz -purged cxxx.purged.gz -index cxxx.index.gz -ker cxxx.kernel

 * or in two steps:

   sqrt -poly cxxx.poly -prefix cxxx.dep.gz -purged cxxx.purged.gz -index cxxx.index.gz -ker cxxx.kernel -ab
   sqrt -poly cxxx.poly -prefix cxxx.dep.gz -sqrt0 -sqrt1 -gcd
 */

#include "cado.h" // IWYU pragma: keep

/* the following avoids the warnings "Unknown pragma" if OpenMP is not
   available, and should come after cado.h, which sets -Werror=all */
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif

#include <cerrno>
#include <cinttypes>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <fstream>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <utility> // pair
#include <vector>

#include <sys/stat.h>

#include <gmp.h>
#include "fmt/base.h"
#include "fmt/format.h"

#include "cado_poly.h"  // cado_poly
#include "cxx_mpz.hpp"   // for cxx_mpz
#include "filter_io.h"  // filter_rels
#include "gmp_aux.h"
#include "getprime.h"  // for getprime_mt, prime_info_clear, prime_info_init
#include "gzip.h"       // fopen_maybe_compressed
#include "memusage.h"   // PeakMemusage
#include "arith/modul_poly.h" // modul_poly
#include "mpz_poly.h"   // mpz_poly
#include "mpz_polymodF.h"   // mpz_poly
#include "mpz_poly_parallel.hpp"
#include "omp_proxy.h"
#include "purgedfile.h" // purgedfile_read_firstline
#include "version_info.h" // cado_revision_string
#include "portability.h" // strndup // IWYU pragma: keep
#include "macros.h"
#include "params.h"
#include "timing.h"

/* define to check the result of cxx_mpz_polymodF_sqrt */
// #define DEBUG

/* frequency of messages "read xxx (a,b) pairs" */
#define REPORT 10000000

/* maximal number of threads when reading dependency files */
#define MAX_IO_THREADS 16

// NOLINTBEGIN(cppcoreguidelines-avoid-non-const-global-variables)
static int verbose = 0;
static double wct0;

static std::mutex stdio_guard;
// NOLINTEND(cppcoreguidelines-avoid-non-const-global-variables)

static std::string
get_depname (std::string const & prefix, std::string const & algrat, unsigned int numdep)
{
    const std::string suffixes[] = {".gz", ".bz2", ".lzma"};
    std::string suffix;

    for (auto const & t : suffixes) {
        if (prefix.substr(prefix.size() - t.size()) == t) {
            suffix = t;
            break;
        }
    }
    auto prefix_base = prefix.substr(0, prefix.size() - suffix.size());
    return fmt::format("{}.{}{:03d}{}", prefix_base, algrat, numdep, suffix);
}

static std::string
get_depsidename (std::string const & prefix, unsigned int numdep, int side)
{
  return get_depname(prefix, fmt::format("side{}.", side), numdep);
}

static FILE*
fopen_maybe_compressed_lock (const char * name, const char * mode)
{
  const std::lock_guard<std::mutex> dummy(stdio_guard);
  return fopen_maybe_compressed (name, mode);
}

static int
fclose_maybe_compressed_lock (FILE * f, const char * name)
{
  const std::lock_guard<std::mutex> dummy(stdio_guard);
  return fclose_maybe_compressed (f, name);
}

/* this function is run sequentially, thus no need to be thread-safe */
static int
check_dep (const char *prefix, unsigned int numdep)
{
  const std::string depname = get_depname (prefix, "", numdep);
  const std::ifstream f(depname);
  return bool(f);
}

/*************************************************************************/
/* Function to handle the initial accumulation (rat & alg) */


/* replace the vector of elements of type T by
 * the product of elements in the range. It's fairly trivial to do this
 * in a DFS manner, but it's slightly less so for a BFS algorithm. The
 * latter is more amenable to parallelization.
 *
 * The most balanced split is of course when the size of the vector v is a
 * power of two. When the size of the vector v is N=2^k+r, the optimal way
 * to fall back to the balanced case is as follows. Form a vector w of
 * length 2^k by moving all elements of the vector v to the vector w, one
 * by one, except that the i-th element of w is created from *two*
 * consecutive elements of v whenever the bit-reversal of i is less than
 * r.
 *
 * Of course we prefer to create w from v in place, which is done easily
 * by having two pointers go through v.
 */
static inline uint64_t bitrev(uint64_t a)
{
    a = (a >> 32U) ^ (a << 32U);
    uint64_t m;
    m = UINT64_C(0x0000ffff0000ffff);
    a = ((a >> 16U) & m) ^ ((a << 16U) & ~m);
    m = UINT64_C(0x00ff00ff00ff00ff);
    a = ((a >> 8U) & m) ^ ((a << 8U) & ~m);
    m = UINT64_C(0x0f0f0f0f0f0f0f0f);
    a = ((a >> 4U) & m) ^ ((a << 4U) & ~m);
    m = UINT64_C(0x3333333333333333);
    a = ((a >> 2U) & m) ^ ((a << 2U) & ~m);
    m = UINT64_C(0x5555555555555555);
    a = ((a >> 1U) & m) ^ ((a << 1U) & ~m);
    return a;
}

/* Modify the range [vb, ve[ so that in the end, *vb contains the product
 * of the range, and [vb+1, ve[ contains only ones.
 *
 * The array A is used internally by this function, but we provide a
 * prototype that exposes it because two consecutive calls to this
 * function with similar-size ranges can profitably re-use the storage in
 * A.
 */
template<typename M>
static void accumulate(std::vector<typename M::T> & A, typename std::vector<typename M::T>::iterator vb, typename std::vector<typename M::T>::iterator ve, M const & m)
{
    /* This is a single-threaded routine. We put the result in *vb, and
     * the caller is reponsible of freeing the range [vb+1,ve[
     */
    if (ve - vb < 16) {
        /* Don't bother with very small vectors. Here we have a vector of
         * size less than 16*2*thr/2, so don't bother. */
        for(auto vi = vb+1; vi < ve ; ++vi) {
            m(*vb, *vb, *vi);
        }
        return;
    }
    constexpr const unsigned int ratio = 2; /* must be an integer > 1 */
    /* create an array of products of exponentially growing size.
     *
     * after each input value considered, we maintain the property that
     * the i-th cell in the vector A[] contains the product of at most
     * ratio^i items in the original range [vb, ve[
     *
     * Hence with n cells int the table A, we can have up to
     * (ratio^n-1)/(ratio-1) values stored.
     */
    const size_t as = 1 + (ratio - 1) * (ve - vb);
    unsigned int ncells = 0;
    for(size_t s = 1 ; s < as ; s *= ratio, ncells++) ;
    /* make sure A is large enough */
    for(; A.size() < ncells; A.emplace_back());
    for(auto& x: A) m.set1(x);

    for(auto vi = vb ; vi < ve ; ++vi) {
        m(A[0], A[0], *vi);
        /* We won't need *vi again. This frees indirect storage */
        *vi = typename M::T();
        /* now maintain the condition */
        unsigned int nprd = vi + 1 - vb;
        for(unsigned int i = 0 ; nprd % ratio == 0 ; ++i, nprd /= ratio) {
            ASSERT_ALWAYS(i + 1 < A.size());
            m(A[i+1], A[i+1], A[i]);
            /* Do not replace by a fresh object. Prefer something that
             * does not trigger a free() yet, since reallocation would
             * occur eventually */
            m.set1(A[i]);
        }
    }
    /* final: accumulate everything in A.back(). Note that by doing it
     * this way, we should normally not trigger reallocation. */
    for(unsigned int i = 0 ; i + 1 < A.size() ; i++) {
        if (m.is1(A[i+1])) {
            // it is tempting to *swap* here, but in fact it would be a
            // terrible idea, since we're being careful to keep the
            // *storage* of A[] grow exponentially. With swap's, we would
            // in effect make all cells stabilize at the maximal storage,
            // which would incur a significant increase of the total
            // storage amount!
            // std::swap(A[i], A[i+1]);
            m.set(A[i+1], A[i]);
            m.set1(A[i]);
        } else {
            m(A[i+1], A[i+1], A[i]);
            // do not free A[i] yet, as we might want to reuse it for another
            // block. Set it to 1 instead.
            // A[i]=T();
            m.set1(A[i]);
        }
    }
    *vb = std::move(A.back());
}
/* not used, just to mention that the simplest prototype for the function
 * above would be as follows, in a way */
template<typename M>
static void accumulate(typename std::vector<typename M::T>::iterator vb, typename std::vector<typename M::T>::iterator ve, M const & m)
{
    std::vector<typename M::T> A;
    accumulate(A, vb, ve, m);
}
/* This does the "floor" level. From a vector of T's or arbitrary size,
 * return a modified vector whose size is a power of two, and whose
 * elements are products of the original elements in sub-ranges of the
 * initial vector */
template<typename M>
static void accumulate_level00(std::vector<typename M::T> & v, M const & m, std::string const & message)
{
    unsigned int nthr;
    /* We want to know the number of threads that we get when starting a
     * parallel region -- but this is probably unneeded. I think that
     * omp_get_max_threads would do an equally good job. */
    #pragma omp parallel
    #pragma omp critical
    nthr = omp_get_num_threads (); /* total number of threads */

    size_t N = v.size();
    /* We're going to split the input in that 1<<n pieces exacly
     * -- we want a power of two, that is a strict condition. We'll
     *  strive to make pieces of balanced size.
     */
    /* min_pieces per_thread is there so that each thread at level 0 has
     * at least that many pieces to work on. Increasing this parameter
     * allows more progress reporting, but more importantly this acts on
     * the balance of the computation. (having "2 or 3" pieces to handle
     * per threads leads to longer idle wait than if we have "16 or 17"
     * pieces to handle, of course)
     */
    constexpr const unsigned int min_pieces_per_thread = 8;
    unsigned int n = 1;
    for( ; (1UL << n) < min_pieces_per_thread * nthr ; n++) ;

    if ((N >> n) < 16) {
        /* Don't bother with very small vectors. Here we have a vector of
         * size less than 16*min_pieces_per_thread*thr/2, so don't bother. */
        for(size_t j = 1 ; j < v.size() ; j++) {
            /* use the instance with no multithreading. This is on
             * purpose, for small ranges and small objects.
             */
            m(v[0], v[0], v[j]);
        }
        v.erase(v.begin() + 1, v.end());
        return;
    }

    const uint64_t r = N % (1UL << n);
    {
        fmt::print (stderr, "{}: starting level 00 at wct={:1.2f}s,"
                " {} -> 2^{}*{}+{}\n",
                message, wct_seconds () - wct0,
                v.size(), n, N>>n, r);
        fflush (stderr);
    }

    /* Interesting job: we want to compute the start and end points of
     * the i-th part among 1<<n of a vector of size N.
     *
     * we know that this i-th part has size (N>>n)+1 if and only if the
     * n-bit reversal of i is less than N mod (1<<n) (denoted by r)
     */
    std::vector<size_t> endpoints;
    endpoints.push_back(0);
    for(uint64_t i = 0 ; i < (UINT64_C(1) << n) ; i++) {
        const uint64_t ir = bitrev(i) >> (64U - n);
        const size_t fragment_size = (N>>n) + (ir < (uint64_t) r);
        endpoints.push_back(endpoints.back() + fragment_size);
    }
    ASSERT_ALWAYS(endpoints.back() == v.size());

    /* Now do the real job */
#pragma omp parallel
    {
        /* This array is thread-private, and we're going to re-use it for
         * all sub-ranges that we multiply. The goal at this point is to
         * avoid hammering the malloc layer.
         */
        std::vector<typename M::T> A;
#pragma omp for
        for(unsigned int i = 0 ; i < (1U << n) ; i++) {
            auto vb = v.begin() + endpoints[i];
            auto ve = v.begin() + endpoints[i+1];
            accumulate(A, vb, ve, m);
            {
                const std::lock_guard<std::mutex> dummy(stdio_guard);
                fmt::print (stderr, "{}: fragment {}/{}"
                        " of level 00 done by thread {} at wct={:1.2f}s\n",
                        message,
                        i, 1U<<n, omp_get_thread_num(),
                        wct_seconds () - wct0);
                fflush (stderr);
            }
        }
    }
    /* This puts pressure at the malloc level, but takes negligible time */
    {
        auto dst = v.begin();
        endpoints.pop_back();
        for(const size_t z : endpoints)
            std::swap(v[z], *dst++);
        v.erase(dst, v.end());
    }

    N = v.size();

    ASSERT_ALWAYS(!(N & (N-1)));
}


template<typename M>
static typename M::T accumulate(std::vector<typename M::T> & v, M const & m, std::string const & message)
{
    accumulate_level00(v, m, message);

    unsigned int nthr;
    /* We want to know the number of threads that we get when starting a
     * parallel region -- but this is probably unneeded. I think that
     * omp_get_max_threads would do an equally good job. */
    #pragma omp parallel
    #pragma omp critical
    nthr = omp_get_num_threads (); /* total number of threads */

    const size_t vs = v.size();
    ASSERT_ALWAYS(!(vs & (vs - 1)));

    /* At this point v has size a power of two */
  for(int level = 0 ; v.size() > 1 ; level++) {
      {
        const std::lock_guard<std::mutex> dummy(stdio_guard);
	fmt::print (stderr, "{}: starting level {} at cpu={:1.2f}s (wct={:1.2f}s), {} values to multiply\n",
		      message, level, seconds (), wct_seconds () - wct0, v.size());
	fflush (stderr);
      }

      const double st = seconds (), wct = wct_seconds ();
      const size_t N = v.size() - 1;
      int local_nthreads;
      /* the loop below performs floor((N+1)/2) products */
      const size_t nloops = (N + 1) / 2;
      if (nthr < 2 * nloops)
	{
	  omp_set_nested (0);
	  local_nthreads = 1;
	}
      else
	{
	  /* we have to set omp_set_nested here and not at the beginning of
	     main(), since it seems that the pthreads reset OpenMP's "nested"
	     value to 0 */
	  omp_set_nested (1);
	  local_nthreads = nthr / nloops;
	}
#pragma omp parallel for
      for(size_t j = 0 ; j < N ; j += 2) {
	  omp_set_num_threads (local_nthreads);
	  m(v[j], v[j], v[j+1]);
          v[j+1] = typename M::T();
      }

      /* reset "nested" to 0 */
      omp_set_nested (0);

      /* shrink (not parallel), takes negligible time */
      for(size_t j = 2 ; j < v.size() ; j += 2) {
          std::swap(v[j], v[j/2]);
      }
      v.erase(v.begin() + (v.size() + 1) / 2, v.end());
      {
        const std::lock_guard<std::mutex> dummy(stdio_guard);
	fmt::print (stderr, "{}: level {} took cpu={:1.2f}s (wct={:1.2f}s)\n",
		      message, level, seconds () - st, wct_seconds () - wct);
	fflush (stderr);
      }
  }
  return std::move(v.front());
}

/*************************************************************************/
/* Parallel I/O for reading the dep file */

template<typename M>
static std::vector<typename M::T>
read_ab_pairs_from_depfile(std::string const & depname, M const & m, std::string const & message, unsigned long & nab, unsigned long & nfree)
{
    nab = nfree = 0;
    FILE * depfile = fopen_maybe_compressed_lock (depname.c_str(), "rb");
    ASSERT_ALWAYS(depfile);
    std::vector<typename M::T> prd;

    if (fseek(depfile, 0, SEEK_END) < 0) {
        fmt::print(stderr, "{}: cannot seek in dependency file, using single-thread I/O\n", message);
        cxx_mpz a, b;
        /* Cannot seek: we have to use serial i/o */
        while (gmp_fscanf(depfile, "%Zd %Zd", (mpz_ptr) a, (mpz_ptr) b) != EOF)
        {
            if(!(nab % REPORT)) {
                fmt::print(stderr, "{}: read {} (a,b) pairs in {:.2f}s (wct {:.2f}s, peak {}M)\n",
                        message, nab, seconds (), wct_seconds () - wct0,
                        PeakMemusage () >> 10U);
                fflush (stderr);
            }
            if (mpz_cmp_ui (a, 0) == 0 && mpz_cmp_ui (b, 0) == 0)
                break;
            prd.emplace_back(m.from_ab(a, b));
            nab++;
            if (mpz_cmp_ui (b, 0) == 0)
                nfree++;
        }
    } else {
        /* We _can_ seek. Good news ! */
        const long endpos = ftell(depfile);
        /* Find accurate starting positions for everyone */
        unsigned int nthreads = omp_get_max_threads();
        /* cap the number of I/O threads */
        if (nthreads > MAX_IO_THREADS)
            nthreads = MAX_IO_THREADS;
        fmt::print(stderr, "{}: Doing I/O with {} threads\n", message, nthreads);
        std::vector<long> spos_tab;
        spos_tab.reserve(nthreads);
        for(unsigned int i = 0 ; i < nthreads ; i++)
            spos_tab.push_back((endpos * i) / nthreads);
        spos_tab.push_back(endpos);
        /* All threads get their private reading head. */
#pragma omp parallel num_threads(nthreads)
        {
            const int i = omp_get_thread_num();
            FILE * fi = fopen_maybe_compressed_lock (depname.c_str(), "rb");
            const int rc = fseek(fi, spos_tab[i], SEEK_SET);
            ASSERT_ALWAYS(rc == 0);
            if (i > 0) {
                /* Except when we're at the end of the stream, read until
                 * we get a newline */
                for( ; fgetc(fi) != '\n' ; ) ;
            }
            spos_tab[i] = ftell(fi);
#pragma omp barrier
            std::vector<typename M::T> loc_prd;
            unsigned long loc_nab = 0;
            unsigned long loc_nfree = 0;
            cxx_mpz a, b;
            for(long pos = ftell(fi) ; pos < spos_tab[i+1] ; ) {
                const int rc = gmp_fscanf(fi, "%Zd %Zd", (mpz_ptr) a, (mpz_ptr) b);
                pos = ftell(fi);
                /* We may be tricked by initial white space, since our
                 * line parsing does not gobble whitespace at end of
                 * line.
                 * Therefore we must check the position _after_ the
                 * parse.
                 */
                if ((rc != 2 && feof(fi)) || pos >= spos_tab[i+1])
                    break;
                ASSERT_ALWAYS(rc == 2);
                if (mpz_cmp_ui (a, 0) == 0 && mpz_cmp_ui (b, 0) == 0)
                    ASSERT_ALWAYS(0);
                loc_prd.emplace_back(m.from_ab(a, b));
                loc_nab++;
                loc_nfree += (mpz_cmp_ui (b, 0) == 0);

                if(!(loc_nab % 100000))
#pragma omp critical
                {
                    nab += loc_nab;
                    nfree += loc_nfree;
                    loc_nab = 0;
                    loc_nfree = 0;
                    /* in truth, a splice() would be best. */
                    for(auto & x: loc_prd)
                        prd.emplace_back(std::move(x));
                    loc_prd.clear();
                    if(!(nab % REPORT)) {
                        fmt::print(stderr, "{}: read {} (a,b) pairs in {:.2f}s (wct {:.2f}s, peak {}M)\n",
                                message, nab, seconds (), wct_seconds () - wct0,
                                PeakMemusage () >> 10U);
                        fflush (stderr);
                    }
                }
            }
#pragma omp critical
            {
                nab += loc_nab;
                nfree += loc_nfree;
                for(auto & x: loc_prd)
                    prd.emplace_back(std::move(x));
                loc_prd.clear();
            }
            fclose_maybe_compressed_lock (fi, depname.c_str());
        }
        {
            fmt::print (stderr, "{} read {} (a,b) pairs, including {} free, in {:1.2f}s (wct {:1.2f}s)\n",
                    message, nab, nfree, seconds (), wct_seconds () - wct0);
            fflush (stderr);
        }
    }
    fclose_maybe_compressed_lock (depfile, depname.c_str());
    return prd;
}

/*************************************************************************/


/********** RATSQRT **********/

/* This is where we store side-relative functions for the rational square
 * root. We deal with integers, here.
 */
struct cxx_mpz_functions {
    cxx_mpz_poly const & P;
    using T = cxx_mpz;
    static void set1(T & x) { mpz_set_ui(x, 1); }
    static void set(T & y, T const & x) { mpz_set(y, x); }
    static bool is1(T & x) { return mpz_cmp_ui(x, 1) == 0; }
    void operator()(T & res, T const & a, T const & b) const {
        mpz_mul(res, a, b);
    }
    ATTRIBUTE_NODISCARD
    T from_ab(cxx_mpz const& a, cxx_mpz const& b) const
    {
        cxx_mpz v;
        /* accumulate g1*a+g0*b */
        mpz_mul (v, mpz_poly_coeff_const(P, 1), a);
        mpz_addmul (v, mpz_poly_coeff_const(P, 0), b);
        return v;
    }
    explicit cxx_mpz_functions(cxx_mpz_poly const & P) : P(P) {}
};


static int
calculateSqrtRat (std::string const & prefix, unsigned int numdep, cxx_cado_poly const & cpoly,
        int side, cxx_mpz const & Np)
{
  const std::string sidename = get_depsidename (prefix, numdep, side);
  FILE *resfile;
  unsigned long nab = 0, nfree = 0;

  ASSERT_ALWAYS (cpoly->pols[side]->deg == 1);

#pragma omp critical
  {
#ifdef __MPIR_VERSION
    fmt::print (stderr, "Using MPIR {}\n", mpir_version);
#else
    fmt::print (stderr, "Using GMP {}\n", gmp_version);
#endif
    fflush (stderr);
  }

  cxx_mpz_poly F(cpoly->pols[side]);
  cxx_mpz prod;

  {
      const cxx_mpz_functions M(F);
      const std::string message = fmt::format("Rat({})", numdep);

      std::vector<cxx_mpz> prd = read_ab_pairs_from_depfile(
              get_depname (prefix, "", numdep),
              M,
              message,
              nab,
              nfree);

      prod = accumulate(prd, M, message);
  }

  /* we must divide by g1^nab: if the number of (a,b) pairs is odd, we
     multiply by g1, and divide by g1^(nab+1) */
  if (nab & 1U)
    mpz_mul (prod, prod, mpz_poly_coeff_const(F, 1));

#pragma omp critical
  {
    fmt::print (stderr, "Rat({}): size of product = {} bits (peak {}M)\n",
	     numdep, mpz_sizeinbase (prod, 2),
	     PeakMemusage () >> 10U);
    fflush (stderr);
  }

  if (mpz_sgn (prod) < 0)
    {
      fmt::print (stderr, "Error, product is negative: try another dependency\n");
      exit(EXIT_FAILURE);
    }

#pragma omp critical
  {
    fmt::print (stderr, "Rat({}): starting rational square root at {:.2f}s (wct {:.2f}s)\n",
	     numdep, seconds (), wct_seconds () - wct0);
    fflush (stderr);
  }

  cxx_mpz v;
  /* since we know we have a square, take the square root */
  mpz_sqrtrem (prod, v, prod);

#pragma omp critical
  {
    fmt::print (stderr, "Rat({}): computed square root at {:.2f}s (wct {:.2f}s)\n",
	     numdep, seconds (), wct_seconds () - wct0);
    fflush (stderr);
  }

  if (mpz_cmp_ui (v, 0) != 0)
    {
      unsigned long p = 2, e, errors = 0;
      cxx_mpz pp;

      fmt::print (stderr, "Error, rational square root remainder is not zero\n");
      /* reconstruct the initial value of prod to debug */
      mpz_mul (prod, prod, prod);
      mpz_add (prod, prod, v);
      prime_info pi;
      prime_info_init (pi);
      while (mpz_cmp_ui (prod, 1) > 0)
        {
          e = 0;
          if (verbose)
            fmt::print ("Removing p={}:", p);
          mpz_set_ui (pp, p);
          e = mpz_remove (prod, prod, pp);
          if (verbose)
            fmt::print (" exponent={}, remaining {} bits\n", e,
                    mpz_sizeinbase (prod, 2));
          if ((e % 2) != 0)
            {
              errors ++;
              fmt::print (stderr, "Prime {} appears to odd power {}\n", p, e);
              if (verbose || errors >= 10)
                break;
            }
          p = getprime_mt (pi);
        }
      prime_info_clear (pi);
      exit(EXIT_FAILURE);
    }

  mpz_mod (prod, prod, Np);

#pragma omp critical
  {
    fmt::print (stderr, "Rat({}): reduced mod n at {:.2f}s (wct {:.2f}s)\n",
	     numdep, seconds (), wct_seconds () - wct0);
    fflush (stderr);
  }

  /* now divide by g1^(nab/2) if nab is even, and g1^((nab+1)/2)
     if nab is odd */

  mpz_powm_ui (v, mpz_poly_coeff_const(F, 1), (nab + 1) / 2, Np);
#pragma omp critical
  {
    fmt::print (stderr, "Rat({}): computed g1^(nab/2) mod n at {:.2f}s (wct {:.2f}s)\n",
	     numdep, seconds (), wct_seconds () - wct0);
    fflush (stderr);
  }
  resfile = fopen_maybe_compressed_lock (sidename.c_str(), "wb");

  mpz_invert (v, v, Np);
  mpz_mul (prod, prod, v);
  mpz_mod (prod, prod, Np);

  fmt::print (resfile, "{}\n", prod);
  fclose_maybe_compressed_lock (resfile, sidename.c_str());

#pragma omp critical
  {
      fmt::print (stderr, "Rat({}): square root is {}\n", numdep, prod);
    fmt::print (stderr, "Rat({}): square root time: {:.2f}s (wct {:.2f}s)\n",
	     numdep, seconds (), wct_seconds () - wct0);
    fflush (stderr);
  }

  return 0;
}

/********** ALGSQRT **********/

/* This is the analogue of cxx_mpz_functions for the algebraic square
 * root. Here, the object must carry a reference to the field polynomial
 */
struct cxx_mpz_polymodF_functions {
    using T = cxx_mpz_polymodF;
    cxx_mpz_poly const & F;
    const mpz_poly_parallel_info * pinf;
    explicit cxx_mpz_polymodF_functions(cxx_mpz_poly & F, const mpz_poly_parallel_info * pinf = nullptr)
        : F(F)
        , pinf(pinf)
    {}
    
    ATTRIBUTE_NODISCARD
    static T from_ab(cxx_mpz const & a, cxx_mpz const & b) {
        cxx_mpz_polymodF ret;
        mpz_polymodF_set_from_ab(ret, a, b);
        return ret;
    }
    static bool is1(T & res) {
        return res->p->deg == 0 && res->v == 0 && mpz_cmp_ui(mpz_poly_coeff_const(res->p, 0), 1) == 0;
    }
    static void set(T & y, T const & x) {
        mpz_polymodF_set(y, x);
    }
    static void set1(T & res) {
        mpz_polymodF_set_ui(res, 1);
    }
    void operator()(T &res, T const & a, T const & b) const {
        pinf->mpz_polymodF_mul(res, a, b, F);
    }
};


/* Reduce the coefficients of R in [-m/2, m/2) */
static void
mpz_poly_mod_center (mpz_poly R, const mpz_t m)
{
#pragma omp parallel for
  for (int i=0; i <= R->deg; i++)
    mpz_ndiv_r (mpz_poly_coeff(R, i), mpz_poly_coeff_const(R, i), m);
}

#if 0
/* Check whether the coefficients of R (that are given modulo m) are in
   fact genuine integers. We assume that x mod m is a genuine integer if
   x or |x-m| is less than m/10^6, i.e., the bit size of x or |x-m| is
   less than that of m minus 20.
   Assumes the coefficients x satisfy 0 <= x < m.
*/
static int
mpz_poly_integer_reconstruction (mpz_poly R, const mpz_t m)
{
  int i;
  size_t sizem = mpz_sizeinbase (m, 2), sizer;

  for (i=0; i <= R->deg; ++i)
    {
      sizer = mpz_sizeinbase (mpz_poly_coeff_const(R, i), 2);
      if (sizer + 20 > sizem)
        {
          mpz_sub (mpz_poly_coeff(R, i), mpz_poly_coeff_const(R, i), m);
          sizer = mpz_sizeinbase (mpz_poly_coeff_const(R, i), 2);
          if (sizer + 20 > sizem)
            return 0;
        }
    }
  return 1;
}
#endif

// compute res := sqrt(a) in Fp[x]/f(x)
static void
TonelliShanks (mpz_poly res, const mpz_poly a, const mpz_poly F, unsigned long p)
{
  const int d = F->deg;
  cxx_mpz q;
  cxx_mpz_poly delta;  // a non quadratic residue
  cxx_mpz_poly auxpol;
  cxx_mpz aux;
  cxx_mpz t;
  unsigned int s;
  cxx_mpz myp = p;

  mpz_ui_pow_ui(q, p, (unsigned long)d);

  // compute aux = (q-1)/2
  // and (s,t) s.t.  q-1 = 2^s*t
  mpz_sub_ui(aux, q, 1);
  mpz_divexact_ui(aux, aux, 2);

  for(t = aux, s = 1 ; mpz_divisible_2exp_p(t, 1) ; s++)
      mpz_fdiv_q_2exp(t, t, 1);

  // find a non quadratic residue delta

  {
    cxx_gmp_randstate state;
    do {
      int i;
      // pick a random delta
      for (i = 0; i < d; ++i)
	mpz_urandomm(mpz_poly_coeff(delta, i), state, myp);
      mpz_poly_cleandeg(delta, d-1);
      // raise it to power (q-1)/2
      mpz_poly_pow_mod_f_mod_ui(auxpol, delta, F, aux, p);
      /* Warning: the coefficients of auxpol might either be reduced in
	 [0, p) or in [-p/2, p/2). This code should work in both cases. */
    } while (auxpol->deg != 0 || (mpz_cmp_ui (mpz_poly_coeff_const(auxpol, 0), p-1) != 0 &&
				  mpz_cmp_si (mpz_poly_coeff_const(auxpol, 0), -1) != 0));
  }

  // follow the description of Crandall-Pomerance, page 94
  {
    cxx_mpz_poly A, D;
    cxx_mpz m = 0;
    mpz_poly_pow_mod_f_mod_ui(A, a, F, t, p);
    mpz_poly_pow_mod_f_mod_ui(D, delta, F, t, p);
    for (unsigned int i = 0; i < s; ++i) {
        mpz_poly_pow_mod_f_mod_ui(auxpol, D, F, m, p);
        mpz_poly_mul_mod_f_mod_mpz(auxpol, auxpol, A, F, myp, nullptr, nullptr);
        mpz_ui_pow_ui(aux, 2, (s-1-i));
        mpz_poly_pow_mod_f_mod_ui(auxpol, auxpol, F, aux, p);
        if ((auxpol->deg == 0) && (mpz_cmp_ui(mpz_poly_coeff_const(auxpol, 0), p-1)== 0))
            mpz_add_ui(m, m, 1UL<<i);
    }
    mpz_add_ui(t, t, 1);
    mpz_divexact_ui(t, t, 2U);
    mpz_poly_pow_mod_f_mod_ui(res, a, F, t, p);
    mpz_divexact_ui(m, m, 2U);
    mpz_poly_pow_mod_f_mod_ui(auxpol, D, F, m, p);

    mpz_poly_mul_mod_f_mod_mpz(res, res, auxpol, F, myp, nullptr, nullptr);
  }
}

// res <- Sqrt(AA) mod F, using p-adic lifting, at prime p.
static unsigned long
cxx_mpz_polymodF_sqrt (cxx_mpz_polymodF & res, cxx_mpz_polymodF & AA, cxx_mpz_poly const & F, unsigned long p,
	       unsigned int numdep,
               const mpz_poly_parallel_info * pinf)
{
  mpz_poly A, *P;
  unsigned long v;
  const int d = F->deg;
  unsigned long k, target_k;
  unsigned long K[65];
  int lk, logk, logk0;
  size_t target_size; /* target bit size for Hensel lifting */

  /* The size of the coefficients of the square root of A should be about half
     the size of the coefficients of A. Here is an heuristic argument: let
     K = Q[x]/(f(x)) where f(x) is the algebraic polynomial. The square root
     r(x) might be considered as a random element of K: it is smooth, not far
     from an integer, but except that has no relationship with the coefficients
     of f(x). When we square r(x), we obtain a polynomial with coefficients
     twice as large, before reduction by f(x). The reduction modulo f(x)
     produces A(x), however that reduction should not decrease the size of
     the coefficients. */
  target_size = mpz_poly_sizeinbase (AA->p, 2);
  target_size = target_size / 2;

  /* note that we scale by the derivative of f_hat, which might increase
   * the coefficient a little bit. I suppose that the bitsize of the
   * discriminant is a safe bet. Better take into account */
  {
      cxx_mpz disc;
      mpz_poly_discriminant(disc, F);
      target_size += mpz_sizeinbase(disc, 2);
  }

  target_size += target_size / 10;
#pragma omp critical
  {
    fmt::print (stderr, "Alg({}): target_size={}\n", numdep,
	     (unsigned long int) target_size);
    fflush (stderr);
  }

  mpz_poly_init(A, d-1);
  // Clean up the mess with denominator: if it is an odd power of fd,
  // then multiply num and denom by fd to make it even.
  mpz_poly_swap(A, AA->p);
  if (((AA->v) & 1U) == 0) {
    v = AA->v / 2;
  } else {
    v = (1+AA->v) / 2;
    pinf->mpz_poly_mul_mpz(A, A, mpz_poly_coeff_const(F, d));
  }

  // Now, we just have to take the square root of A (without denom) and
  // divide by fd^v.

  // Variables for the lifted values
  mpz_poly invsqrtA;
  // variables for A and F modulo pk
  mpz_poly a;
  mpz_poly_init(invsqrtA, d-1);
  mpz_poly_init(a, d-1);
  // variable for the current pk
  mpz_t pk;
  mpz_init (pk);

  /* Jason Papadopoulos's trick: since we will lift the square root of A to at
     most target_size bits, we can reduce A accordingly */
  double st = seconds (), wct = wct_seconds ();
  target_k = (unsigned long) ((double) target_size * log ((double) 2) / log((double) p));
  mpz_ui_pow_ui (pk, p, target_k);
  while (mpz_sizeinbase (pk, 2) <= target_size)
    {
      mpz_mul_ui (pk, pk, p);
      target_k ++;
    }
  pinf->mpz_poly_mod_mpz (A, A, pk, nullptr);
  for (k = target_k, logk = 0; k > 1; k = (k + 1) / 2, logk ++)
    K[logk] = k;
  K[logk] = 1;
#pragma omp critical
  {
    fmt::print (stderr, "Alg({}): reducing A mod p^{} took {:.2f}s (wct {:.2f}s)\n",
	     numdep, target_k, seconds () - st, wct_seconds () - wct);
    fflush (stderr);
  }

  // Initialize things modulo p:
  mpz_set_ui (pk, p);
  k = 1; /* invariant: pk = p^k */
  lk = 0; /* k = 2^lk */
  st = seconds ();
  wct = wct_seconds ();
  P = pinf->mpz_poly_base_modp_init (A, p, K, logk0 = logk);
#pragma omp critical
  {
    fmt::print (stderr, "Alg({}): mpz_poly_base_modp_init took {:.2f}s (wct {:.2f}s)\n",
	     numdep, seconds () - st, wct_seconds () - wct);
    fflush (stderr);
  }

  /* A is not needed anymore, thus we can clear it now */
  mpz_poly_clear (A);

  mpz_poly_set (a, P[0]);

  // First compute the inverse square root modulo p
  {
    mpz_t q, aux;
    mpz_init(q);
    mpz_init(aux);
    mpz_ui_pow_ui(q, p, (unsigned long)d);

#if 0
    // compute (q-2)(q+1)/4   (assume q == 3 mod 4, here !!!!!)
    // where q = p^d, the size of the finite field extension.
    // since we work mod q-1, this gives (3*q-5)/4
    mpz_mul_ui(aux, q, 3);
    mpz_sub_ui(aux, aux, 5);
    mpz_divexact_ui(aux, aux, 4);               // aux := (3q-5)/4
    mpz_poly_pow_mod_f_mod_ui(invsqrtA, a, F, aux, p);
#else
    TonelliShanks(invsqrtA, a, F, p);
    mpz_sub_ui(aux, q, 2);
    mpz_poly_pow_mod_f_mod_ui(invsqrtA, invsqrtA, F, aux, p);
#endif

    mpz_clear(aux);
    mpz_clear(q);
  }

  // Now, the lift begins
  // When entering the loop, invsqrtA contains the inverse square root
  // of A computed modulo p.

  mpz_poly tmp;
  mpz_t invpk;
  mpz_init (invpk);
  mpz_poly_init (tmp, 2*d-1);
  do {
    double st;

    if (mpz_sizeinbase (pk, 2) > target_size)
      {
        fmt::print (stderr, "Failed to reconstruct an integer polynomial\n");
        fmt::print ("Failed\n");
        exit(EXIT_FAILURE);
      }

    /* invariant: invsqrtA = 1/sqrt(A) bmod p^k */

    lk += 1;
    st = seconds ();
    wct = wct_seconds ();
    /* a <- a + pk*P[lk] */
    pinf->mpz_poly_base_modp_lift (a, P, lk, pk);
    /* free P[lk] which is no longer needed */
    mpz_poly_clear (P[lk]);
    if (verbose)
#pragma omp critical
      {
	fmt::print (stderr, "Alg({}):    mpz_poly_base_modp_lift took {:.2f}s (wct {:.2f}s, peak {}M)\n",
		 numdep, seconds () - st, wct_seconds () - wct,
		 PeakMemusage () >> 10U);
	fflush (stderr);
      }

    mpz_mul (pk, pk, pk);   // double the current precision
    k = k + k;
    logk --;
    if (K[logk] & 1U)
      {
        mpz_div_ui (pk, pk, p);
        k --;
      }
    barrett_precompute_inverse (invpk, pk);

    /* check the invariant k = K[logk] */
    ASSERT_ALWAYS(k == K[logk]);

#pragma omp critical
    {
      fmt::print (stderr, "Alg({}): start lifting mod p^{} ({} bits) at {:.2f}s (wct {:.2f}s)\n",
	       numdep, k, (unsigned long int) mpz_sizeinbase (pk, 2),
	       seconds (), wct_seconds () - wct0);
      fflush (stderr);
    }

    // now, do the Newton operation x <- 1/2(3*x-a*x^3)
    st = seconds ();
    wct = wct_seconds ();
    pinf->mpz_poly_sqr_mod_f_mod_mpz (tmp, invsqrtA, F, pk, nullptr, invpk); /* tmp = invsqrtA^2 */
    if (verbose)
#pragma omp critical
      {
        fmt::print (stderr, "Alg({}):    mpz_poly_sqr_mod_f_mod_mpz took {:.2f}s (wct {:.2f}s, peak {}M)\n",
		 numdep, seconds () - st, wct_seconds () - wct,
		 PeakMemusage () >> 10U);
        fflush (stderr);
      }

    /* Faster version which computes x <- x + x/2*(1-a*x^2).
       However I don't see how to use the fact that the coefficients
       if 1-a*x^2 are divisible by p^(k/2). */
    st = seconds ();
    wct = wct_seconds ();
    pinf->mpz_poly_mul_mod_f_mod_mpz (tmp, tmp, a, F, pk, nullptr, invpk); /* tmp=a*invsqrtA^2 */
    if (verbose)
#pragma omp critical
      {
        fmt::print (stderr, "Alg({}):    mpz_poly_mul_mod_f_mod_mpz took {:.2f}s (wct {:.2f}s, peak {}M)\n",
		 numdep, seconds () - st, wct_seconds () - wct,
		 PeakMemusage () >> 10U);
        fflush (stderr);
      }
    mpz_poly_sub_ui (tmp, tmp, 1); /* a*invsqrtA^2-1 */
    pinf->mpz_poly_div_2_mod_mpz (tmp, tmp, pk); /* (a*invsqrtA^2-1)/2 */
    st = seconds ();
    wct = wct_seconds ();
    pinf->mpz_poly_mul_mod_f_mod_mpz (tmp, tmp, invsqrtA, F, pk, nullptr, invpk);
    if (verbose)
#pragma omp critical
      {
        fmt::print (stderr, "Alg({}):    mpz_poly_mul_mod_f_mod_mpz took {:.2f}s (wct {:.2f}s, peak {}M)\n",
		 numdep, seconds () - st, wct_seconds () - wct,
		 PeakMemusage () >> 10U);
        fflush (stderr);
      }
    /* tmp = invsqrtA/2 * (a*invsqrtA^2-1) */
    pinf->mpz_poly_sub_mod_mpz (invsqrtA, invsqrtA, tmp, pk);
  } while (k < target_k);

  /* multiply by a to get an approximation of the square root */
  st = seconds ();
  wct = wct_seconds ();
  pinf->mpz_poly_mul_mod_f_mod_mpz (tmp, invsqrtA, a, F, pk, nullptr, invpk);
  mpz_clear (invpk);
  if (verbose)
#pragma omp critical
    {
      fmt::print (stderr, "Alg({}):    final mpz_poly_mul_mod_f_mod_mpz took {:.2f}s (wct {:.2f}s, peak {}M)\n",
	       numdep, seconds () - st, wct_seconds () - wct,
	       PeakMemusage () >> 10U);
      fflush (stderr);
    }
  mpz_poly_mod_center (tmp, pk);

  mpz_poly_base_modp_clear (P, logk0);

  mpz_poly_set(res->p, tmp);
  res->v = v;

  mpz_clear (pk);
  mpz_poly_clear(tmp);
  mpz_poly_clear (invsqrtA);
  mpz_poly_clear (a);

  const size_t sqrt_size = mpz_poly_sizeinbase (res->p, 2);
#pragma omp critical
  {
    fmt::print (stderr, "Alg({}): maximal sqrt bit-size = {} ({:.0f}%% of target size)\n",
	     numdep, sqrt_size, 100.0 * (double) sqrt_size / (double) target_size);
    fflush (stderr);
  }

  return target_k;
}

static unsigned long
FindSuitableModP (cxx_mpz_poly const & F, cxx_mpz const & N)
{
  unsigned long p = 2;
  const int dF = F->deg;

  modul_poly_t fp;

  modul_poly_init (fp, dF);
  prime_info pi;
  prime_info_init (pi);
  while (true)
    {
    int d;

    p = getprime_mt (pi);

    /* check p does not divide N */
    if (mpz_gcd_ui (nullptr, N, p) != 1)
      continue;

    /* check the leading coefficient of F does not vanish mod p */
    d = modul_poly_set_mod (fp, F, &p);
    if (d != dF)
      continue;

    /* check that F is irreducible mod p */
    if (modul_poly_is_irreducible (fp, &p))
      break;

#define MAXP 1000000
    if (p > MAXP)
      {
	fmt::print (stderr, "Error, found no suitable prime up to {}\n", MAXP);
	fmt::print (stderr, "See paragraph \"Factoring with SNFS\" in README\n");
	exit(EXIT_FAILURE);
      }
    }
  modul_poly_clear (fp);
  prime_info_clear (pi);

  return p;
}

/*
   Process dependencies numdep to numdep + nthreads - 1.
*/
static int
calculateSqrtAlg (std::string const & prefix, unsigned int numdep,
                  cxx_cado_poly const & cpoly, int side, cxx_mpz const & Np,
                  const mpz_poly_parallel_info * pinf)
{
  FILE *resfile;
  unsigned long p;
  unsigned long nab = 0, nfree = 0;

  ASSERT_ALWAYS(side == 0 || side == 1);

  const std::string sidename = get_depsidename (prefix, numdep, side);

  // Init F to be the corresponding polynomial
  cxx_mpz_poly F(cpoly->pols[side]);
  cxx_mpz_poly F_hat;
  cxx_mpz_polymodF prod;

  /* create F_hat, the minimal polynomial of alpha_hat = lc(F) * alpha */
  mpz_poly_to_monic(F_hat, F);

  // Accumulate product with a subproduct tree
  {
      const cxx_mpz_polymodF_functions M(F, pinf);

      const std::string message = fmt::format("Alg({})", numdep);
      std::vector<cxx_mpz_polymodF> prd = read_ab_pairs_from_depfile(
              get_depname (prefix, "", numdep),
              M,
              message,
              nab,
              nfree);

      {
          cxx_mpz_poly scale, alpha_hat;
          mpz_poly_set_xi(alpha_hat, 1);
          mpz_poly_mul_mpz(alpha_hat, alpha_hat, mpz_poly_lc(F));
          mpz_poly_eval_diff_poly(scale, F_hat, alpha_hat);
          /* Add these two to the subproduct tree. ok, they're slightly
           * larger than the other elements, but it won't make a big
           * difference in the long run, and is certainly not worse than
           * doing the full unbalanced mul after the accumulation */
          prd.emplace_back(scale,0);
          prd.emplace_back(scale,0);
      }

      ASSERT_ALWAYS ((nab & 1U) == 0);
      ASSERT_ALWAYS ((nfree & 1U) == 0);
      /* nfree being even is forced by a specific character column added
       * by character.c. The correspond assert should not fail.
       *
       * nab being even is kind of a mystery. There is a character column
       * that gives the sign of the rational side. It could well be that
       * with the parameters we usually use, it is negative for all
       * relations; that would force the overall number of relations to
       * be even. Another possibility is when f_d contains a big prime
       * number that does not occur anywhere else, so that the power of
       * this prime is controlled by the usual characters, and since f_d
       * is always present...
       *
       * But: I wouldn't be surprised if the assert(even(nab)) fails.
       * Then, a patch would be:
       *    - remove the assert (!)
       *    - in the numerator of the big product, eliminate powers of
       *       f_d that divides all coefficients.
       *    - this should finally give an even power of f_d in the
       *      denominator, and the algorithm can continue.
       */

      prod = accumulate(prd, M, message);
  }

    p = FindSuitableModP(F, Np);
#pragma omp critical
    {
      fmt::print (stderr, "Alg({}): finished accumulating product at {:.2f}s (wct {:.2f}s)\n",
	       numdep, seconds(), wct_seconds () - wct0);
      fmt::print (stderr, "Alg({}): nab = {}, nfree = {}, v = {}\n", numdep,
	       nab, nfree, prod->v);
      fmt::print (stderr, "Alg({}): maximal polynomial bit-size = {}\n", numdep,
	       (unsigned long) mpz_poly_sizeinbase (prod->p, 2));
      fmt::print (stderr, "Alg({}): using p={} for lifting\n", numdep, p);
      fflush (stderr);
    }


#ifdef DEBUG
    cxx_mpz_poly prod0;
    unsigned long v0 = prod.v;
    ASSERT_ALWAYS(prod.p->deg == F->deg - 1);
    mpz_poly_set (prod0, prod.p);
#endif
    unsigned long target_k MAYBE_UNUSED;
    target_k = cxx_mpz_polymodF_sqrt (prod, prod, F, p, numdep, pinf);
#ifdef DEBUG
    unsigned long v = prod.v;
    /* we should have prod.p/fd^v = sqrt(prod0/fd^v0) mod (F,p^k)
       thus prod.p^2/fd^(2v) = prod0/fd^v0 mod (F,p^k)
       thus fd^v0*prod.p^2 = fd^(2v)*prod0 mod (F,p^k) */
    cxx_mpz_poly q;
    cxx_mpz pk, fdv0, fd2v;
    mpz_ui_pow_ui (pk, p, target_k);
    mpz_poly_sqr_mod_f_mod_mpz (q, prod.p, F, pk, nullptr, nullptr);
    /* we should have fd^v0*q = fd^(2v)*prod0 mod p^k */
    mpz_powm_ui (fdv0, F->coeff[F->deg], v0, pk);
    mpz_poly_mul_mpz (q, q, fdv0);
    mpz_poly_mod_mpz (q, q, pk, nullptr);
    /* now we should have q = fd^(2v)*prod0 mod p^k */
    mpz_powm_ui (fd2v, F->coeff[F->deg], 2 * v, pk);
    mpz_poly_mul_mpz (prod0, prod0, fd2v);
    mpz_poly_mod_mpz (prod0, prod0, pk, nullptr);
    /* now we should have q = prod0 */
    ASSERT_ALWAYS(mpz_poly_cmp (q, prod0) == 0);
#endif

#pragma omp critical
    {
        fmt::print (stderr,
                "Alg({}): square root lifted at {:.2f}s (wct {:.2f}s)\n",
                numdep, seconds(), wct_seconds () - wct0);
        fflush (stderr);
    }

    cxx_mpz algsqrt, aux;
    cxx_mpz m;
    int ret;

    cxx_mpz Np1 = Np;

    do {
        ret = cado_poly_getm(m, cpoly, Np1);
        if (!ret) {
            fmt::print(stderr, "When trying to compute m, got the factor {}\n", m);
            mpz_divexact(Np1, Np1, m);
        }
    } while (!ret);
    mpz_poly_eval_mod_mpz(algsqrt, prod->p, m, Np1);

    /* now divide by f'(alpha_hat), but mod N. alpha_hat mod N is lc*m
     */
    {
        cxx_mpz scale_modN, alpha_hat_modN;
        mpz_mul(alpha_hat_modN, mpz_poly_lc(F), m);
        mpz_poly_eval_diff(scale_modN, F_hat, alpha_hat_modN);
        mpz_invert(scale_modN, scale_modN, Np1);
        mpz_mul(algsqrt, algsqrt, scale_modN);
        mpz_mod(algsqrt, algsqrt, Np1);
    }

    mpz_invert(aux, mpz_poly_lc(F), Np1);  // 1/fd mod n
    mpz_powm_ui(aux, aux, prod->v, Np1);      // 1/fd^v mod n
    mpz_mul(algsqrt, algsqrt, aux);
    mpz_mod(algsqrt, algsqrt, Np1);

    resfile = fopen_maybe_compressed_lock (sidename.c_str(), "wb");
    fmt::print (resfile, "{}\n", algsqrt);
    fclose_maybe_compressed_lock (resfile, sidename.c_str());

#pragma omp critical
    {
        fmt::print (stderr, "Alg({}): square root is: {}\n",
		   numdep, algsqrt);
        fmt::print (stderr, "Alg({}): square root completed at {:.2f}s (wct {:.2f}s)\n",
	       numdep, seconds(), wct_seconds () - wct0);
        fflush (stderr);
    }
    return 0;
}

/****** Sqrt for quadratic sieve ******/

/* This is where we store side-relative functions for the square root for the
 * quadratic sieve.
 */
struct cxx_mpz_qs_functions {
    cxx_mpz_poly const & P;
    cxx_mpz const & Np;
    using T = std::pair<cxx_mpz, cxx_mpz>;
    static void set1(T & x) {
        mpz_set_ui(x.first, 1u);
        mpz_set_ui(x.second, 1u);
    }
    static void set(T & y, T const & x) {
        mpz_set(y.first, x.first);
        mpz_set(y.second, x.second);
    }
    static bool is1(T & x) {
        return mpz_cmp_ui(x.first, 1u) == 0 && mpz_cmp_ui(x.second, 1u) == 0;
    }
    void operator()(T & res, T const & u, T const & v) const {
        mpz_mul(res.first, u.first, v.first);
        mpz_mul(res.second, u.second, v.second);
        mpz_mod(res.second, res.second, Np);
    }
    ATTRIBUTE_NODISCARD
    T from_ab(cxx_mpz const& a, cxx_mpz const& b) const
    {
        /* return P(a) and a */
        ASSERT_ALWAYS(mpz_cmp_ui(b, 1u) == 0);
        cxx_mpz v;
        mpz_poly_eval(v, P, a);
        return { v, a };
    }
};
static int
calculateSqrtQS(
        std::string const & prefix,
        unsigned int numdep,
        cxx_cado_poly const & cpoly,
        cxx_mpz const & Np)
{
    ASSERT_ALWAYS(cpoly->nb_polys == 1 && mpz_poly_is_monic(cpoly->pols[0]));
    unsigned long nab = 0, nfree = 0;

#pragma omp critical
    {
#ifdef __MPIR_VERSION
        fmt::print (stderr, "Using MPIR {}\n", mpir_version);
#else
        fmt::print (stderr, "Using GMP {}\n", gmp_version);
#endif
        fflush (stderr);
    }

    cxx_mpz_poly F(cpoly->pols[0]);
    std::pair<cxx_mpz, cxx_mpz> prod;

    {
        const cxx_mpz_qs_functions M(F, Np);
        const std::string message = fmt::format("QS({})", numdep);

        std::vector<decltype(M)::T> prd = read_ab_pairs_from_depfile(
            get_depname (prefix, "", numdep),
            M,
            message,
            nab,
            nfree);

        prod = accumulate(prd, M, message);
    }

#pragma omp critical
    {
        fmt::print (stderr, "QS({}): size of product = {} bits (peak {}M)\n",
                            numdep, mpz_sizeinbase (prod.first, 2),
                            PeakMemusage () >> 10U);
        fflush (stderr);
    }

    if (mpz_sgn(prod.first) < 0)
    {
        fmt::print(stderr, "Error, product is negative: try another "
                            "dependency\n");
        exit(EXIT_FAILURE);
    }

#pragma omp critical
    {
        fmt::print(stderr, "QS({}): starting rational square root at "
                           "{:.2f}s (wct {:.2f}s)\n",
                           numdep, seconds (), wct_seconds () - wct0);
        fflush(stderr);
    }

    cxx_mpz v;
    /* since we know we have a square, take the square root */
    mpz_sqrtrem (prod.first, v, prod.first);

#pragma omp critical
    {
        fmt::print(stderr, "QS({}): computed square root at {:.2f}s "
                           "(wct {:.2f}s)\n",
                           numdep, seconds (), wct_seconds () - wct0);
        fflush(stderr);
    }

    if(mpz_cmp_ui (v, 0) != 0)
    {
        unsigned long p = 2, e, errors = 0;
        cxx_mpz pp;

        fmt::print(stderr, "Error, rational square root remainder is not "
                           "zero\n");
        /* reconstruct the initial value of prod to debug */
        mpz_mul(prod.first, prod.first, prod.first);
        mpz_add(prod.first, prod.first, v);
        prime_info pi;
        prime_info_init(pi);
        while (mpz_cmp_ui(prod.first, 1) > 0)
        {
            e = 0;
            if (verbose)
                fmt::print("Removing p={}:", p);
            mpz_set_ui (pp, p);
            e = mpz_remove (prod.first, prod.first, pp);
            if (verbose)
                fmt::print(" exponent={}, remaining {} bits\n", e,
                      mpz_sizeinbase (prod.first, 2));
            if ((e % 2) != 0)
            {
                errors ++;
                fmt::print(stderr, "Prime {} appears to odd power {}\n", p, e);
                if (verbose || errors >= 10)
                  break;
            }
            p = getprime_mt(pi);
        }
        prime_info_clear(pi);
        exit(EXIT_FAILURE);
    }

    mpz_mod(prod.first, prod.first, Np);

#pragma omp critical
    {
        fmt::print(stderr, "QS({}): reduced mod n at {:.2f}s (wct {:.2f}s)\n",
                           numdep, seconds (), wct_seconds () - wct0);
        fflush(stderr);
    }

    for (size_t i = 0; i < 2; ++i) {
        const std::string sidename = get_depsidename(prefix, numdep, i);
        FILE *resfile = fopen_maybe_compressed_lock(sidename.c_str(), "wb");
        fmt::print(resfile, "{}\n", i == 0 ? prod.first : prod.second);
        fclose_maybe_compressed_lock (resfile, sidename.c_str());
    }

#pragma omp critical
    {
        fmt::print(stderr, "QS({}): square root is {}\n", numdep, prod.first);
        fmt::print(stderr, "QS({}): square root time: {:.2f}s (wct {:.2f}s)\n",
                           numdep, seconds (), wct_seconds () - wct0);
        fmt::print(stderr, "QS({}): product of a is {}\n", numdep, prod.second);
        fflush(stderr);
    }

  return 0;
}

/*
 * Try to factor input using trial division up to bound B.
 * Found factors are printed (one per line).
 * Returns 1 if input is completely factored, otherwise, returns
 * remaining factor.
 */
static unsigned long
trialdivide_print(unsigned long N, unsigned long B)
{
    ASSERT(N != 0);
    if (N == 1) return 1;
    unsigned long p;
    prime_info pi;
    prime_info_init (pi);
    for (p = 2; p <= B; p = getprime_mt (pi)) {
        while ((N%p) == 0) {
            N /= p;
            fmt::print("{}\n", p);
            if (N == 1) {
                prime_info_clear (pi);
                return N;
            }
        }
    }
    prime_info_clear (pi);
    return N;
}

static void print_nonsmall(cxx_mpz const & zx)
{
    if (mpz_probab_prime_p(zx, 10))
        fmt::print("{}\n", zx);
    else {
        int pp = mpz_perfect_power_p(zx);
        if (pp) {
            pp = mpz_sizeinbase(zx, 2);
            cxx_mpz roo;
            while (!mpz_root(roo, zx, pp))
                pp--;
            int i;
            for (i = 0; i < pp; ++i)
                fmt::print("{}\n", roo);
        } else
            fmt::print("{}\n", zx);
    }
    fflush (stdout);
}

static void print_factor(cxx_mpz const & N)
{
    unsigned long xx = mpz_get_ui(N);
    if (mpz_cmp_ui(N, xx) == 0) {
        xx = trialdivide_print(xx, 1000000);
        if (xx != 1) {
            mpz_t zx;
            mpz_init(zx);
            mpz_set_ui(zx, xx);
            print_nonsmall(zx);
            mpz_clear(zx);
        }
    } else
        print_nonsmall(N);
}


/********** GCD **********/
static int
calculateGcd (std::string const & prefix, unsigned int numdep, cxx_mpz const & Np)
{
    std::string sidename[2];
    FILE *sidefile[2] = { nullptr, nullptr };
    cxx_mpz sidesqrt[2];
    int found = 0;

    for (int side = 0; side < 2; ++side) {
        sidename[side] = get_depsidename (prefix, numdep, side);
        sidefile[side] = fopen_maybe_compressed_lock (sidename[side].c_str(), "rb");
        if (!sidefile[side]) {
            const std::lock_guard<std::mutex> dummy(stdio_guard);
            fmt::print(stderr, "Error, cannot open file {} for reading\n",
                    sidename[side]);
            exit(EXIT_FAILURE);
        }
        gmp_fscanf (sidefile[side], "%Zd", (mpz_ptr) sidesqrt[side]);
        fclose_maybe_compressed_lock (sidefile[side], sidename[side].c_str());
    }

    cxx_mpz g1, g2;

    // reduce mod Np
    mpz_mod(sidesqrt[0], sidesqrt[0], Np);
    mpz_mod(sidesqrt[1], sidesqrt[1], Np);

    // First check that the squares agree
    mpz_mul(g1, sidesqrt[0], sidesqrt[0]);
    mpz_mod(g1, g1, Np);

    mpz_mul(g2, sidesqrt[1], sidesqrt[1]);
    mpz_mod(g2, g2, Np);

    if (mpz_cmp(g1, g2)!=0) {
      const std::lock_guard<std::mutex> dummy(stdio_guard);
      fmt::print(stderr, "Bug: the squares do not agree modulo n!\n");
      ASSERT_ALWAYS(0);
      //      fmt::print("g1:={};\ng2:={};\n", g1, g2);
    }

    mpz_sub(g1, sidesqrt[0], sidesqrt[1]);
    mpz_gcd(g1, g1, Np);
    if (mpz_cmp(g1,Np)) {
      if (mpz_cmp_ui(g1,1)) {
        found = 1;
        const std::lock_guard<std::mutex> dummy(stdio_guard);
        print_factor(g1);
      }
    }

    mpz_add(g2, sidesqrt[0], sidesqrt[1]);
    mpz_gcd(g2, g2, Np);
    if (mpz_cmp(g2,Np)) {
      if (mpz_cmp_ui(g2,1)) {
        found = 1;
        const std::lock_guard<std::mutex> dummy(stdio_guard);
        print_factor(g2);
      }
    }

    if (!found) {
      const std::lock_guard<std::mutex> dummy(stdio_guard);
      fmt::print ("Failed\n");
    }

    return 0;
}

typedef struct
{
  uint64_t *abs;
  uint64_t *dep_masks;
  unsigned int *dep_counts;
  unsigned int nonzero_deps;
  FILE **dep_files;
} sqrt_data_t;

template<filter_io_config cfg>
static void *
thread_sqrt (void * context_data, typename cfg::rel_ptr rel)
{
  auto data = (sqrt_data_t *) context_data;
  for(unsigned int j = 0 ; j < data->nonzero_deps ; j++)
  {
    if (data->abs[rel->num] & data->dep_masks[j])
    {
      if constexpr (std::is_same_v<cfg, filter_io_large_ab_cfg>) {
        fmt::print(data->dep_files[j], "{} {}\n",
                   cxx_mpz(rel->a), cxx_mpz(rel->b));
      } else {
        fmt::print(data->dep_files[j], "{} {}\n", rel->a, rel->b);
      }
      data->dep_counts[j]++;
    }
  }
  return nullptr;
}

template<filter_io_config cfg>
static void create_dependencies(const char * prefix, const char * indexname, const char * purgedname, const char * kername)
{
    FILE * ix = fopen_maybe_compressed(indexname, "r");
    uint64_t small_nrows;
    int ret;

    ret = fscanf(ix, "%" SCNu64 "\n", &small_nrows);
    ASSERT(ret == 1);

    FILE * ker;
    size_t ker_stride;
    /* Check that kername has consistent size */
    {
        ker = fopen(kername, "rb");
        if (!ker) {
            perror(kername);
            exit(errno);
        }
        struct stat sbuf[1];
        ret = fstat(fileno(ker), sbuf);
        if (ret < 0) {
            perror(kername);
            exit(errno);
        }
        ASSERT_ALWAYS(sbuf->st_size % small_nrows == 0);
        const unsigned int ndepbytes = sbuf->st_size / small_nrows;
        fmt::print(stderr, "{} contains {} dependencies (including padding)\n",
                kername, 8 * ndepbytes);
        ker_stride = ndepbytes - sizeof(uint64_t);
        if (ker_stride)
            fmt::print(stderr, "Considering only the first 64 dependencies\n");
    }

    /* Read the number of (a,b) pairs */
    uint64_t nrows, ncols;
    purgedfile_read_firstline (purgedname, &nrows, &ncols);

    std::unique_ptr<uint64_t[]> const abs(new uint64_t[nrows]);
    std::fill(abs.get(), abs.get() + nrows, 0);

    for(uint64_t i = 0 ; i < small_nrows ; i++) {
        uint64_t v;
        ret = (int) fread(&v, sizeof(uint64_t), 1, ker);
        if (ker_stride) fseek(ker, ker_stride, SEEK_CUR);

        /* read the index row */
        int nc;
        ret = (int) fscanf(ix, "%d", &nc); ASSERT_ALWAYS(ret == 1);
        for(int k = 0 ; k < nc ; k++) {
            uint64_t col;
            ret = (int) fscanf(ix, "%" SCNx64 "", &col); ASSERT_ALWAYS(ret == 1);
            ASSERT_ALWAYS(col < nrows);
            abs[col] ^= v;
        }
    }
    fclose_maybe_compressed(ix, indexname);
    fclose(ker);

    unsigned int nonzero_deps = 0;
    uint64_t sanity = 0;
    for(uint64_t i = 0 ; i < nrows ; i++) {
        sanity |= abs[i];
    }
    uint64_t dep_masks[64]={0,};
    std::string dep_names[64];
    FILE * dep_files[64];
    unsigned int dep_counts[64]={0,};

    for(unsigned int i = 0 ; i < 64U ; i++) {
        uint64_t const m = UINT64_C(1) << i;
        if (sanity & m)
            dep_masks[nonzero_deps++] = m;
    }
    fmt::print(stderr, "Total: {} non-zero dependencies\n", nonzero_deps);
    for(unsigned int i = 0 ; i < nonzero_deps ; i++) {
        dep_names[i] = get_depname (prefix, "", i);
        dep_files[i] = fopen_maybe_compressed (dep_names[i].c_str(), "wb");
        ASSERT_ALWAYS(dep_files[i]);
    }

    sqrt_data_t data = {.abs = abs.get(), .dep_masks = dep_masks,
                        .dep_counts = dep_counts, .nonzero_deps = nonzero_deps,
                        .dep_files = dep_files};
    std::vector<std::string> fic {purgedname};
    filter_rels<cfg>(fic, thread_sqrt<cfg>, &data, EARLYPARSE_NEED_AB_HEXA,
                     nullptr, nullptr);

    fmt::print(stderr, "Written {} dependencies files\n", nonzero_deps);
    for(unsigned int i = 0 ; i < nonzero_deps ; i++) {
        fmt::print(stderr, "{} : {} (a,b) pairs\n", dep_names[i], dep_counts[i]);
        fclose_maybe_compressed (dep_files[i], dep_names[i].c_str());
    }
}


#define TASK_SQRT 0
#define TASK_GCD  2
#define TASK_QS 4
/* perform one task (rat or alg or gcd) on one dependency */
static void
one_thread (std::string const & prefix, int task, unsigned int numdep, cxx_cado_poly const & cpoly, int side, cxx_mpz const & Np, const mpz_poly_parallel_info * pinf)
{
  if (task == TASK_SQRT) {
      if (cpoly->pols[side]->deg == 1) {
          calculateSqrtRat (prefix, numdep, cpoly, side, Np);
      } else {
          calculateSqrtAlg (prefix, numdep, cpoly, side, Np, pinf);
      }
  } else if (task == TASK_QS) {
    ASSERT_ALWAYS(side == 0);
    calculateSqrtQS(prefix, numdep, cpoly, Np);
  } else /* gcd */
    calculateGcd (prefix, numdep, Np);
}

/* process task (0=sqrt, 2=gcd) in parallel for
   dependencies numdep to numdep + nthreads - 1 */
static void
calculateTaskN (int task, std::string const & prefix, unsigned int numdep, int nthreads,
                cxx_cado_poly const & cpoly, int side, cxx_mpz const & Np)
{
  /* This descriptor will hold info about "how" we parallelize the
   * mpz_poly operations. For the moment, there's not a lot in there.
   * Only the fact that the pointer that we pass is not NULL is
   * significant! But eventually, we may find it useful to add more
   * stuff.
   */
  const mpz_poly_parallel_info pinf;

  omp_set_num_threads(iceildiv(omp_get_max_threads(), nthreads));

  std::vector<std::thread> threads;
  threads.reserve(nthreads);

#if defined(__OpenBSD__) || defined(HAVE_MUSL)
  /* On openbsd, we have obscure failures that seem to be triggered
   * by multithreading. So let's play it simple.
   *
   * We seem to experience exactly a similar problem with musl libc (used
   * by alpine linux). Note that HAVE_MUSL is a custom flag that we define
   * ourselves.
   */
  for (int j = 0; j < nthreads; j++)
      one_thread(prefix, task, numdep + j, cpoly, side, Np, &pinf);
#else
  for (int j = 0; j < nthreads; j++)
      threads.emplace_back(one_thread, prefix, task, numdep + j, cpoly, side, Np, & pinf);

  for(auto & t : threads) t.join();
#endif
}

static void declare_usage(param_list pl)
{
    param_list_decl_usage(pl, "poly", "Polynomial file");
    param_list_decl_usage(pl, "purged", "Purged relations file, as produced by 'purge'");
    param_list_decl_usage(pl, "index", "Index file, as produced by 'merge'");
    param_list_decl_usage(pl, "ker", "Kernel file, as produced by 'characters'");
    param_list_decl_usage(pl, "prefix", "File name prefix used for output files");
    param_list_decl_usage(pl, "ab", "For each dependency, create file with the a,b-values of the relations used in that dependency");
    param_list_decl_usage(pl, "side0", "Compute square root for side 0 and store in file");
    param_list_decl_usage(pl, "side1", "Compute square root for side 1 and store in file");
    param_list_decl_usage(pl, "qs", "Compute square root for quadratic sieve and store in file");
    param_list_decl_usage(pl, "gcd", "Compute gcd of the two square roots. Requires square roots on both sides");
    param_list_decl_usage(pl, "dep", "The initial dependency for which to compute square roots");
    param_list_decl_usage(pl, "t",   "The number of dependencies to process (default 1)");
    param_list_decl_usage(pl, "v", "More verbose output");
    param_list_decl_usage(pl, "large-ab", "enable support for a and b larger "
                                          "than 64 bits");
    param_list_decl_usage(pl, "force-posix-threads", "force the use of posix threads, do not rely on platform memory semantics");
}

static void usage(param_list pl, const char * argv0, FILE *f)
{
    param_list_print_usage(pl, argv0, f);
    fmt::print(f, "Usage: {} [-ab || -side0 || -side1 || -gcd] -poly polyname -prefix prefix -dep numdep -t ndep", argv0);
    fmt::print(f, " -purged purgedname -index indexname -ker kername\n");
    fmt::print(f, "or {} (-side0 || -side1 || -gcd) -poly polyname -prefix prefix -dep numdep -t ndep\n\n", argv0);
    fmt::print(f, "(a,b) pairs of dependency relation 'numdep' will be r/w in file 'prefix.numdep',");
    fmt::print(f, " side0 sqrt in 'prefix.side0.numdep' ...\n");
    exit(EXIT_FAILURE);
}

// coverity[root_function]
int main(int argc, char const *argv[])
{
    unsigned int numdep = UINT_MAX;
    int nthreads = 1, ret MAYBE_UNUSED, i;
    int largeab = 0;

    char const * me = *argv;
    /* print the command line */
    fmt::print (stderr, "{}.r{}", argv[0], cado_revision_string);
    for (i = 1; i < argc; i++)
      fmt::print (stderr, " {}", argv[i]);
    fmt::print (stderr, "\n");

    cxx_param_list pl;
    declare_usage(pl);

    int opt_ab = 0;    /* create dependency files */
    int opt_side0 = 0; /* compute square root on side 0 */
    int opt_side1 = 0; /* compute square root on side 1 */
    int opt_gcd = 0;   /* compute gcd */
    int opt_qs = 0; /* compute root for quadratic sieve */
    param_list_configure_switch(pl, "ab", &opt_ab);
    param_list_configure_switch(pl, "side0", &opt_side0);
    param_list_configure_switch(pl, "side1", &opt_side1);
    param_list_configure_switch(pl, "qs", &opt_qs);
    param_list_configure_switch(pl, "gcd", &opt_gcd);
    param_list_configure_switch(pl, "-v", &verbose);
    param_list_configure_switch(pl, "force-posix-threads", &filter_rels_force_posix_threads);
    param_list_configure_switch(pl, "large-ab", &largeab);
    argc--,argv++;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) continue;
        if (strcmp(*argv, "--help") == 0) {
            usage(pl, me, stderr);
            return EXIT_SUCCESS;
        } else {
            fmt::print(stderr, "unexpected argument: {}\n", *argv);
            usage(pl, me, stderr);
            return EXIT_FAILURE;
        }
    }
    const char * tmp;
    if(!(tmp = param_list_lookup_string(pl, "poly"))) {
        fmt::print(stderr, "Parameter -poly is missing\n");
        usage(pl, me, stderr);
        return EXIT_FAILURE;
    }

    cxx_cado_poly cpoly;

    ret = cado_poly_read(cpoly, tmp);
    if (ret == 0) {
        fmt::print(stderr, "Could not read polynomial file\n");
        return EXIT_FAILURE;
    }

    if (cpoly->nb_polys < 1 || cpoly->nb_polys > 2) {
        fmt::print(stderr, "Error: number of polys should be 1 or 2, got {}\n",
                           cpoly->nb_polys);
        exit (EXIT_FAILURE);
    } else if (cpoly->nb_polys == 1 && opt_side1) {
        fmt::print(stderr, "Error: -side1 is not compatible with only one "
                           "side\n");
        exit (EXIT_FAILURE);
    }

    param_list_parse_uint (pl, "dep", &numdep);
    param_list_parse_int (pl, "t", &nthreads);
    const char * purgedname = param_list_lookup_string(pl, "purged");
    const char * indexname = param_list_lookup_string(pl, "index");
    const char * kername = param_list_lookup_string(pl, "ker");
    const char * prefix = param_list_lookup_string(pl, "prefix");
    if (!prefix) {
        fmt::print(stderr, "Parameter -prefix is missing\n");
        return EXIT_FAILURE;
    }
    if (param_list_warn_unused(pl))
        return EXIT_FAILURE;

    /* if no options then -ab -side0 -side1 -gcd */
    if (!(opt_ab || opt_side0 || opt_side1 || opt_qs || opt_gcd)) {
        opt_ab = opt_side0 = opt_gcd = 1;
        opt_side1 = (cpoly->nb_polys == 2);
    } else if (opt_qs) {
        if (cpoly->nb_polys != 1) {
            fmt::print(stderr, "Error: -qs is only valid for one-sided poly\n");
            exit (EXIT_FAILURE);
        } else if (opt_side0 || opt_side1) {
            fmt::print(stderr, "Error: -qs is not compatible with -side0 "
                               "and/or -side1\n");
            exit (EXIT_FAILURE);
        }
    }

    const double cpu0 = seconds ();
    wct0 = wct_seconds();

    /*
     * In the case where the number N to factor has a prime factor that
     * divides the leading coefficient of f or g, the reduction modulo N
     * will fail. Let's compute N', the factor of N that is coprime to
     * those leading coefficients.
     */
    cxx_mpz Np;
    {
        cxx_mpz gg;
        mpz_set(Np, cpoly->n);
        for (int side = 0; side < cpoly->nb_polys; ++side) {
            do {
                mpz_gcd(gg, Np, mpz_poly_lc(cpoly->pols[side]));
                if (mpz_cmp_ui(gg, 1) != 0) {
                    fmt::print(stderr, "Warning: found the following factor of N as a factor of g: {}\n", gg);
                    print_factor(gg);
                    mpz_divexact(Np, Np, gg);
                }
            } while (mpz_cmp_ui(gg, 1) != 0);
        }
        /* Trial divide Np, to avoid bug if a stupid input is given */
        {
            unsigned long p;
            prime_info pi;
            prime_info_init (pi);
            for (p = 2; p <= 1000000; p = getprime_mt (pi)) {
                while (mpz_tdiv_ui(Np, p) == 0) {
                    fmt::print("{}\n", p);
                    mpz_divexact_ui(Np, Np, p);
                }
            }
            prime_info_clear (pi);
        }
        if (mpz_cmp(cpoly->n, Np) != 0)
            fmt::print(stderr, "Now factoring N' = {}\n", Np);
        if (mpz_cmp_ui(Np, 1) == 0) {
            fmt::print(stderr, "Hey N' is 1! Stopping\n");
            return 0;
        }
        if (mpz_probab_prime_p(Np, 10) || mpz_perfect_power_p(Np)) {
            fmt::print(stderr, "Hey N' is (power of) prime! Stopping\n");
            print_factor(Np);
            return 0;
        }
    }

    if (opt_ab) {
        /* Computing (a,b) pairs is now done in batch for 64 dependencies
         * together -- should be enough for our purposes, even if we do
         * have more dependencies !
         */
        if (!indexname) {
            fmt::print(stderr, "Parameter -index is missing\n");
            return EXIT_FAILURE;
        }
        if (!purgedname) {
            fmt::print(stderr, "Parameter -purged is missing\n");
            return EXIT_FAILURE;
        }
        if (!kername) {
            fmt::print(stderr, "Parameter -ker is missing\n");
            return EXIT_FAILURE;
        }
        if (!largeab) {
            create_dependencies<filter_io_default_cfg>(prefix, indexname,
                                                       purgedname, kername);
        } else {
            create_dependencies<filter_io_large_ab_cfg>(prefix, indexname,
                                                        purgedname, kername);
        }
    }

#ifdef __OpenBSD__
    if (nthreads > 1) {
        fmt::print(stderr, "Warning: reducing number of threads to 1 for openbsd ; unexplained failure https://ci.inria.fr/cado/job/compile-openbsd-59-amd64-random-integer/2775/console\n");
        /* We'll still process everything we've been told to. But in a
         * single-threaded fashion */
    }
#endif

    if (opt_side0 || opt_side1 || opt_qs || opt_gcd)
      {
        int i;

        for (i = 0; i < nthreads; i++)
          if (check_dep (prefix, numdep + i) == 0)
            {
              fmt::print (stderr, "Warning: dependency {} does not exist, reducing the number of threads to {}\n",
                       numdep + i, i);
              nthreads = i;
              break;
            }
      }

    if (nthreads == 0)
      {
        fmt::print (stderr, "Error, no more dependency\n");
        return 1;
      }

    if (opt_side0) {
        ASSERT_ALWAYS(numdep != UINT_MAX);
        calculateTaskN(TASK_SQRT, prefix, numdep, nthreads, cpoly, 0, Np);
    }

    if (opt_side1) {
        ASSERT_ALWAYS(numdep != UINT_MAX);
        calculateTaskN(TASK_SQRT, prefix, numdep, nthreads, cpoly, 1, Np);
    }

    if (opt_qs) {
        ASSERT_ALWAYS(numdep != UINT_MAX);
        calculateTaskN(TASK_QS, prefix, numdep, nthreads, cpoly, 0, Np);
    }

    if (opt_gcd) {
        ASSERT_ALWAYS(numdep != UINT_MAX);
        calculateTaskN(TASK_GCD, prefix, numdep, nthreads, cpoly, 0, Np);
    }

    print_timing_and_memory (stderr, cpu0, wct0);
    return 0;
}
