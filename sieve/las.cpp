#define __STDCPP_MATH_SPEC_FUNCS__ 201003L
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1       /* for expint() */
/* the macro above is for #include <cmath> -- however it must happen
 * first, because it may well be that one of the intermediary headers
 * pull stuff that is dependent on this flag.
 */
#include "cado.h"
#include <stdint.h>     /* AIX wants it first (it's a bug) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h> /* for PRIx64 macro and strtoumax */
#include <cmath>   // for ceiling, floor in cfrac
#include <ctype.h>
#include <float.h>
#include <fcntl.h>  /* for _O_BINARY */
#include <stdarg.h> /* Required so that GMP defines gmp_vfprintf() */
#include <algorithm>
#include <vector>
#include <sstream>  /* for c++ string handling */
#include <iterator> /* ostream_iterator */
#include <thread>   /* we use std::thread for las sub-jobs */
#include <iomanip>  /* std::setprecision, std::fixed */
#include "threadpool.hpp"
#include "fb.hpp"
#include "portability.h"
#include "utils.h"           /* lots of stuff */
#include "relation.hpp"
#include "ecm/facul.hpp"
#include "bucket.hpp"
#include "trialdiv.hpp"
#include "las-config.h"
#include "las-info.hpp"
#include "las-coordinates.hpp"
#include "las-debug.hpp"
#include "las-duplicate.hpp"
#include "las-report-stats.hpp"
#include "las-norms.hpp"
#include "las-unsieve.hpp"
#include "las-arith.hpp"
#include "las-plattice.hpp"
#include "las-qlattice.hpp"
#include "las-smallsieve.hpp"
#include "las-descent-trees.hpp"
#include "las-cofactor.hpp"
#include "las-fill-in-buckets.hpp"
#include "las-threads.hpp"
#include "las-todo-entry.hpp"
#include "las-choose-sieve-area.hpp"
#include "las-auxiliary-data.hpp"
#include "las-process-bucket-region.hpp"
#include "las-galois.hpp"
#include "las-parallel.hpp"

#include "memusage.h"
#include "tdict.hpp"
#ifdef  DLP_DESCENT
#include "las-dlog-base.hpp"
#endif

#ifdef HAVE_SSE41
/* #define SSE_SURVIVOR_SEARCH 1 */
#include <smmintrin.h>
#endif

size_t base_memory = 0;
int recursive_descent = 0;
int prepend_relation_time = 0;
int exit_after_rel_found = 0;
int allow_largesq = 0;
int sync_at_special_q = 0;
int trialdiv_first_side = 0;

double general_grace_time_ratio = DESCENT_DEFAULT_GRACE_TIME_RATIO;

double tt_qstart;

/*****************************/

/*}}}*/


/*{{{ stuff related to las output: -out, -stats-stderr, and so on. */
struct las_augmented_output_channel {
    int verbose = 1;
    /* outputname is owned by pl, and output by the libc. We want to
     * control when these get released, with the release() function */
    FILE *output = NULL;
    const char * outputname = NULL; /* keep track of whether it's gzipped or not */
    void set(cxx_param_list & pl);
    void release();
    void fflush();
    static void declare_usage(cxx_param_list & pl);
    static void configure_aliases(cxx_param_list &) { }
    static void configure_switches(cxx_param_list & pl);
};

las_augmented_output_channel las_output;

static void las_verbose_enter(cxx_param_list & pl, FILE * output, int verbose)
{
    verbose_interpret_parameters(pl);
    verbose_output_init(NR_CHANNELS);
    verbose_output_add(0, output, verbose + 1);
    verbose_output_add(1, stderr, 1);
    /* Channel 2 is for statistics. We always print them to las' normal output */
    verbose_output_add(2, output, 1);
    if (param_list_parse_switch(pl, "-stats-stderr")) {
        /* If we should also print stats to stderr, add stderr to channel 2 */
        verbose_output_add(2, stderr, 1);
    }
#ifdef TRACE_K
    const char *trace_file_name = param_list_lookup_string(pl, "traceout");
    FILE *trace_file = stderr;
    if (trace_file_name != NULL) {
        trace_file = fopen(trace_file_name, "w");
        DIE_ERRNO_DIAG(trace_file == NULL, "fopen", trace_file_name);
    }
    verbose_output_add(TRACE_CHANNEL, trace_file, 1);
#endif
}

static void las_verbose_leave()
{
    verbose_output_clear();
}

void las_augmented_output_channel::fflush()
{
    verbose_output_start_batch();
    ::fflush(output);
    verbose_output_end_batch();
}

void las_augmented_output_channel::set(cxx_param_list & pl)
{
    ASSERT_ALWAYS(output == NULL);
    output = stdout;
    outputname = param_list_lookup_string(pl, "out");
    if (outputname) {
	if (!(output = fopen_maybe_compressed(outputname, "w"))) {
	    fprintf(stderr, "Could not open %s for writing\n", outputname);
	    exit(EXIT_FAILURE);
	}
    }
    verbose = param_list_parse_switch(pl, "-v");
    setvbuf(output, NULL, _IOLBF, 0);      /* mingw has no setlinebuf */
    las_verbose_enter(pl, output, verbose);

    param_list_print_command_line(output, pl);
    las_display_config_flags();
}

/* This cannot be a dtor because las_output is a global, and it holds
 * references to pl which is local -- Thus we want to control the exact
 * time where fclose is called.
 */
void las_augmented_output_channel::release()
{
    if (outputname)
        fclose_maybe_compressed(output, outputname);
    las_verbose_leave();
    outputname = NULL;
}


void las_augmented_output_channel::configure_switches(cxx_param_list & pl)
{
    param_list_configure_switch(pl, "-stats-stderr", NULL);
    param_list_configure_switch(pl, "-v", NULL);
}

void las_augmented_output_channel::declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "v",    "verbose mode, also prints sieve-area checksums");
    param_list_decl_usage(pl, "out",  "filename where relations are written, instead of stdout");
#ifdef TRACE_K
    param_list_decl_usage(pl, "traceout", "Output file for trace output, default: stderr");
#endif
    param_list_decl_usage(pl, "stats-stderr", "print stats to stderr in addition to stdout/out file");
}
/*}}}*/



/* {{{ apply_buckets */
template <typename HINT>
#ifndef TRACE_K
/* backtrace display can't work for static symbols (see backtrace_symbols) */
NOPROFILE_STATIC
#endif
void
apply_one_update (unsigned char * const S,
        bucket_update_t<1, HINT> const & u,
        const unsigned char logp, where_am_I & w)
{
  WHERE_AM_I_UPDATE(w, h, u.hint_for_where_am_i());
  WHERE_AM_I_UPDATE(w, x, u.x);
  sieve_increase(S + u.x, logp, w);
}

template <typename HINT>
#ifndef TRACE_K
/* backtrace display can't work for static symbols (see backtrace_symbols) */
NOPROFILE_STATIC
#endif
void
apply_one_bucket (unsigned char *S,
        bucket_array_t<1, HINT> const &BA, const int i,
        fb_factorbase::slicing::part const & fbp, where_am_I & w)
{
  WHERE_AM_I_UPDATE(w, p, 0);

  /* The slice indices are only fake indices in this case. So we must not
   * instantiate this code ! */
  static_assert(!std::is_same<HINT, logphint_t>::value, "This code must not be instantiated with hint type logphint_t");
  static_assert(!std::is_same<HINT, longhint_t>::value, "This code must not be instantiated with hint type longhint_t");

  for (slice_index_t i_slice = 0; i_slice < BA.get_nr_slices(); i_slice++) {
      auto sl = BA.slice_range(i, i_slice);
      auto it = sl.begin();
      auto it_end = sl.end();

    const slice_index_t slice_index = BA.get_slice_index(i_slice);
    const unsigned char logp = fbp[slice_index].get_logp();

    /* TODO: the code below is quite possibly correct, except perhaps for the
     * treatment of where_am_I stuff. I get inconsistent
     * reports, esp when I vary the number of threads.
     */
#if 1
    const bucket_update_t<1, HINT> *next_align;
    if (sizeof(bucket_update_t<1, HINT>) == 4) {
      next_align = (bucket_update_t<1, HINT> *) (((size_t) it + 0x3F) & ~((size_t) 0x3F));
      if (UNLIKELY(next_align > it_end)) next_align = it_end;
    } else {
      next_align = it_end;
    }

    while (it != next_align)
      apply_one_update<HINT> (S, *it++, logp, w);

    while (it + 16 <= it_end) {
      uint64_t x0, x1, x2, x3, x4, x5, x6, x7;
      uint16_t x;
#ifdef HAVE_SSE2
#if defined(__ICC) || defined(__INTEL_COMPILER)
      /* https://gcc.gnu.org/bugzilla/show_bug.cgi?id=45414 */
      _mm_prefetch(((const char *) it)+256, _MM_HINT_NTA);
#else
      _mm_prefetch(((unsigned char *) it)+256, _MM_HINT_NTA);
#endif
#endif
      x0 = ((uint64_t *) it)[0];
      x1 = ((uint64_t *) it)[1];
      x2 = ((uint64_t *) it)[2];
      x3 = ((uint64_t *) it)[3];
      x4 = ((uint64_t *) it)[4];
      x5 = ((uint64_t *) it)[5];
      x6 = ((uint64_t *) it)[6];
      x7 = ((uint64_t *) it)[7];
      it += 16;
      __asm__ __volatile__ (""); /* To be sure all x? are read together in one operation */
#ifdef CADO_LITTLE_ENDIAN
#define INSERT_2_VALUES(X) do {						\
	(X) >>= 16; x = (uint16_t) (X);					\
	WHERE_AM_I_UPDATE(w, x, x); sieve_increase(S + x, logp, w);	\
	(X) >>= 32;							\
	WHERE_AM_I_UPDATE(w, x, X); sieve_increase(S + (X), logp, w);	\
      } while (0);
#else
#define INSERT_2_VALUES(X) do {						\
	x = (uint16_t) (X);						\
	WHERE_AM_I_UPDATE(w, x, x); sieve_increase(S + x, logp, w);	\
	(X) >>= 32; x = (uint16_t) (X);					\
	WHERE_AM_I_UPDATE(w, x, x); sieve_increase(S + x, logp, w);	\
      } while (0);

#endif
      INSERT_2_VALUES(x0);
      INSERT_2_VALUES(x1);
      INSERT_2_VALUES(x2);
      INSERT_2_VALUES(x3);
      INSERT_2_VALUES(x4);
      INSERT_2_VALUES(x5);
      INSERT_2_VALUES(x6);
      INSERT_2_VALUES(x7);
    }
#endif
    while (it != it_end)
      apply_one_update<HINT> (S, *it++, logp, w);
  }
}

// Create the four instances, longhint_t and logphint_t are specialized.
template 
void apply_one_bucket<shorthint_t> (unsigned char *S,
        const bucket_array_t<1, shorthint_t> &BA, const int i,
        fb_factorbase::slicing::part const & fbp, where_am_I & w);

template 
void apply_one_bucket<emptyhint_t> (unsigned char *S,
        const bucket_array_t<1, emptyhint_t> &BA, const int i,
        fb_factorbase::slicing::part const & fbp, where_am_I & w);

template <>
void apply_one_bucket<longhint_t> (unsigned char *S,
        const bucket_array_t<1, longhint_t> &BA, const int i,
        fb_factorbase::slicing::part const & fbp, where_am_I & w)
{
    WHERE_AM_I_UPDATE(w, p, 0);
    // There is only one fb_slice. Slice indices are embedded in the hints.
    ASSERT(BA.get_nr_slices() == 1);
    ASSERT(BA.get_slice_index(0) == std::numeric_limits<slice_index_t>::max());
    for(auto const & it : BA.slice_range(i, 0)) {
        slice_index_t index = it.index;
        const unsigned char logp = fbp[index].get_logp();
        apply_one_update<longhint_t> (S, it, logp, w);
    }
}
template <>
void apply_one_bucket<logphint_t> (unsigned char *S,
        const bucket_array_t<1, logphint_t> &BA, const int i,
        fb_factorbase::slicing::part const &, where_am_I & w)
{
    WHERE_AM_I_UPDATE(w, p, 0);
    // There is only one fb_slice. logp's are embedded in the hints.
    ASSERT(BA.get_nr_slices() == 1);
    ASSERT(BA.get_slice_index(0) == std::numeric_limits<slice_index_t>::max());
    for(auto const & it : BA.slice_range(i, 0)) {
        apply_one_update<logphint_t> (S, it, it.logp, w);
    }
}
/* }}} */

/* {{{ Trial division */
typedef std::vector<uint64_t> factor_list_t;

static void 
factor_list_add(factor_list_t & fl, const uint64_t p)
{
    fl.push_back(p);
}

bool
factor_list_contains(factor_list_t const & fl, const uint64_t p)
{
    return std::find(fl.begin(), fl.end(), p) != fl.end();
}

// print a comma-separated list of factors.
// returns the number of factor printed (in particular, a comma is needed
// after this output only if the return value is non zero)
int factor_list_fprint(FILE *f, factor_list_t const & fl) {
    for(size_t i = 0; i < fl.size(); ++i) {
        if (i) fprintf(f, ",");
        fprintf(f, "%" PRIx64, fl[i]);
    }
    return fl.size();
}


static const int bucket_prime_stats = 0;
static long nr_bucket_primes = 0;
static long nr_bucket_longhints = 0;
static long nr_div_tests = 0;
static long nr_composite_tests = 0;
static long nr_wrap_was_composite = 0;
/* The entries in BP must be sorted in order of increasing x */
static void
divide_primes_from_bucket (factor_list_t & fl, mpz_t norm, const unsigned int N, const unsigned int x,
                           bucket_primes_t *BP, const int very_verbose)
{
  while (!BP->is_end()) {
      const bucket_update_t<1, primehint_t> prime = BP->get_next_update();
      if (prime.x > x)
        {
          BP->rewind_by_1();
          break;
        }
      if (prime.x == x) {
          if (bucket_prime_stats) nr_bucket_primes++;
          const unsigned long p = prime.p;
          if (very_verbose) {
              verbose_output_vfprint(0, 1, gmp_vfprintf,
                                     "# N = %u, x = %d, dividing out prime hint p = %lu, norm = %Zd\n",
                                     N, x, p, norm);
          }
          /* If powers of a prime p get bucket-sieved and more than one such
              power hits, then the second (and later) hints will find a
              cofactor that already had all powers of p divided out by the
              loop below that removes multiplicities. Thus, if a prime does
              not divide, we check whether it was divided out before (and thus
              is in the factors list) */
          if (UNLIKELY(!mpz_divisible_ui_p (norm, p))) {
              if(!factor_list_contains(fl, p)) {
                  verbose_output_print(1, 0,
                           "# Error, p = %lu does not divide at (N,x) = (%u,%d)\n",
                           p, N, x);
                  abort();
              } else {
                  verbose_output_print(0, 2,
                           "# Note (harmless): p = %lu does not divide at (N,x) = (%u,%d), was divided out before\n",
                           p, N, x);
              }
          } else do {
              /* Remove powers of prime divisors */
              factor_list_add(fl, p);
              mpz_divexact_ui (norm, norm, p);
              /* Lacking bucket-sieving for powers, we have to check for
               * divisibility once again */
          } while (mpz_divisible_ui_p (norm, p));
      }
  }
}


/* The entries in BP must be sorted in order of increasing x */
static void
divide_hints_from_bucket (factor_list_t &fl, mpz_t norm, const unsigned int N, const unsigned int x,
                          bucket_array_complete *purged,
                          fb_factorbase::slicing const & fbs,
                          const int very_verbose)
{
  while (!purged->is_end()) {
      const bucket_update_t<1, longhint_t> complete_hint = purged->get_next_update();
      if (complete_hint.x > x)
        {
          purged->rewind_by_1();
          break;
        }
      if (complete_hint.x == x) {
          if (bucket_prime_stats) nr_bucket_longhints++;
          fb_slice_interface const & fb_slice = fbs[complete_hint.index];
          const unsigned long p = fb_slice.get_prime(complete_hint.hint);
          if (very_verbose) {
              const unsigned char k = fb_slice.get_k(complete_hint.hint);
              verbose_output_print(0, 1,
                                   "# N = %u, x = %d, dividing out fb_slice hint, "
                                   "index = %lu offset = %lu ",
                                   N, x, (unsigned long) complete_hint.index,
                                   (unsigned long) complete_hint.hint);
              if (fb_slice.is_general()) {
                  verbose_output_print(0, 1, "(general)");
              } else {
                  verbose_output_print(0, 1, "(%d roots)", fb_slice.get_nr_roots());
              }
              verbose_output_vfprint(0, 1, gmp_vfprintf,
                                     ", q = %lu^%hhu, norm = %Zd\n",
                                     p, k, norm);
          }
          if (UNLIKELY(!mpz_divisible_ui_p (norm, p))) {
              if (!factor_list_contains(fl, p)) {
                  verbose_output_print(1, 0,
                           "# Error, p = %lu does not divide at (N,x) = (%u,%d)\n",
                           p, N, x);
                  abort();
              } else {
                  verbose_output_print(0, 2,
                           "# Note (harmless): p = %lu (from hint) does not divide at (N,x) = (%u,%d), was divided out before\n",
                           p, N, x);
              }
          } else do {
              factor_list_add(fl, p);
              mpz_divexact_ui (norm, norm, p);
          } while (mpz_divisible_ui_p (norm, p));
      }
  }
}

/* Extract all known primes (from bucket and small sieves) from the norm.
 * It also removes all the tiny factors that were not resieved and are
 * therefore trial-divided. (see -ththresh parameter)
 *
 * Note: there is another function trialdiv() without underscore that
 * does just the second step.
 */
NOPROFILE_STATIC void
divide_known_primes (std::vector<uint64_t> & fl, cxx_mpz & norm, const unsigned int N, unsigned int x,
           const bool handle_2, bucket_primes_t *primes,
           bucket_array_complete *purged,
	   trialdiv_data const & td,
           int64_t a, uint64_t b,
           fb_factorbase::slicing const & fbs)
{
#ifdef TRACE_K
    const int trial_div_very_verbose = trace_on_spot_ab(a,b);
#else
    const int trial_div_very_verbose = 0;
#endif

    if (trial_div_very_verbose) {
        verbose_output_start_batch();
        verbose_output_print(TRACE_CHANNEL, 0, "# divide_known_primes() entry, N = %u, x = %d, a = %" PRId64 ", b = %" PRIu64 ", norm = ", N, x, a, b);
        verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf, "%Zd\n", (mpz_srcptr) norm);
    }

    // handle 2 separately, if it is in fb
    if (handle_2) {
        int bit = mpz_scan1(norm, 0);
        for (int i = 0; i < bit; ++i)
            fl.push_back(2);
        if (trial_div_very_verbose)
            verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf, "# x = %d, dividing out 2^%d, norm = %Zd\n", x, bit, (mpz_srcptr) norm);
        mpz_tdiv_q_2exp(norm, norm, bit);
    }

    // remove primes in "primes" that map to x
    divide_primes_from_bucket (fl, norm, N, x, primes, trial_div_very_verbose);
    size_t nf_divide_primes = fl.size();

    // now remove prime hints in "purged". If we had no factor base, then
    // we really should have an empty list here.
    divide_hints_from_bucket (fl, norm, N, x, purged, fbs, trial_div_very_verbose);
    size_t nf_divide_hints = fl.size();

    if (trial_div_very_verbose)
        verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf, "# x = %d, after dividing out bucket/resieved norm = %Zd\n", x, (mpz_srcptr) norm);

    /* Trial divide primes with precomputed tables */

    if (trial_div_very_verbose) {
        std::ostringstream os;
        for(auto p : td) os << " " << p.p;
        verbose_output_print(TRACE_CHANNEL, 0, "# Trial division by%s\n",
                os.str().c_str());
    }

    td.trial_divide(fl, norm);
    size_t nf_td = fl.size();

    if (trial_div_very_verbose) {
        std::ostringstream os;
        size_t i = 0;
        if (i < nf_divide_primes) {
            os << " [resieved:";
            for( ; i < nf_divide_primes ; ++i) os << " " << fl[i];
            os << "]";
        }
        if (i < nf_divide_hints) {
            os << " [hints:";
            for( ; i < nf_divide_hints ; ++i) os << " " << fl[i];
            os << "]";
        }
        if (i < nf_td) {
            os << " [trialdiv:";
            for( ; i < nf_td ; ++i) os << " " << fl[i];
            os << "]";
        }
        verbose_output_print(TRACE_CHANNEL, 0, "# %zu factors found:%s\n",
                fl.size(), os.str().c_str());
        verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf, "# After trialdiv(): norm = %Zd\n", (mpz_srcptr) norm);
    }

    if (trial_div_very_verbose)
        verbose_output_end_batch();
}
/* }}} */

MAYBE_UNUSED static inline void subusb(unsigned char *S1, unsigned char *S2, ssize_t offset)
{
    int ex = (unsigned int) S1[offset] - (unsigned int) S2[offset];
    if (UNLIKELY(ex < 0)) S1[offset] = 0; else S1[offset] = ex;	     
}

/* S1 = S1 - S2, with "-" in saturated arithmetic,
 * and memset(S2, 0, EndS1-S1).
 */
void SminusS (unsigned char *S1, unsigned char *EndS1, unsigned char *S2) {/*{{{*/
#ifndef HAVE_SSE2
    ssize_t mysize = EndS1 - S1;
    unsigned char *cS2 = S2;
    while (S1 < EndS1) {
        subusb(S1,S2,0);
        subusb(S1,S2,1);
        subusb(S1,S2,2);
        subusb(S1,S2,3);
        subusb(S1,S2,4);
        subusb(S1,S2,5);
        subusb(S1,S2,6);
        subusb(S1,S2,7);
        S1 += 8; S2 += 8;
    }
    memset(cS2, 0, mysize);
#else
    __m128i *S1i = (__m128i *) S1, *EndS1i = (__m128i *) EndS1, *S2i = (__m128i *) S2,
            z = _mm_setzero_si128();
    while (S1i < EndS1i) {
        __m128i x0, x1, x2, x3;
        __asm__ __volatile__
            ("prefetcht0 0x1000(%0)\n"
             "prefetcht0 0x1000(%1)\n"
             "movdqa (%0),%2\n"
             "movdqa 0x10(%0),%3\n"
             "movdqa 0x20(%0),%4\n"
             "movdqa 0x30(%0),%5\n"
             "psubusb (%1),%2\n"
             "psubusb 0x10(%1),%3\n"
             "psubusb 0x20(%1),%4\n"
             "psubusb 0x30(%1),%5\n"
             "movdqa %6,(%1)\n"
             "movdqa %6,0x10(%1)\n"
             "movdqa %6,0x20(%1)\n"
             "movdqa %6,0x30(%1)\n"
             "movdqa %2,(%0)\n"
             "movdqa %3,0x10(%0)\n"
             "movdqa %4,0x20(%0)\n"
             "movdqa %5,0x30(%0)\n"
             "add $0x40,%0\n"
             "add $0x40,%1\n"
             : "+&r"(S1i), "+&r"(S2i), "=&x"(x0), "=&x"(x1), "=&x"(x2), "=&x"(x3) : "x"(z));
        /* I prefer use ASM than intrinsics to be sure each 4
         * instructions which use exactly a cache line are together. I'm
         * 99% sure it's not useful...  but it's more beautiful :-)
         */
        /*
           __m128i x0, x1, x2, x3;
           _mm_prefetch(S1i + 16, _MM_HINT_T0); _mm_prefetch(S2i + 16, _MM_HINT_T0);
           x0 = _mm_load_si128(S1i + 0);         x1 = _mm_load_si128(S1i + 1);
           x2 = _mm_load_si128(S1i + 2);         x3 = _mm_load_si128(S1i + 3);
           x0 = _mm_subs_epu8(S2i[0], x0);       x1 = _mm_subs_epu8(S2i[1], x1);
           x2 = _mm_subs_epu8(S2i[2], x2);       x3 = _mm_subs_epu8(S2i[3], x3);
           _mm_store_si128(S2i + 0, z);          _mm_store_si128(S1i + 1, z);
           _mm_store_si128(S2i + 2, z);          _mm_store_si128(S1i + 3, z);
           _mm_store_si128(S1i + 0, x0);         _mm_store_si128(S1i + 1, x1);
           _mm_store_si128(S1i + 2, x2);         _mm_store_si128(S1i + 3, x3);
           S1i += 4; S2i += 4;
           */
    }
#endif 
}/*}}}*/

#ifdef  DLP_DESCENT
/* This returns true only if this descent node is now done, either based
 * on the new relation we have registered, or because the previous
 * relation is better anyway */
bool register_contending_relation(las_info const & las, las_todo_entry const & doing, relation & rel)/*{{{*/
{
    if (las.tree.must_avoid(rel)) {
        verbose_output_vfprint(0, 1, gmp_vfprintf, "# [descent] Warning: we have already used this relation, avoiding\n");
        return true;
    }

    /* compute rho for all primes, even on the rational side */
    rel.fixup_r(true);

    descent_tree::candidate_relation contender;
    contender.rel = rel;
    double time_left = 0;

    for(int side = 0 ; side < 2 ; side++) {
        for(unsigned int i = 0 ; i < rel.sides[side].size() ; i++) {
            relation::pr const & v(rel.sides[side][i]);
            if (mpz_cmp(doing.p, v.p) == 0)
                continue;
            unsigned long p = mpz_get_ui(v.p);
            if (mpz_fits_ulong_p(v.p)) {
                unsigned long r = mpz_get_ui(v.r);
                if (las.dlog_base.is_known(side, p, r))
                    continue;
            }

            unsigned int n = mpz_sizeinbase(v.p, 2);
            siever_config_pool::key_type K(side, n);
            double e = las.config_pool.hint_expected_time(K);
            if (e < 0) {
                /* This is not worrysome per se. We just do
                 * not have the info in the descent hint table,
                 * period.
                 */
                verbose_output_vfprint(0, 1, gmp_vfprintf, "# [descent] Warning: cannot estimate refactoring time for relation involving %d@%d (%Zd,%Zd)\n", n, side, (mpz_srcptr) v.p, (mpz_srcptr) v.r);
                time_left = INFINITY;
            } else {
                if (std::isfinite(time_left))
                    time_left += e;
            }
            contender.outstanding.push_back(std::make_pair(side, v));
        }
    }
    verbose_output_print(0, 1, "# [descent] This relation entails an additional time of %.2f for the smoothing process (%zu children)\n",
            time_left, contender.outstanding.size());

    /* when we're re-examining this special-q because of a previous
     * failure, there's absolutely no reason to hurry up on a relation */
    contender.set_time_left(time_left, doing.iteration ? INFINITY : general_grace_time_ratio);

    return las.tree.new_candidate_relation(contender);
}/*}}}*/
#endif /* DLP_DESCENT */

struct process_bucket_region_run : public process_bucket_region_spawn {
    worker_thread * worker;
    nfs_aux::thread_data & taux;
    nfs_work::thread_data & tws;
    timetree_t & timer;
    int bucket_relative_index;
    las_report& rep;
    unsigned char * S[2];
    /* We will have this point to the thread's where_am_I data member.
     * (within nfs_aux::th). However it might be just as easy to let this
     * field be defined here, and drop the latter.
     */
    where_am_I & w;
    bool do_resieve;

    /* A note on SS versus S[side]
     *
     * SS is temp data. It's only used here, and it could well be defined
     * here only. We declare it at the thread_data level to avoid
     * constant malloc()/free().
     *
     * S[side] is where we compute the norm initialization. Some
     * tolerance is subtracted from these lognorms to account for
     * accepted cofactors.
     *
     * SS is the bucket region where we apply the buckets, and also later
     * where we do the small sieve.
     *
     * as long as SS[x] >= S[side][x], we are good.
     */

    unsigned char *SS;
    
    /* FIXME: Having the "primes" array allocated to BUCKET_REGION seems
     * grossly wrong. See e-mail exchange ET/PG 20170622141358 */
    struct side_data {/*{{{*/
        bucket_array_complete purged;   /* for purge_buckets */
        bucket_primes_t primes;         /* for resieving */
        side_data() :
            purged(bucket_array_complete(BUCKET_REGION)),
            primes(bucket_primes_t(BUCKET_REGION))
        {}
    };/*}}}*/

    std::array<side_data, 2> sides;

    process_bucket_region_run(process_bucket_region_spawn const & p, worker_thread * worker, int id);

    /* will be passed as results of functions
    std::vector<uint32_t> survivors;
    std::vector<bucket_update_t<1, shorthint_t>::br_index_t> survivors2;
     * */

    /* most probably useless, I guess
    int N;
    int cpt;
    int copr;
    */

    void init_norms(int side);

    template<bool with_hints>
    void apply_buckets_inner(int side);

    void apply_buckets(int side);
    void small_sieve(int side);
    void SminusS(int side);
    typedef std::vector<bucket_update_t<1, shorthint_t>::br_index_t> survivors_t;
    survivors_t search_survivors();
    void purge_buckets(int side);
    void resieve(int side);
    void cofactoring_sync (survivors_t & survivors2);
    void operator()();
};


/*{{{ process_bucket_region, split into pieces. */
    process_bucket_region_run::process_bucket_region_run(process_bucket_region_spawn const & p, worker_thread * worker, int id): /* {{{ */
            process_bucket_region_spawn(p),
            worker(worker),
            taux(aux_p->th[worker->rank()]),
            tws(ws.th[worker->rank()]),
            timer(aux_p->get_timer(worker)),
            bucket_relative_index(id),
            rep(taux.rep),
            w(taux.w)
    {
        w = w_saved;
        WHERE_AM_I_UPDATE(w, N, first_region0_index + already_done + bucket_relative_index);

        /* This is local to this thread */
        for(int side = 0 ; side < 2 ; side++) {
            nfs_work::side_data & wss(ws.sides[side]);
            if (wss.no_fb())
                S[side] = NULL;
            else
                S[side] = tws.sides[side].bucket_region;
        }

        SS = tws.SS;
        memset(SS, 0, BUCKET_REGION);

        /* see comment in process_bucket_region_run::operator()() */
        do_resieve = ws.conf.sides[0].lim && ws.conf.sides[1].lim;

        /* we're ready to go ! processing is in the operator() method.
         */
    }/*}}}*/
void process_bucket_region_spawn::operator()(worker_thread * worker, int id) /*{{{{*/
{
    timetree_t & timer(aux_p->get_timer(worker));
    ENTER_THREAD_TIMER(timer);
    /* create a temp object with more fields, and dispose it shortly
     * afterwards once we're done.  */
    process_bucket_region_run(*this, worker, id)();
}/*}}}*/
void process_bucket_region_run::init_norms(int side)/*{{{*/
{
    CHILD_TIMER(timer, "init norms");
    TIMER_CATEGORY(timer, norms(side));

    int N = first_region0_index + already_done + bucket_relative_index;

    rep.tn[side] -= seconds_thread ();
    ws.sides[side].lognorms.fill(S[side], N);
    rep.tn[side] += seconds_thread ();

#if defined(TRACE_K) 
    if (trace_on_spot_N(w.N))
        verbose_output_print(TRACE_CHANNEL, 0, "# After side %d init_norms_bucket_region, N=%u S[%u]=%u\n",
                side, w.N, trace_Nx.x, S[side][trace_Nx.x]);
#endif
}/*}}}*/

template<bool with_hints> void process_bucket_region_run::apply_buckets_inner(int side)/*{{{*/
{
    nfs_work::side_data & wss(ws.sides[side]);

    typedef typename hints_proxy<with_hints>::l my_longhint_t;
    typedef typename hints_proxy<with_hints>::s my_shorthint_t;
    rep.ttbuckets_apply -= seconds_thread();
    {
        CHILD_TIMER(timer, "apply buckets");
        TIMER_CATEGORY(timer, sieving(side));
        for (auto const & BA : wss.bucket_arrays<1, my_shorthint_t>())
            apply_one_bucket(SS, BA, already_done + bucket_relative_index, wss.fbs->get_part(1), w);
    }

    /* Apply downsorted buckets, if necessary. */
    if (ws.toplevel > 1) {
        CHILD_TIMER(timer, "apply downsorted buckets");
        TIMER_CATEGORY(timer, sieving(side));

        for (auto const & BAd : wss.bucket_arrays<1, my_longhint_t>()) {
            // FIXME: the updates could come from part 3 as well,
            // not only part 2.
            ASSERT_ALWAYS(ws.toplevel <= 2);
            apply_one_bucket(SS, BAd, already_done + bucket_relative_index, wss.fbs->get_part(2), w);
        }
    }
    rep.ttbuckets_apply += seconds_thread();
}/*}}}*/
void process_bucket_region_run::apply_buckets(int side)
{
    if (do_resieve) {
        apply_buckets_inner<true>(side);
    } else {
        apply_buckets_inner<false>(side);
    }
}

void process_bucket_region_run::small_sieve(int side)/*{{{*/
{
    CHILD_TIMER(timer, "small sieve");
    TIMER_CATEGORY(timer, sieving(side));

    nfs_work::side_data & wss(ws.sides[side]);

    sieve_small_bucket_region(SS,
            first_region0_index + already_done + bucket_relative_index,
            wss.ssd,
            wss.ssd.ssdpos_many[bucket_relative_index],
            ws.conf.logI, ws.Q.sublat,
            w);
}/*}}}*/
void process_bucket_region_run::SminusS(int side)/*{{{*/
{
    /* compute S[side][x] = max(S[side][x] - SS[x], 0),
     * and clear SS.  */
    CHILD_TIMER(timer, "S minus S (2)");
    TIMER_CATEGORY(timer, cofactoring(side));

    ::SminusS(S[side], S[side] + BUCKET_REGION, SS);
#if defined(TRACE_K) 
    if (trace_on_spot_N(w.N))
        verbose_output_print(TRACE_CHANNEL, 0,
                "# Final value on side %d, N=%u S[%u]=%u\n",
                side, w.N, trace_Nx.x, S[side][trace_Nx.x]);
#endif
}/*}}}*/
process_bucket_region_run::survivors_t process_bucket_region_run::search_survivors() /*{{{*/
{
    typedef std::vector<uint32_t> surv1_t;

    surv1_t temp_sv;

    CHILD_TIMER(timer, __func__);
    TIMER_CATEGORY(timer, search_survivors());

    int N = first_region0_index + already_done + bucket_relative_index;

    /* change N, which is a bucket number, to
     * (i0, i1, j0, j1) */
    int logI = ws.conf.logI;
    /* This bit of code is replicated from las-smallsieve.cpp */
    const unsigned int log_lines_per_region = MAX(0, LOG_BUCKET_REGION - logI);
    const unsigned int log_regions_per_line = MAX(0, logI - LOG_BUCKET_REGION);
    const unsigned int regions_per_line = 1 << log_regions_per_line;           
    const unsigned int region_rank_in_line = N & (regions_per_line - 1);       
    const unsigned int j0 = (N >> log_regions_per_line) << log_lines_per_region;    
    const unsigned int j1 MAYBE_UNUSED = j0 + (1 << log_lines_per_region);    
    const int I = 1 << logI;                                            
    const int i0 = (region_rank_in_line << LOG_BUCKET_REGION) - I/2;          
    const int i1 = i0 + (1 << MIN(LOG_BUCKET_REGION, logI));     

    ASSERT(j1 > j0); /* even when we have a line fragment */


#ifdef TRACE_K /* {{{ */
    if (trace_on_spot_Nx(N, trace_Nx.x)) {
        verbose_output_print(TRACE_CHANNEL, 0,
                "# When entering factor_survivors for bucket %u, "
                "S[0][%u]=%u, S[1][%u]=%u\n",
                trace_Nx.N, trace_Nx.x,
                S[0] ? S[0][trace_Nx.x] : ~0u,
                trace_Nx.x,
                S[1] ? S[1][trace_Nx.x] : ~0u);
        verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf,
                "# Remaining norms which have not been accounted for in sieving: (%Zd, %Zd)\n",
                (mpz_srcptr) traced_norms[0],
                (mpz_srcptr) traced_norms[1]);
    }
#endif  /* }}} */

#ifdef TRACE_K /* {{{ */
    nfs_work::side_data & side0(ws.sides[0]);
    nfs_work::side_data & side1(ws.sides[1]);
    for (int x = 0; x < 1 << LOG_BUCKET_REGION; x++) {
        if (trace_on_spot_Nx(N, x)) {
            verbose_output_print(TRACE_CHANNEL, 0,
                    "# side0.Bound[%u]=%u, side1.Bound[%u]=%u\n",
                    S[0] ? S[0][trace_Nx.x] : ~0u,
                    S[0] ? (S[0][x] <= side0.lognorms.bound ? 0 : side0.lognorms.bound) : ~0u,
                    S[1] ? S[1][trace_Nx.x] : ~0u,
                    S[1] ? (S[1][x] <= side1.lognorms.bound ? 0 : side1.lognorms.bound) : ~0u);
        }
    }
#endif /* }}} */

    rep.survivors.before_sieve += 1U << LOG_BUCKET_REGION;

    temp_sv.reserve(128);

    for (unsigned int j = j0; j < j1; j++)
    {
        int offset = (j-j0) << logI;

        unsigned char * const both_S[2] = {
            S[0] ? S[0] + offset : NULL,
            S[1] ? S[1] + offset : NULL,
        };
        /* TODO FIXME XXX that's weird. How come don't we merge that with
         * the lognorm computation that goes in the ws.sides[side]
         * regions before apply_buckets + small_sieve ?? Could it help
         * save a bit of time in search_survivors_in_line ?
         */
        const unsigned char both_bounds[2] = {
            ws.sides[0].lognorms.bound,
            ws.sides[1].lognorms.bound,
        };
        size_t old_size = temp_sv.size();

        ASSERT(j < ws.J);

        search_survivors_in_line(both_S, both_bounds,
                j,
                i0, i1,
                N,
                *ws.jd,
                ws.conf.unsieve_thresh,
                *ws.us,
                temp_sv,
                ws.Q.sublat);

        /* Survivors written by search_survivors_in_line() have index
         * relative to their j-line. We need to convert to index within
         * the bucket region by adding line offsets.
         */


        /* When several bucket regions are in a line, we have nothing to
         * change. Note in particular that we must not adjust with
         * respect to the starting i -- this info is already encoded with
         * the bucket number N.
         */

        if (!offset) continue;

        for (size_t i_surv = old_size; i_surv < temp_sv.size(); i_surv++)
            temp_sv[i_surv] += offset;
    }

    /* This used to be called convert_survivors */
    return survivors_t(begin(temp_sv), end(temp_sv));
}/*}}}*/
void process_bucket_region_run::purge_buckets(int side)/*{{{*/
{
    nfs_work::side_data & wss(ws.sides[side]);

    SIBLING_TIMER(timer, "purge buckets");
    TIMER_CATEGORY(timer, cofactoring(side));

    unsigned char * Sx = S[0] ? S[0] : S[1];

    for (auto & BA : wss.bucket_arrays<1, shorthint_t>()) {
#if defined(HAVE_SSE2) && defined(SMALLSET_PURGE)
        sides[side].purged.purge(BA, already_done + bucket_relative_index, Sx, survivors);
#else
        sides[side].purged.purge(BA, already_done + bucket_relative_index, Sx);
#endif
    }

    /* Add entries coming from downsorting, if any */
    for (auto const & BAd : wss.bucket_arrays<1, longhint_t>()) {
        sides[side].purged.purge(BAd, already_done + bucket_relative_index, Sx);
    }

    /* Sort the entries to avoid O(n^2) complexity when looking for
       primes during trial division */
    sides[side].purged.sort();
}/*}}}*/
void process_bucket_region_run::resieve(int side)/*{{{*/
{
    nfs_work::side_data & wss(ws.sides[side]);
    SIBLING_TIMER(timer, "resieve");
    TIMER_CATEGORY(timer, sieving(side));

    unsigned char * Sx = S[0] ? S[0] : S[1];

    /* Resieve small primes for this bucket region and store them 
       together with the primes recovered from the bucket updates */
    resieve_small_bucket_region (&sides[side].primes,
            Sx,
            first_region0_index + already_done + bucket_relative_index,
            wss.ssd,
            wss.ssd.ssdpos_many[bucket_relative_index],
            ws.conf.logI, ws.Q.sublat,
            w);

    /* same reason as above */
    sides[side].primes.sort();
}/*}}}*/
/*{{{ asynchronous cofactorization */
/* This is one input to the late cofactoring process (aka ECM). Here, we
 * mean the stuff that is done detached from the rest of the siever
 * stuff: we no longer care about purging buckets and so on, these may
 * safely be used for later work.
 */
struct cofac_standalone {
    std::array<uint8_t, 2> S;
    std::array<cxx_mpz, 2> norm;
    std::array<factor_list_t, 2> factors;
    std::array<std::vector<cxx_mpz>, 2> lps;
    int64_t a;
    uint64_t b;
#ifdef SUPPORT_LARGE_Q
    cxx_mpz az, bz;
#endif
    cofac_standalone() : a(0), b(0) {/*{{{*/
        S[0] = S[1] = 0;
#ifdef SUPPORT_LARGE_Q
        mpz_set_ui(az, 0);
        mpz_set_ui(bz, 0);
#endif
    }/*}}}*/
    cofac_standalone(int N, size_t x, int logI, qlattice_basis const & Q) {/*{{{*/
        S[0] = S[1] = 0;
        NxToAB (a, b, N, x, logI, Q);
#ifdef SUPPORT_LARGE_Q
        NxToABmpz (az, bz, N, x, logI, Q);
#endif
    }/*}}}*/
#ifdef TRACE_K
    bool trace_on_spot() const {/*{{{*/
        return trace_on_spot_ab(a, b);
    }/*}}}*/
#endif
    inline bool both_even() const {/*{{{*/
#ifndef SUPPORT_LARGE_Q
        return ((((a | b) & 1) == 0));
#else
        return ((mpz_even_p(az) && mpz_even_p(bz)));
#endif
    }/*}}}*/
    bool gcd_coprime_with_q(las_todo_entry const & E) {/*{{{*/
        /* Since the q-lattice is exactly those (a, b) with
           a == rho*b (mod q), q|b  ==>  q|a  ==>  q | gcd(a,b) */
        /* In case of composite sq, have to check all factors... */
        /* FIXME: fast divisibility test here! */
        if (E.is_prime()) {
#ifndef SUPPORT_LARGE_Q
            if (b == 0 || (mpz_cmp_ui(E.p, b) <= 0 && b % mpz_get_ui(E.p) == 0))
#else
            if ((mpz_cmp_ui(bz, 0) == 0) || 
                (mpz_cmp(E.p, bz) <= 0 &&
                 mpz_divisible_p(bz, E.p)))
#endif
                return false;
        } else {
#ifdef SUPPORT_LARGE_Q
            if (mpz_cmp_ui(bz, 0) == 0)
                return false;
            for (auto const& facq : E.prime_factors) {
                if ((mpz_cmp_ui(bz, facq) >= 0) && (mpz_divisible_ui_p(bz, facq))) {
                    return false;
                }
            }
#else
            if (b == 0)
                return false;
            for (auto const& facq : E.prime_factors) {
                if (facq <= b && b % facq == 0) {
                    return false;
                }
            }
#endif
        }
        return true;
    }/*}}}*/
    bool ab_coprime() const {/*{{{*/
#ifndef SUPPORT_LARGE_Q
        return bin_gcd_int64_safe(a,b) == 1;
#else
        cxx_mpz g;
        mpz_gcd(g, az, bz);
        return mpz_cmp_ui(g, 1) == 0;
#endif
    }/*}}}*/
    void print_as_survivor(FILE * f) {/*{{{*/
#ifndef SUPPORT_LARGE_Q
        gmp_fprintf(f, "%" PRId64 " %" PRIu64 " %Zd %Zd\n", a, b,
                (mpz_srcptr) norm[0],
                (mpz_srcptr) norm[1]);
#else
        gmp_fprintf(f, "%Zd %Zd %Zd %Zd\n",
                (mpz_srcptr) az,
                (mpz_srcptr) bz,
                (mpz_srcptr) norm[0],
                (mpz_srcptr) norm[1]);
#endif
    }/*}}}*/
    relation get_relation(las_todo_entry const & doing) {/*{{{*/
#ifndef SUPPORT_LARGE_Q
        relation rel(a, b);
#else
        relation rel(az, bz);
#endif

        /* Note that we explicitly do not bother about storing r in
         * the relations below */
        for (int side = 0; side < 2; side++) {
            for (auto const& z : factors[side])
                rel.add(side, z, 0);
            for (auto const& z : lps[side])
                rel.add(side, z, 0);
        }
        if (doing.is_prime()) {
            rel.add(doing.side, doing.p, 0);
        } else {
            for (auto const& facq : doing.prime_factors)
                rel.add(doing.side, facq, 0);
        }

        rel.compress();
        return rel;
    }/*}}}*/
    void transfer_to_cofac_list(lock_guarded_container<cofac_list> & L, las_todo_entry const & doing) {/*{{{*/
        std::lock_guard<std::mutex> foo(L.mutex());
        /* "doing" must be an object that lives in the main todo list,
         * and will stay alive until the end of the program. Yes, it's a
         * bit dangerous. */
        L.emplace_back(a, b, norm, &doing);
#if 0
        /* make sure threads don't write the cofactor list at the
         * same time !!! */
#warning "possible contention"
        static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
        pthread_mutex_lock(&lock);
        cofac_list_add (L, a, b, norm, doing.side, doing.p);
        pthread_mutex_unlock(&lock);
#endif
    }/*}}}*/
    int factor_both_leftover_norms(nfs_work_cofac & wc) {/*{{{*/
        /* This proxies to las-cofactor.cpp */
        return ::factor_both_leftover_norms(norm,
                lps,
                {{ wc.sc.sides[0].lim, wc.sc.sides[1].lim }},
                wc.strategies);
    }/*}}}*/
};

struct detached_cofac_parameters : public cofac_standalone, public task_parameters {
    std::shared_ptr<nfs_work_cofac> wc_p;
    std::shared_ptr<nfs_aux> aux_p;
    detached_cofac_parameters(std::shared_ptr<nfs_work_cofac> wc_p, std::shared_ptr<nfs_aux> aux_p, cofac_standalone&& C) : cofac_standalone(C), wc_p(wc_p), aux_p(aux_p) {}
};

struct detached_cofac_result : public task_result {
    std::shared_ptr<relation> rel_p;
};

task_result * detached_cofac(worker_thread * worker, task_parameters * _param, int) /* {{{ */
{
    auto clean_param = call_dtor([_param]() { delete _param; });

    /* We must exit by cleaning the param structure we've been given. But
     * everything we do with the objects whose life is dependent on our
     * param structure must of course be completed at this point. This
     * holds as well for the timer. Yet ACTIVATE_TIMER below registers
     * some stuff to be done at dtor time, so it's important that we
     * clean up the parameters *after* the timer cleans up.
     */
    detached_cofac_parameters *param = static_cast<detached_cofac_parameters *>(_param);

    /* Import some contextual stuff. */
    int id = worker->rank();
    nfs_work_cofac & wc(*param->wc_p);
    nfs_aux & aux(*param->aux_p);
    nfs_aux::thread_data & taux(aux.th[id]);
    las_info const & las(wc.las);
    las_report & rep(taux.rep);
    timetree_t & timer(aux.get_timer(worker));
    /* The timer is normally not running, as we're in a thread task.
     * However, in descent mode, this is called synchronously, and then
     * the situation is different since the timer has already been
     * activated above.
     */
#ifndef DLP_DESCENT
    ENTER_THREAD_TIMER(timer);
#else
    CHILD_TIMER(timer, __func__);
#endif
    nfs_aux::rel_hash_t& rel_hash(aux.get_rel_hash());

    cofac_standalone & cur(*param);

    std::array<int, 2> cof_bitsize {{ 0,0 }}; /* placate compiler */
    las.cofac_stats.call(cur.norm, cof_bitsize);

    SIBLING_TIMER(timer, "cofactoring"); // aka factor_both_leftover_norms
    TIMER_CATEGORY(timer, cofactoring_mixed());

    rep.ttcof -= seconds_thread ();
    int pass = cur.factor_both_leftover_norms(wc);
    rep.survivors.cofactored += (pass != 0);
    rep.ttcof += seconds_thread ();

    auto res = new detached_cofac_result;

#ifdef TRACE_K
    if (cur.trace_on_spot() && pass == 0) {
        verbose_output_print(TRACE_CHANNEL, 0,
                "# factor_both_leftover_norm failed for (%" PRId64 ",%" PRIu64 "), ", cur.a, cur.b);
        verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf,
                "remains %Zd, %Zd unfactored\n",
                (mpz_srcptr) cur.norm[0],
                (mpz_srcptr) cur.norm[1]);
    }
#endif
    if (pass <= 0) {
        /* a factor was > 2^lpb, or some
           factorization was incomplete */
        return res;
    }

    rep.survivors.smooth++;

    /* yippee: we found a relation! */
    SIBLING_TIMER(timer, "print relations");
    TIMER_CATEGORY(timer, bookkeeping());

    las.cofac_stats.success(cof_bitsize);

    relation rel = cur.get_relation(aux.doing);

#ifdef TRACE_K
    if (cur.trace_on_spot()) {
        verbose_output_print(TRACE_CHANNEL, 0, "# Relation for (%"
                PRId64 ",%" PRIu64 ") printed\n", cur.a, cur.b);
    }
#endif

    {
        int do_check = las.suppress_duplicates;

        /* note that if we have large primes which don't fit in
         * an unsigned long, then the duplicate check will
         * quickly return "no".
         */

        nfs_aux::abpair_t ab(cur.a, cur.b);
        bool is_new_rel;
        {
            std::lock_guard<std::mutex> foo(rel_hash.mutex());
            is_new_rel = rel_hash.insert(ab).second;
        }

        const char * dup_comment = NULL;

        if (do_check && relation_is_duplicate(rel, wc.doing, wc.las)) {
            dup_comment = "# DUPE ";
        } else {
            if (!is_new_rel) {
                /* the relation was already printed in a prior attempt,
                 * that was aborted because of an exception. */
                dup_comment = "# DUP ";
            }
            /* Even if the relation was already printed, the las_report
             * object is (now) specific to the current attempt, and has
             * no memory of the number of reports for the failed
             * attempt. Therefore we need to count the relation as a
             * report no matter what.
             */
            rep.reports ++;
            /* Not clear what gives when we have Galois relations.  */
        }

        /* In some rare cases, the norm on one side is exactly 1, which
         * creates undefined behaviour later on. (bug # 21707) */
        if (rel.nb_polys > 2)
            for (int i = 0; i < rel.nb_polys; ++i)
                if (rel.sides[i].size() == 0)
                    dup_comment = "# NORM1 ";

        if (!dup_comment) dup_comment = "";

        std::ostringstream os;

        if (prepend_relation_time)
            os << "(" << seconds() - tt_qstart << ") ";

        // verbose_output_print(0, 3, "# i=%d, j=%u, lognorms = %hhu, %hhu\n", i, j, cur.S[0], cur.S[1]);

        os << dup_comment << rel << "\n";

        if(las.galois != NULL) {
            // adding relations on the fly in Galois cases
            // once filtering is ok for all Galois cases, 
            // this entire block would have to disappear
            add_relations_with_galois(las.galois, os, dup_comment,
                    &rep.reports, rel);
        }

        /* print all in one go */
        verbose_output_start_batch();     /* unlock I/O */
        verbose_output_print(0, 1, "%s", os.str().c_str());
        verbose_output_end_batch();     /* unlock I/O */
#ifdef DLP_DESCENT
        res->rel_p = std::make_shared<relation>(std::move(rel));
#endif
    }

    /* Build histogram of lucky S[x] values */
    rep.mark_report(cur.S[0], cur.S[1]);

    return (task_result*) res;
}

/* }}} */
/*}}}*/
void process_bucket_region_run::cofactoring_sync (survivors_t & survivors)/*{{{*/
{
    CHILD_TIMER(timer, __func__);
    TIMER_CATEGORY(timer, cofactoring_mixed());

    int N = first_region0_index + already_done + bucket_relative_index;
    unsigned char * Sx = S[0] ? S[0] : S[1];

    cofac_standalone cur;


    for (size_t i_surv = 0 ; i_surv < survivors.size(); i_surv++) {
#ifdef DLP_DESCENT
        if (ws.las.tree.must_take_decision())
            break;
#endif
        const size_t x = survivors[i_surv];
        ASSERT_ALWAYS (Sx[x] != 255);
        ASSERT(x < ((size_t) 1 << LOG_BUCKET_REGION));

        rep.survivors.after_sieve++;

        if (S[0] && S[1])
            rep.mark_survivor(S[0][x], S[1][x]);

        /* For factor_leftover_norm, we need to pass the information of the
         * sieve bound. If a cofactor is less than the square of the sieve
         * bound, it is necessarily prime. we implement this by keeping the
         * log to base 2 of the sieve limits on each side, and compare the
         * bitsize of the cofactor with their double.
         */

        SIBLING_TIMER(timer, "check_coprime");

        /* start building a new object. This is a swap operation */
        cur = cofac_standalone(N, x, ws.conf.logI, ws.Q);

        for(int side = 0 ; side < 2 ; side++) {
            if (ws.sides[side].no_fb()) continue;
            cur.S[side] = S[side][x];
        }

#ifdef TRACE_K/*{{{*/
        if (cur.trace_on_spot())
            verbose_output_print(TRACE_CHANNEL, 0, "# about to start cofactorization for (%"
                    PRId64 ",%" PRIu64 ")  %zu %u\n", cur.a, cur.b, x, Sx[x]);
#endif/*}}}*/

        /* since a,b both even were not sieved, either a or b should
         * be odd. However, exceptionally small norms, even without
         * sieving, may fall below the report bound (see tracker
         * issue #15437). Therefore it is safe to continue here. */
        // ASSERT((a | b) & 1);
        if (UNLIKELY(cur.both_even()))
            continue;

        rep.survivors.not_both_even++;

        if (!cur.gcd_coprime_with_q(ws.Q.doing))
            continue;

        rep.survivors.not_both_multiples_of_p++;

        BOOKKEEPING_TIMER(timer);
        int pass = 1;

        int i;
        unsigned int j;
        // Note that are, (i,j) must be true coordinates, not the
        // ones reduced to (-I/2, I/2) using sublattices.
        NxToIJ (i, j, N, x, ws.conf.logI);
        adjustIJsublat(i, j, ws.Q.sublat);

        if (do_resieve) {

            for(int pside = 0 ; pass && pside < 2 ; pside++) {
                int side = trialdiv_first_side ^ pside;
                nfs_work::side_data & wss(ws.sides[side]);

                CHILD_TIMER_PARAMETRIC(timer, "side ", side, " pre-cofactoring checks");
                TIMER_CATEGORY(timer, cofactoring(side));


                SIBLING_TIMER(timer, "recompute complete norm");

                // Trial divide norm on side 'side'
                /* Compute the norms using the polynomials transformed to 
                   i,j-coordinates. The transformed polynomial on the 
                   special-q side is already divided by q */
                wss.lognorms.norm(cur.norm[side], i, j);

#ifdef TRACE_K/*{{{*/
                if (cur.trace_on_spot()) {
                    verbose_output_vfprint(TRACE_CHANNEL, 0,
                            gmp_vfprintf, "# start trial division for norm=%Zd ", (mpz_srcptr) cur.norm[side]);
                    verbose_output_print(TRACE_CHANNEL, 0,
                            "on side %d for (%" PRId64 ",%" PRIu64 ")\n", side, cur.a, cur.b);
                }
#endif/*}}}*/

                if (wss.no_fb()) {
                    /* This is a shortcut. We're probably replacing sieving
                     * by a product tree, there's no reason to bother doing
                     * trial division at this point (or maybe there is ?
                     * would that change the bit size significantly ?) */
                    rep.survivors.check_leftover_norm_on_side[side] ++;
                    continue;
                }

                SIBLING_TIMER(timer, "trial division");

                verbose_output_print(1, 2, "FIXME %s, line %d\n", __FILE__, __LINE__);
                const bool handle_2 = true; /* FIXME */
                rep.survivors.trial_divided_on_side[side]++;

                divide_known_primes (cur.factors[side], cur.norm[side], N, x,
                        handle_2,
                        &sides[side].primes,
                        &sides[side].purged,
                        *wss.td,
                        cur.a, cur.b,
                        *wss.fbs);

                /* if q is composite, its prime factors have not been sieved.
                 * Check if they divide. They probably don't, since we
                 * have computed the norm with the polynomials adapted to
                 * the (i,j) plane, and q divided out. But still,
                 * valuations are not desired here.
                 */
                if ((side == ws.Q.doing.side) && (!ws.Q.doing.is_prime())) {
                    for (const auto &x : ws.Q.doing.prime_factors) {
                        if (mpz_divisible_uint64_p(cur.norm[side], x)) {
                            mpz_divexact_uint64(cur.norm[side], cur.norm[side], x);
                            cur.factors[side].push_back(x);
                        }
                    }
                }

                SIBLING_TIMER(timer, "check_leftover_norm");

                pass = check_leftover_norm (cur.norm[side], ws.conf.sides[side]);
#ifdef TRACE_K
                if (cur.trace_on_spot()) {
                    verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf,
                            "# checked leftover norm=%Zd", (mpz_srcptr) cur.norm[side]);
                    verbose_output_print(TRACE_CHANNEL, 0,
                            " on side %d for (%" PRId64 ",%" PRIu64 "): %d\n",
                            side, cur.a, cur.b, pass);
                }
#endif
                rep.survivors.check_leftover_norm_on_side[side] += pass;
            }
        } else {
            ASSERT_ALWAYS(ws.las.batch || ws.las.batch_print_survivors);

            /* no resieve, so no list of prime factors to divide. No
             * point in doing trial division anyway either.
             */
            for(int side = 0 ; side < 2 ; side++) {
                CHILD_TIMER_PARAMETRIC(timer, "side ", side, " pre-cofactoring checks");
                TIMER_CATEGORY(timer, cofactoring(side));

                SIBLING_TIMER(timer, "recompute complete norm");

                nfs_work::side_data & wss(ws.sides[side]);

                /* factor() in batch.cpp recomputes the complete norm, so
                 * there's no need to compute the norm right now for the
                 * side we've sieved with.
                 */
                if (wss.no_fb()) {
                    wss.lognorms.norm(cur.norm[side], i, j);
                } else {
                    /* This is recognized specially in the
                     * factor_simple_minded() code in batch.cpp
                     */
                    mpz_set_ui(cur.norm[side], 0);
                }

                /* We don't even bother with q and its prime factors.
                 * We're expecting to recover just everything after the
                 * game anyway */

                /* Note that we're *NOT* doing the equivalent of
                 * check_leftover_norm here. This is explained by two
                 * things:
                 *
                 *  - while the "red zone" of post-sieve values that we
                 *  know can't yield relations is quite wide (from L to
                 *  B^2), it's only a marginal fraction of the total
                 *  number of reports. Even more so if we take into
                 *  account the necessary tolerance near the boundaries
                 *  of the red zone.
                 *
                 *  - we don't have the complete norm (with factors taken
                 *  out) at this point, so there's no way we can do a
                 *  primality check -- which is, in fact, the most
                 *  stringent check because it applies to the bulk of the
                 *  candidates.
                 *
                 * Bottom line: we just hand over *everything* to the
                 * batch cofactorization.
                 */
            }
        }

        if (!pass) continue;

#if 0
        if (ws.las.batch || ws.las.batch_print_survivors) {
            /* in these cases, we won't go through detached_cofac, hence
             * we won't check for potential (a,b) duplicates. Those can
             * and will happen if we encounter overflowing buckets at
             * level 2.
             *
             * However, in the long term we expect that this won't happen
             * for real, as overflowing buckets are unlikely to occur
             * unless very early in the process. Therefore we can omit
             * the check here, and do it only at relation printing time.
             */
            nfs_aux::rel_hash_t& rel_hash(aux_p->get_rel_hash());
            nfs_aux::abpair_t ab(cur.a, cur.b);
            std::lock_guard<std::mutex> foo(rel_hash.mutex());
            if (!rel_hash.insert(ab).second) continue;
        }
#endif

        rep.survivors.enter_cofactoring++;

        if (ws.las.batch_print_survivors) {
            // In sublat mode, some non-primitive survivors can exist.
            // The cofactoring via ECM is made aware of this, but not the
            // batch mode, so we have to ensure it.
            rep.reports++;
            if (ws.conf.sublat_bound && !cur.ab_coprime()) continue;
            verbose_output_start_batch ();
            cur.print_as_survivor(ws.las.batch_print_survivors);
            verbose_output_end_batch ();
            continue;
        }

        if (ws.las.batch)
        {
            /* see above */
            rep.reports++;
            if (ws.conf.sublat_bound && !cur.ab_coprime()) continue;
            /* make sure threads don't write the cofactor list at the
             * same time !!! */
            cur.transfer_to_cofac_list(ws.cofac_candidates, aux_p->doing);
            continue; /* we deal with all cofactors at the end of las */
        }

        auto D = new detached_cofac_parameters(wc_p, aux_p, std::move(cur));

#ifndef  DLP_DESCENT
        worker->get_pool().add_task(detached_cofac, D, N, 1); /* id N, queue 1 */
#else
        /* We must proceed synchronously for the descent */
        auto res = dynamic_cast<detached_cofac_result*>(detached_cofac(worker, D, N));
        bool cc = false;
        if (res->rel_p) {
            cc = register_contending_relation(ws.las, ws.Q.doing, *res->rel_p);
        }
        delete res;
        if (cc)
            break;
#endif  /* DLP_DESCENT */
    }
}/*}}}*/
void process_bucket_region_run::operator()() {/*{{{*/

    // This is too verbose.
    // fprintf(stderr, "=== entering PBR for report id %lu\n", rep.id);

    /* first check some early abort conditions. */
    if (recursive_descent) {
        /* For the descent mode, we bail out as early as possible. We
         * need to do so in a multithread-compatible way, though.
         * Therefore the following access is mutex-protected within
         * las.tree. */
        if (ws.las.tree.must_take_decision())
            return;
    } else if (exit_after_rel_found) {
        if (rep.reports)
            return;
    }

    for (int side = 0; side < 2; side++) {
        WHERE_AM_I_UPDATE(w, side, side);
        nfs_work::side_data & wss(ws.sides[side]);
        if (wss.no_fb()) {
            ASSERT_ALWAYS(S[side] == NULL);
            continue;
        }

        MARK_TIMER_FOR_SIDE(timer, side);
        TIMER_CATEGORY(timer, sieving(side));

        /* Compute norms in S[side] */
        init_norms(side);

        /* Accumulate sieve contributions in SS */
        apply_buckets(side);
        small_sieve(side);

        /* compute S[side][x] = max(S[side][x] - SS[x], 0),
         * and clear SS.  */
        SminusS(side);

        ws.sides[side].dumpfile.write(S[side], BUCKET_REGION);

        BOOKKEEPING_TIMER(timer);
    }

    if (las_output.verbose >= 2)
        taux.update_checksums(tws);

    /* rep.ttf does not count the asynchronous time spent in
     * detached_cofac. */
    rep.ttf -= seconds_thread ();

    auto survivors = search_survivors();

    /* The "do_resieve" flag (set above in the ctor) checks if one of the
     * factor bases is empty. This means that we may have decided to
     * *not* sieve on that side, refusing to pay a per-area time a second
     * time. This is based on the rationale that we expect the *other*
     * side to be such a selective test that it isn't worth the trouble.
     * But then, it means that purge_buckets and resieving are not worth
     * the trouble either.
     */

    /* These two steps used to be called "prepare_cofactoring" */
    for(int side = 0 ; do_resieve && side < 2 ; side++) {
        MARK_TIMER_FOR_SIDE(timer, side);
        purge_buckets(side);
        resieve(side);
    }

#ifdef TRACE_K
    int N = first_region0_index + already_done + bucket_relative_index;
    if (trace_on_spot_Nx(N, trace_Nx.x)) {
        unsigned char * Sx = S[0] ? S[0] : S[1];
        verbose_output_print(TRACE_CHANNEL, 0, "# Slot [%u] in bucket %u has value %u\n",
                trace_Nx.x, trace_Nx.N, Sx[trace_Nx.x]);
    }
#endif

#ifdef  DLP_DESCENT
    /* we're doing detached_cofac right away in the descent case, so we
     * must untangle the ttf and ttcof timings.
     */
    double tt = rep.ttcof;
#endif

    cofactoring_sync(survivors);

    rep.ttf += seconds_thread ();
#ifdef  DLP_DESCENT
    /* hack, see above. The situation where cofactoring_sync also induces
     * ttcof is now exceptional, and taken care of here.
     */
    rep.ttf -= rep.ttcof - tt;
#endif

}/*}}}*/
/*}}}*/

/*************************** main program ************************************/

static void configure_aliases(cxx_param_list & pl)
{
    las_info::configure_aliases(pl);
    param_list_configure_alias(pl, "log-bucket-region", "B");
    las_augmented_output_channel::configure_aliases(pl);
    tdict::configure_aliases(pl);
}

static void configure_switches(cxx_param_list & pl)
{
    las_info::configure_switches(pl);
    las_augmented_output_channel::configure_switches(pl);
    tdict::configure_switches(pl);

    param_list_configure_switch(pl, "-allow-largesq", &allow_largesq);
#ifdef  DLP_DESCENT
    param_list_configure_switch(pl, "-recursive-descent", &recursive_descent);
#endif

    param_list_configure_switch(pl, "-prepend-relation-time", &prepend_relation_time);
    param_list_configure_switch(pl, "-sync", &sync_at_special_q);
    param_list_configure_switch(pl, "-never-discard", &never_discard);
}

static void declare_usage(cxx_param_list & pl)/*{{{*/
{
    param_list_usage_header(pl,
            "In the names and in the descriptions of the parameters, below there are often\n"
            "aliases corresponding to the convention that 0 is the rational side and 1\n"
            "is the algebraic side. If the two sides are algebraic, then the word\n"
            "'rational' just means the side number 0. Note also that for a rational\n"
            "side, the factor base is recomputed on the fly (or cached), and there is\n"
            "no need to provide a fb0 parameter.\n"
            );

    las_info::declare_usage(pl);
    las_parallel_desc::declare_usage(pl);
    las_todo_list::declare_usage(pl);
    las_augmented_output_channel::declare_usage(pl);
    tdict::declare_usage(pl);

    param_list_decl_usage(pl, "trialdiv-first-side", "begin trial division on this side");
    param_list_decl_usage(pl, "allow-largesq", "allows large special-q, e.g. for a DL descent");

    param_list_decl_usage(pl, "sublat", "modulus for sublattice sieving");

    param_list_decl_usage(pl, "log-bucket-region", "set bucket region to 2^x");

    siever_config::declare_usage(pl);

    param_list_decl_usage(pl, "exit-early", "once a relation has been found, go to next special-q (value==1), or exit (value==2)");
    param_list_decl_usage(pl, "file-cofact", "provide file with strategies for the cofactorization step");
    param_list_decl_usage(pl, "prepend-relation-time", "prefix all relation produced with time offset since beginning of special-q processing");
    param_list_decl_usage(pl, "sync", "synchronize all threads at each special-q");
#ifdef TRACE_K
    param_list_decl_usage(pl, "traceab", "Relation to trace, in a,b format");
    param_list_decl_usage(pl, "traceij", "Relation to trace, in i,j format");
    param_list_decl_usage(pl, "traceNx", "Relation to trace, in N,x format");
#endif
#ifdef  DLP_DESCENT
    param_list_decl_usage(pl, "recursive-descent", "descend primes recursively");
    param_list_decl_usage(pl, "grace-time-ratio", "Fraction of the estimated further descent time which should be spent processing the current special-q, to find a possibly better relation");

#endif /* DLP_DESCENT */
    /* given that this option is dangerous, we used to enable it only for
     * las_descent
     */
    param_list_decl_usage(pl, "never-discard", "Disable the discarding process for special-q's. This is dangerous. See bug #15617");

    verbose_decl_usage(pl);
}/*}}}*/

double nprimes_interval(double p0, double p1)
{
#ifdef HAVE_STDCPP_MATH_SPEC_FUNCS
    return std::expint(log(p1)) - std::expint(log(p0));
#else
    /* that can't be sooo wrong... */
    double l0 = log(p0);
    double l1 = log(p1);
    double s1 = p1*(1/l1+1/pow(l1,2)+2/pow(l1,3)+6/pow(l1,4));
    double s0 = p0*(1/l0+1/pow(l0,2)+2/pow(l0,3)+6/pow(l0,4));
    return s1 - s0;
#endif
}


/* Our fetching of the siever_config fields is definitely wrong here. We
 * should only use logA, logI, and the siever thresholds.
 *
 * Note that we are sowewhat limited in our queries of the las_info type
 * here, as it has not been *completely* initialized yet. Well, almost,
 * but the shared structure cache has not loaded the factor base yet.
 */
size_t expected_memory_usage_per_binding_zone(siever_config const & sc,/*{{{*/
        las_info const & las,
        int print)
{
    int hush = print ? 0 : 3;
    cado_poly_srcptr cpoly = las.cpoly;

    /*
    int logImin = (1+sc.logA)/2;
    if (las.adjust_strategy == 2) {
        logImin -= ADJUST_STRATEGY2_MAX_SQUEEZE;
    sc.instantiate_thresholds(0);
    sc.instantiate_thresholds(1);
    */

    /* do the estimate based on the typical config stuff provided.
     * This is most often going to give a reasonable rough idea anyway.
     */
    size_t memory = 0;

    for(int side = 0 ; side < 2 ; side++) {
        if (!sc.sides[side].lim) continue;
        double p1 = sc.sides[side].lim;
        double p0 = 2;
        /* in theory this should depend on the galois group and so on.
         * Here we're counting only with respect to a full symmetric
         * Galois group.
         * 
         * The average number of roots modulo primes is the average
         * number of fixed points of a permutation, and that is 1. If we
         * average over permutations with at least one fixed point, then
         * we have n! / (n! - D_n), and D_n/n! = \sum_{0\leq k\leq
         * n}(-1)^k/k! (which tends to 1/e, so the average number of
         * roots whenever there's at least one root tends to
         * 1/(1-1/e) = 1.5819767...)
         */
        int d = cpoly->pols[side]->deg;
        double ideals_per_prime = 1;
        double fac=1;
        for(int k = 1 ; k <= d ; k++) {
            fac *= -k;
            ideals_per_prime += 1/fac;
        }
        ideals_per_prime = 1/(1-ideals_per_prime);
        size_t nideals = nprimes_interval(p0, p1);
        size_t nprimes = nideals / ideals_per_prime;
        size_t more = 0;
        /* we have nideals/ideals_per_prime prime numbers, totalling
         * nideals roots.
         * Per prime, we have:
         *      fbprime_t
         *      redc_invp_t
         *      double  (for the weight_cdf table).
         */
        more += sizeof(fbprime_t) * nprimes;
        more += sizeof(redc_invp_t) * nprimes;
        more += sizeof(double) * nprimes;
        /* Per root we have:
         *      fbroot_t
         */
        more += sizeof(fbroot_t) * nideals;

        verbose_output_print(0, 3 + hush,
                "# side %d, lim=%lu, %zu fb primes"
                " (d=%d, %f roots per p if G=S_d): %zuMB\n",
                side,
                sc.sides[side].lim,
                nideals,
                d, ideals_per_prime,
                more >> 20);
        memory += more;

        /* Maybe also count the cost of storing a few slicings ? */
    }
    return memory;
}/*}}}*/
/* This does not count the footprint per binding zone */
size_t expected_memory_usage_per_subjob(siever_config const & sc,/*{{{*/
        las_info const & las,
        int nthreads,
        int print)
{
    char buf[20];
    int hush = print ? 0 : 3;
    bkmult_specifier bkmult = las.get_bk_multiplier();

    /* FIXME: I think that this code misses the case of sublat. */

    /* sc.instantiate_thresholds() depends on sc.logI */
    fb_factorbase::key_type K[2] {
        sc.instantiate_thresholds(0),
        sc.instantiate_thresholds(1)
    };

    bool do_resieve = sc.sides[0].lim && sc.sides[1].lim;

    size_t memory = 0;
    size_t more;

#if 0 && defined(__linux__) && defined(__GLIBC__) && defined(__x86_64__)
    /* count threads. Each costs 8M+4k for the stack, 64MB for the
     * private heap. However, this is only virtual address space.
     * Therefore it's not really clear how we should count it. Maybe "not
     * at all" is a good start, in fact.
     *
     * 64MB is actually transiently 128MB, then 64MB.
     */
    if (0) {
        verbose_output_print(0, 3 + hush, "# %d threads: %zuMB\n",
                nthreads,
                (more = nthreads * 0x4801000) >> 20);
        memory += more;
    }
#endif

    // toplevel is computed by fb_factorbase::slicing::slicing, based on
    // thresholds in fbK
    int toplevel = -1;
    for(int side = 0 ; side < 2 ; side++) {
        int m;
        for(m = 0 ; m < FB_MAX_PARTS && K[side].thresholds[m] < sc.sides[side].lim; ++m);
        if (m > toplevel)
            toplevel = m;
    }

    ASSERT_ALWAYS(toplevel == 1 || toplevel == 2);

    /* the code path is radically different depending on toplevel. */

    int nba = NUMBER_OF_BAS_FOR_THREADS(nthreads);

    double m1s, m1l, m2s;
    size_t s1s, s1l, s2s;

    struct round_me {
        slice_index_t initial;
        slice_index_t increase;
        slice_index_t operator()(slice_index_t y) const {
            return std::max(initial, increase * iceildiv(y, increase));
        }
    };

    round_me round1s, round2s, round1l;

    if (do_resieve) {
        typedef bucket_update_t<1, shorthint_t> T1s;
        typedef bucket_update_t<2, shorthint_t> T2s;
        typedef bucket_update_t<1, longhint_t> T1l;
        typedef bucket_slice_alloc_defaults<1, shorthint_t> W1s;
        typedef bucket_slice_alloc_defaults<2, shorthint_t> W2s;
        typedef bucket_slice_alloc_defaults<1, longhint_t> W1l;
        s2s=sizeof(T2s); m2s=bkmult.get<T2s>();
        s1s=sizeof(T1s); m1s=bkmult.get<T1s>();
        s1l=sizeof(T1l); m1l=bkmult.get<T1l>();
        round1s = round_me { W1s::initial, W1s::increase };
        round2s = round_me { W2s::initial, W2s::increase };
        round1l = round_me { W1l::initial, W1l::increase };
    } else {
        typedef bucket_update_t<1, emptyhint_t> T1s;
        typedef bucket_update_t<2, emptyhint_t> T2s;
        typedef bucket_update_t<1, logphint_t> T1l;
        typedef bucket_slice_alloc_defaults<1, emptyhint_t> W1s;
        typedef bucket_slice_alloc_defaults<2, emptyhint_t> W2s;
        typedef bucket_slice_alloc_defaults<1, logphint_t> W1l;
        s2s=sizeof(T2s); m2s=bkmult.get<T2s>();
        s1s=sizeof(T1s); m1s=bkmult.get<T1s>();
        s1l=sizeof(T1l); m1l=bkmult.get<T1l>();
        round1s = round_me { W1s::initial, W1s::increase };
        round2s = round_me { W2s::initial, W2s::increase };
        round1l = round_me { W1l::initial, W1l::increase };
    }

    if (toplevel == 2) {
        // very large factor base primes, between bkthresh1 and lim = we
        // compute all bucket updates in one go. --> those updates are
        // bucket_update_t<2, shorthint_t>, that is, an XSIZE2 position
        // (should be 24 bits, is actually 32) and a short hint.
        //
        // For each big bucket region (level-2), we then transform these
        // updates into updates for the lower-level buckets. We thus
        // create bucket_update_t<1, longhint_t>'s with the downsort<>
        // function. The long hint is because we have the full fb_slice
        // index. the position in such a bucket update is shorter, only
        // XSIZE1.

        // For moderately large factor base primes (between bkthresh and
        // bkthresh1), we precompute the FK lattices (at some cost), and we
        // fill the buckets locally with short hints (and short positions:
        // bucket_update_t<1, shorthint_t>)

        for(int side = 0 ; side < 2 ; side++) {
            if (!sc.sides[side].lim) continue;
            /* In truth, I sort of know it isn't valid. We've built most
             * of the stuff on the idea that there's a global "toplevel"
             * notion, but that barely applies when one of the factor
             * bases happens to be much smaller than the other one */
            ASSERT_ALWAYS(K[side].thresholds[2] == sc.sides[side].lim);
            double p1 = K[side].thresholds[2];
            double p0 = K[side].thresholds[1];
            p0 = std::min(p1, p0);
            size_t nprimes = nprimes_interval(p0, p1);
            double w = (std::log(std::log(p1)) - std::log(std::log(p0)));

            /* we duplicate code that is found in allocate_memory. TODO:
             * refactor that */
            size_t nreg = 1UL << (sc.logA - LOG_BUCKET_REGIONS[2]);
            size_t nup_per_reg = 0.25 * w * BUCKET_REGIONS[2] / nba;
            /* assume LOG_BUCKET_REGIONS[2] > logI */
            nup_per_reg *= 3;
            nup_per_reg += NB_DEVIATIONS_BUCKET_REGIONS * sqrt(nup_per_reg);
            size_t nupdates = nup_per_reg * nreg * nba;
            {
                verbose_output_print(0, 3 + hush,
                        "# level 2, side %d:"
                        " %zu primes,"
                        " room for %zu 2-updates [2s] in %d arrays:"
                        " %s\n",
                        side, nprimes, nupdates,
                        nba,
                        size_disp(more = m2s * nupdates * s2s, buf));
                memory += more;
            }
            {
                /* Count the slice_start pointers as well. We need to know
                 * how many slices will be processed in each bucket
                 * array. A rough rule of thumb probably works.
                 */
                size_t nslices_estim = iceildiv(nprimes >> 16, nba);
                size_t nslices_alloc = round2s(nslices_estim);
                std::ostringstream os;
                size_t waste = (nslices_alloc - nslices_estim) * nba * nreg * sizeof(void*);
                if (waste > (100<<20))
                    os << " [note: using coarse-grain value of " << nslices_alloc << " slices instead; " << 100.0*(nslices_alloc-nslices_estim)/nslices_alloc << "% waste ("<<(waste>>20)<<" MB) !]";
                verbose_output_print(0, 3 + hush,
                        "# level 2, side %d:"
                        " expect %zu slices per array,"
                        " %zu pointers each, in %d arrays:"
                        " %s%s\n",
                        side,
                        nslices_estim,
                        nreg,
                        nba,
                        size_disp(more = nba * nreg * MAX(nslices_estim, nslices_alloc) * sizeof(void*), buf),
                        os.str().c_str());
                memory += more;
            }
            {
                // how many downsorted updates are alive at a given point in
                // time ?
                size_t nupdates_D = nupdates >> 8;
                verbose_output_print(0, 3 + hush,
                        "# level 1, side %d:"
                        " %zu downsorted 1-updates [1l]: %s\n",
                        side, nupdates_D,
                        size_disp(more = m1l * nupdates_D * s1l, buf));
                memory += more;
            }
            {
                size_t nslices_estim = 1;
                size_t nslices_alloc = round1l(nslices_estim);
                size_t nreg = 1 << (LOG_BUCKET_REGIONS[2] - LOG_BUCKET_REGIONS[1]);
                std::ostringstream os;
                size_t waste = (nslices_alloc - nslices_estim) * nba * nreg * sizeof(void*);
                if (waste > (100<<20))
                    os << " [note: using coarse-grain value of " << nslices_alloc << " slices instead; " << 100.0*(nslices_alloc-nslices_estim)/nslices_alloc << "% waste ("<<(waste>>20)<<" MB) !]";
                verbose_output_print(0, 3 + hush,
                        "# level 1, side %d:"
                        " expect %zu slices per array,"
                        " %zu pointers each, in %d arrays: %s%s\n",
                        side,
                        nslices_estim,
                        nreg,
                        nba,
                        size_disp(more = nba * nreg * MAX(nslices_estim, nslices_alloc) * sizeof(void*), buf),
                        os.str().c_str());
                memory += more;
            }
        }

        for(int side = 0 ; side < 2 ; side++) {
            if (!sc.sides[side].lim) continue;
            double p1 = K[side].thresholds[1];
            double p0 = K[side].thresholds[0];
            size_t nprimes = nprimes_interval(p0, p1);
            double w = (std::log(std::log(p1)) - std::log(std::log(p0)));

            size_t nreg = 1 << (LOG_BUCKET_REGIONS[2] - LOG_BUCKET_REGIONS[1]);
            size_t nup_per_reg = 0.25 * w * BUCKET_REGIONS[1] / nba;
            /* assume LOG_BUCKET_REGIONS[1] > logI -- if it's not the
             * case, the count will not be too wrong anyway. */
            nup_per_reg *= 3;
            nup_per_reg += NB_DEVIATIONS_BUCKET_REGIONS * sqrt(nup_per_reg);
            size_t nupdates = nup_per_reg * nreg * nba;
            verbose_output_print(0, 3 + hush,
                    "# level 1, side %d:"
                    " %zu primes,"
                    " %zu 1-updates [1s] in %d arrays:"
                    " %s\n",
                    side, nprimes, nupdates, nba,
                    size_disp(more = m1s * nupdates * s1s, buf));
            memory += more;
            verbose_output_print(0, 3 + hush,
                    "# level 1, side %d:"
                    " %zu primes => precomp_plattices: %s\n",
                    side, nprimes,
                    size_disp(more = nprimes * sizeof(plattice_enumerator<1>), buf));
            memory += more;

            {
                /* Count the slice_start pointers as well. We need to know
                 * how many slices will be processed in each bucket
                 * array. A rough rule of thumb probably works.
                 */
                size_t nslices_estim = iceildiv(nprimes >> 16, nba);
                size_t nslices_alloc = round1s(nslices_estim);
                std::ostringstream os;
                size_t waste = (nslices_alloc - nslices_estim) * nba * nreg * sizeof(void*);
                if (waste > (100<<20))
                    os << " [note: using coarse-grain value of " << nslices_alloc << " slices instead; " << 100.0*(nslices_alloc-nslices_estim)/nslices_alloc << "% waste ("<<(waste>>20)<<" MB) !]";
                verbose_output_print(0, 3 + hush,
                        "# level 1, side %d:"
                        " expect %zu slices per array,"
                        " %zu pointers each, in %d arrays:"
                        " %s%s\n",
                        side,
                        nslices_estim,
                        nreg,
                        nba,
                        size_disp(more = nba * nreg * MAX(nslices_estim, nslices_alloc) * sizeof(void*), buf),
                        os.str().c_str());
                memory += more;
            }
        }
    } else if (toplevel == 1) {
        // *ALL* bucket updates are computed in one go as
        // bucket_update_t<1, shorthint_t>
        for(int side = 0 ; side < 2 ; side++) {
            if (!sc.sides[side].lim) continue;
            ASSERT_ALWAYS(K[side].thresholds[1] == sc.sides[side].lim);
            double p1 = K[side].thresholds[1];
            double p0 = K[side].thresholds[0];
            size_t nprimes = nprimes_interval(p0, p1);
            double w = (std::log(std::log(p1)) - std::log(std::log(p0)));

            /* we duplicate code that is found in allocate_memory. TODO:
             * refactor that */
            size_t nreg = 1UL << (sc.logA - LOG_BUCKET_REGIONS[1]);
            size_t nup_per_reg = 0.25 * w * BUCKET_REGIONS[1] / nba;
            /* assume LOG_BUCKET_REGIONS[1] > logI -- if it's not the
             * case, the count will not be too wrong anyway. */
            nup_per_reg *= 3;
            nup_per_reg += NB_DEVIATIONS_BUCKET_REGIONS * sqrt(nup_per_reg);
            size_t nupdates = nup_per_reg * nreg * nba;
            nupdates += NB_DEVIATIONS_BUCKET_REGIONS * sqrt(nupdates);
            verbose_output_print(0, 3 + hush,
                    "# level 1, side %d:"
                    " %zu primes,"
                    " %zu 1-updates [1s] in %d arrays:"
                    " %s\n",
                    side, nprimes, nupdates, nba,
                    size_disp(more = m1s * nupdates * s1s, buf));
            memory += more;
            {
                /* Count the slice_start pointers as well. We need to know
                 * how many slices will be processed in each bucket
                 * array. A rough rule of thumb probably works.
                 */
                size_t nslices_estim = iceildiv(nprimes >> 16, nba);
                size_t nslices_alloc = round1s(nslices_estim);
                std::ostringstream os;
                size_t waste = (nslices_alloc - nslices_estim) * nba * nreg * sizeof(void*);
                if (waste > (100<<20))
                    os << " [note: using coarse-grain value of " << nslices_alloc << " slices instead; " << 100.0*(nslices_alloc-nslices_estim)/nslices_alloc << "% waste ("<<(waste>>20)<<" MB) !]";
                verbose_output_print(0, 3 + hush,
                        "# level 1, side %d:"
                        " expect %zu slices per array,"
                        " %zu pointers each, in %d arrays:"
                        " %s%s\n",
                        side,
                        nslices_estim,
                        nreg,
                        nba,
                        size_disp(more = nba * nreg * MAX(nslices_estim, nslices_alloc) * sizeof(void*), buf),
                        os.str().c_str());
                memory += more;
            }
        }
    }

    return memory;
}/*}}}*/
size_t expected_memory_usage_per_subjob_worst_logI(siever_config const & sc0, las_info const & las, int nthreads, int print)/*{{{*/
{
    /* We're not getting number_of_threads_per_subjob() from las, because
     * this function gets called before the las_parallel_desc layer of
     * the las_info structure is set.
     */
    int hush = print ? 0 : 3;
    char buf[20];
    /* do the estimate based on the typical config stuff provided.
     * This is most often going to give a reasonable rough idea anyway.
     */
    siever_config sc = sc0;
    /* See for which I we'll have the most expensive setting */
    int logImin, logImax;
    if (las.adjust_strategy == 2) {
        logImin = (1+sc.logA)/2 - ADJUST_STRATEGY2_MAX_SQUEEZE;
        logImax = (1+sc.logA)/2 - ADJUST_STRATEGY2_MIN_SQUEEZE;
    } else {
        logImin = logImax = (1+sc.logA)/2;
    }
    size_t max_memory = 0;
    int logI_max_memory = 0;
    for(int logI = logImin ; logI <= logImax ; logI++) {
        sc.logI = logI;
        verbose_output_print(0, 3 + hush,
                "# Expected memory usage per subjob for logI=%d [%d threads]:\n",
                sc.logI, nthreads);

        size_t memory = expected_memory_usage_per_subjob(sc, las, nthreads, print);

        verbose_output_print(0, 2 + hush,
                "# Expected memory usage per subjob for logI=%d: %s\n",
                sc.logI, size_disp(memory, buf));

        if (memory > max_memory) {
            logI_max_memory = sc.logI;
            max_memory = memory;
        }
    }
    if (logImin != logImax || las_output.verbose < 2 + hush)
        verbose_output_print(0, 0 + hush,
                "# Expected memory use per subjob (max reached for logI=%d):"
                " %s\n",
                logI_max_memory, size_disp(max_memory, buf));
    return max_memory;
}/*}}}*/

size_t expected_memory_usage(siever_config const & sc,/*{{{*/
        las_info & las,
        int print,
        size_t base_memory = 0)
{
    /* Contrary to the previous functions which are used very early on in
     * order to decide on the parallel setting, this function can safely
     * assume that this parallel setting is now complete.
     */
    int hush = print ? 0 : 3;
    char buf[20];

    size_t fb_memory = expected_memory_usage_per_binding_zone(sc, las, print);
    verbose_output_print(0, 0 + hush,
            "# Expected memory usage per binding zone for the factor base: %s\n",
            size_disp(fb_memory, buf));

    size_t subjob_memory = expected_memory_usage_per_subjob_worst_logI(sc, las, las.number_of_threads_per_subjob(), print);

    size_t memory;
    memory = subjob_memory;
    memory *= las.number_of_subjobs_per_memory_binding_zone();
    memory += fb_memory;
    memory *= las.number_of_memory_binding_zones();
    memory += base_memory;

    verbose_output_print(0, 0 + hush,
            "# Expected memory use for %d binding zones and %d %d-threaded jobs per zone, counting %zu MB of base footprint: %s\n",
            las.number_of_memory_binding_zones(),
            las.number_of_subjobs_per_memory_binding_zone(),
            las.number_of_threads_per_subjob(),     /* per subjob, always */
            base_memory >> 20,
            size_disp(memory, buf));

    return memory;

}/*}}}*/

#ifdef  DLP_DESCENT
void postprocess_specialq_descent(las_info & las, las_todo_list & todo, las_todo_entry const & doing, timetree_t & timer_special_q)/*{{{*/
{
    SIBLING_TIMER(timer_special_q, "descent");
    TIMER_CATEGORY(timer_special_q, bookkeeping());

    descent_tree::candidate_relation const & winner(las.tree.current_best_candidate());
    if (winner) {
        /* Even if not going for recursion, store this as being a
         * winning relation. This is useful for preparing the hint
         * file, and also for the initialization of the descent.
         */
        las.tree.take_decision();
        verbose_output_start_batch();
        FILE * output;
        for (size_t i = 0;
                (output = verbose_output_get(0, 0, i)) != NULL;
                i++) {
            winner.rel.print(output, "Taken: ");
        }
        verbose_output_end_batch();
        {
            unsigned int n = mpz_sizeinbase(doing.p, 2);
            verbose_output_start_batch();
            verbose_output_print (0, 1, "# taking path: ");
            for(int i = 0 ; i < doing.depth ; i++) {
                verbose_output_print (0, 1, " ");
            }
            verbose_output_print (0, 1, "%d@%d ->", n, doing.side);
            for(unsigned int i = 0 ; i < winner.outstanding.size() ; i++) {
                int side = winner.outstanding[i].first;
                relation::pr const & v(winner.outstanding[i].second);
                unsigned int n = mpz_sizeinbase(v.p, 2);
                verbose_output_print (0, 1, " %d@%d", n, side);
            }
            if (winner.outstanding.empty()) {
                verbose_output_print (0, 1, " done");
            }
            verbose_output_vfprint (0, 1, gmp_vfprintf, " \t%d %Zd %Zd\n", doing.side,
                    (mpz_srcptr) doing.p,
                    (mpz_srcptr) doing.r);
            verbose_output_end_batch();
        }
        if (recursive_descent) {
            /* reschedule the possibly still missing large primes in the
             * todo list */
            for(unsigned int i = 0 ; i < winner.outstanding.size() ; i++) {
                int side = winner.outstanding[i].first;
                relation::pr const & v(winner.outstanding[i].second);
                unsigned int n = mpz_sizeinbase(v.p, 2);
                verbose_output_vfprint(0, 1, gmp_vfprintf, "# [descent] " HILIGHT_START "pushing side-%d (%Zd,%Zd) [%d@%d]" HILIGHT_END " to todo list\n", side, (mpz_srcptr) v.p, (mpz_srcptr) v.r, n, side);
                todo.push_withdepth(v.p, v.r, side, doing.depth + 1);
            }
        }
    } else {
        las.tree.mark_try_again(doing.iteration + 1);
        unsigned int n = mpz_sizeinbase(doing.p, 2);
        verbose_output_print (0, 1, "# taking path: %d@%d -> loop (#%d)", n, doing.side, doing.iteration + 1);
        verbose_output_vfprint (0, 1, gmp_vfprintf, " \t%d %Zd %Zd\n", doing.side,
                (mpz_srcptr) doing.p,
                (mpz_srcptr) doing.r);
        verbose_output_vfprint(0, 1, gmp_vfprintf, "# [descent] Failed to find a relation for " HILIGHT_START "side-%d (%Zd,%Zd) [%d@%d]" HILIGHT_END " (iteration %d). Putting back to todo list.\n", doing.side,
                (mpz_srcptr) doing.p,
                (mpz_srcptr) doing.r, n, doing.side, doing.iteration);
        todo.push_withdepth(doing.p, doing.r, doing.side, doing.depth + 1, doing.iteration + 1);
    }
}/*}}}*/
#endif  /* DLP_DESCENT */

void process_many_bucket_regions(nfs_work & ws, std::shared_ptr<nfs_work_cofac> wc_p, std::shared_ptr<nfs_aux> aux_p, thread_pool & pool, int first_region0_index, where_am_I const & w)/*{{{*/
{
    /* first_region0_index is always 0 when toplevel == 1, but the
     * present function is also called from within downsort_tree when
     * toplevel > 1, and then first_region0_index may be larger.
     */
    auto P = thread_pool::make_shared_task<process_bucket_region_spawn>(ws, wc_p, aux_p, w);

    /* Make sure we don't schedule too many tasks when J was truncated
     * anyway */

    int first_skipped_br = ws.J;

    if (ws.conf.logI >= LOG_BUCKET_REGION)
        first_skipped_br <<= ws.conf.logI - LOG_BUCKET_REGION;
    else
        first_skipped_br >>= LOG_BUCKET_REGION - ws.conf.logI;

    size_t small_sieve_regions_ready = std::min(SMALL_SIEVE_START_POSITIONS_MAX_ADVANCE, ws.nb_buckets[1]);

    for(int done = 0, ready = small_sieve_regions_ready ; done < ws.nb_buckets[1] ; ) {

        /* yes, it's a bit ugly */
        P->first_region0_index = first_region0_index;
        P->already_done = done;

        for(int i = 0 ; i < ready ; i++) {
            if (first_region0_index + done + i >= first_skipped_br) {
                /* Hmm, then we should also make sure that we truncated
                 * fill_in_buckets, right ? */
                break;
            }
            pool.add_shared_task(P, i, 0);
        }

        /* it's only really done when we do drain_queue(0), of course */
        done += ready;

        if (done < ws.nb_buckets[1]) {

            /* We need to compute more init positions */
            int more = std::min(SMALL_SIEVE_START_POSITIONS_MAX_ADVANCE, ws.nb_buckets[1] - done);

            for(int side = 0 ; side < 2 ; side++) {
                nfs_work::side_data & wss(ws.sides[side]);
                if (wss.no_fb()) continue;
                pool.add_task_lambda([=,&ws](worker_thread * worker, int){
                        timetree_t & timer(aux_p->get_timer(worker));
                        ENTER_THREAD_TIMER(timer);
                        MARK_TIMER_FOR_SIDE(timer, side);
                        SIBLING_TIMER(timer, "prepare small sieve");
                        nfs_work::side_data & wss(ws.sides[side]);
                        // if (wss.no_fb()) return;
                        SIBLING_TIMER(timer, "small sieve start positions");
                        /* When we're doing 2-level sieving, there is probably
                         * no real point in doing ssdpos initialization in
                         * several passes.
                         */
                        small_sieve_prepare_many_start_positions(
                                wss.ssd,
                                first_region0_index + done,
                                more,
                                ws.conf.logI, ws.Q.sublat);
                        },0);
            }

            pool.drain_queue(0);

            ready = more;

            /* Now these new start positions are ready to be used */
            for(int side = 0 ; side < 2 ; side++) {
                nfs_work::side_data & wss(ws.sides[side]);
                if (wss.no_fb()) continue;
                small_sieve_activate_many_start_positions(wss.ssd);
            }
        }
    }
}/*}}}*/


/* This is the core of the sieving routine. We do fill-in-buckets,
 * downsort, apply-buckets, lognorm computation, small sieve computation,
 * and survivor search and detection, all from here.
 */
void do_one_special_q_sublat(nfs_work & ws, std::shared_ptr<nfs_work_cofac> wc_p, std::shared_ptr<nfs_aux> aux_p, thread_pool & pool)/*{{{*/
{
    nfs_aux & aux(*aux_p);
    timetree_t & timer_special_q(aux.timer_special_q);
    las_report& rep(aux.rep);
    where_am_I & w(aux.w);

    /* essentially update the fij polynomials and the max log bounds */
    if (las_output.verbose >= 2) {
        verbose_output_print (0, 1, "# f_0'(x) = ");
        mpz_poly_fprintf(las_output.output, ws.sides[0].lognorms.fij);
        verbose_output_print (0, 1, "# f_1'(x) = ");
        mpz_poly_fprintf(las_output.output, ws.sides[1].lognorms.fij);
    }

#ifdef TRACE_K
    trace_per_sq_init(ws);
#endif

    /* TODO: is there a way to share this in sublat mode ? */

    multityped_array<precomp_plattice_t, 1, FB_MAX_PARTS> precomp_plattice;

    {
        rep.ttbuckets_fill -= seconds();
        auto dummy = call_dtor([&rep]() { rep.ttbuckets_fill += seconds(); });

        {
            /* allocate_bucket_regions is probably ridiculously cheap in
             * comparison to allocate_buckets */
            ws.allocate_bucket_regions();

            ws.allocate_buckets(*aux_p, pool);
        }

        for(int side = 0 ; side < 2 ; side++) {
            nfs_work::side_data & wss(ws.sides[side]);
            if (wss.no_fb()) continue;

            fill_in_buckets_toplevel(ws, aux, pool, side, w);

            fill_in_buckets_prepare_plattices(ws, pool, side, precomp_plattice);

        }

        /*
         * Mixing-and-matching threads here with the fill-in-buckets threads
         * might lead to unbalance.
         */
        BOOKKEEPING_TIMER(timer_special_q);

        for(int side = 0 ; side < 2 ; side++) {
            nfs_work::side_data & wss(ws.sides[side]);
            if (wss.no_fb()) continue;
            pool.add_task_lambda([&ws,aux_p,side](worker_thread * worker,int){
                    timetree_t & timer(aux_p->get_timer(worker));
                    ENTER_THREAD_TIMER(timer);
                    MARK_TIMER_FOR_SIDE(timer, side);

                    SIBLING_TIMER(timer, "prepare small sieve");

                    nfs_work::side_data & wss(ws.sides[side]);
                    // if (wss.no_fb()) return;

                    small_sieve_init(wss.ssd,
                            wss.fbs->small_sieve_entries.resieved,
                            wss.fbs->small_sieve_entries.rest,
                            ws.conf.logI,
                            side,
                            wss.fbK,
                            ws.Q,
                            wss.lognorms.scale);

                    small_sieve_info("small sieve", side, wss.ssd);

                    if (ws.toplevel == 1) {
                        /* when ws.toplevel > 1, this start_many call is done
                         * several times.
                         */
                        SIBLING_TIMER(timer, "small sieve start positions ");
                        small_sieve_prepare_many_start_positions(wss.ssd,
                                0,
                                std::min(SMALL_SIEVE_START_POSITIONS_MAX_ADVANCE, ws.nb_buckets[1]),
                                ws.conf.logI, ws.Q.sublat);
                        small_sieve_activate_many_start_positions(wss.ssd);
                    }
            },0);
        }

        /* Note: we haven't done any downsorting yet ! */

        pool.drain_queue(0);

        ws.check_buckets_max_full<shorthint_t>(ws.toplevel);
        ws.check_buckets_max_full<emptyhint_t>(ws.toplevel);
        auto exc = pool.get_exceptions<buckets_are_full>(0);
        if (!exc.empty())
            throw *std::max_element(exc.begin(), exc.end());

        /* trigger the increase of ttbuckets_fill, even if we encountered an
         * exception. See the dummy object above.
         */
    }

    {
        CHILD_TIMER(timer_special_q, "process_bucket_region outer container");
        TIMER_CATEGORY(timer_special_q, sieving_mixed());
        if (ws.toplevel == 1) {
            /* Process bucket regions in parallel */
            process_many_bucket_regions(ws, wc_p, aux_p, pool, 0, w);
        } else {
            // Prepare plattices at internal levels

            // Visit the downsorting tree depth-first.
            // If toplevel = 1, then this is just processing all bucket
            // regions.
            size_t (&BRS)[FB_MAX_PARTS] = BUCKET_REGIONS;
            for (int i = 0; i < ws.nb_buckets[ws.toplevel]; i++) {
                switch (ws.toplevel) {
                    case 2:
                        downsort_tree<1>(ws, wc_p, aux_p, pool,
                                i, i*BRS[2]/BRS[1],
                                precomp_plattice, w);
                        break;
                    case 3:
                        downsort_tree<2>(ws, wc_p, aux_p, pool, i,
                                i*BRS[3]/BRS[1],
                                precomp_plattice, w);
                        break;
                    default:
                        ASSERT_ALWAYS(0);
                }
            }
        }
    }

    BOOKKEEPING_TIMER(timer_special_q);
    /* This ensures proper serialization of stuff that is in queue 0.
     * Maybe we could be looser about this.
     */
    if (sync_at_special_q) {
        pool.drain_all_queues();
    } else {
        pool.drain_queue(0);
    }
}/*}}}*/

/* This returns false if the special-q was discarded */
bool do_one_special_q(las_info & las, nfs_work & ws, std::shared_ptr<nfs_aux> aux_p, thread_pool & pool)/*{{{*/
{
    nfs_aux & aux(*aux_p);
    ws.Q.doing = aux.doing;     /* will be set by choose_sieve_area anyway */
    cxx_mpz const & p(aux.doing.p);
    cxx_mpz const & r(aux.doing.r);

    ASSERT_ALWAYS(mpz_poly_is_root(las.cpoly->pols[aux.doing.side], r, p));

    timetree_t& timer_special_q(aux.timer_special_q);
    las_report& rep(aux.rep);

    // arrange so that we don't have the same header line as the one
    // which prints the q-lattice basis
    verbose_output_print(0, 2, "#\n");
    verbose_output_vfprint(0, 1, gmp_vfprintf,
                         "# "
                         "Now sieving side-%d q=%Zd; rho=%Zd\n",
                         aux.doing.side,
                         (mpz_srcptr) aux.doing.p,
                         (mpz_srcptr) aux.doing.r);

    SIBLING_TIMER(timer_special_q, "skew Gauss");
    TIMER_CATEGORY(timer_special_q, bookkeeping());

    /* By choosing the sieve area, we trigger the setup of most of the
     * per-special-q fields.
     */
    if (!choose_sieve_area(las, aux_p, aux.doing, ws.conf, ws.Q, ws.J))
        return false;

    /* Check whether q is larger than the large prime bound.
     * This can create some problems, for instance in characters.
     * By default, this is not allowed, but the parameter
     * -allow-largesq is a by-pass to this test.
     */
    if (!allow_largesq) {
        if ((int)mpz_sizeinbase(aux.doing.p, 2) >
                ws.conf.sides[aux.doing.side].lpb) {
            fprintf(stderr, "ERROR: The special q (%d bits) is larger than the "
                    "large prime bound on side %d (%d bits).\n",
                    (int) mpz_sizeinbase(aux.doing.p, 2),
                    aux.doing.side,
                    ws.conf.sides[aux.doing.side].lpb);
            fprintf(stderr, "       You can disable this check with "
                    "the -allow-largesq argument,\n");
            fprintf(stderr, "       It is for instance useful for the "
                    "descent.\n");
	    fprintf(stderr, "       Use tasks.sieve.allow_largesq=true.\n");
            exit(EXIT_FAILURE);
        }
    }

    BOOKKEEPING_TIMER(timer_special_q);

    ws.prepare_for_new_q(las);

    /* the where_am_I structure is store in nfs_aux. We have a few
     * adjustments to make, and we want to make sure that the threads,
     * which also have their where_am_I object in nfs_aux, have it too.
     */
    WHERE_AM_I_UPDATE(aux.w, logI, ws.conf.logI);
    WHERE_AM_I_UPDATE(aux.w, pQ, &ws.Q);
    WHERE_AM_I_UPDATE(aux.w, sides[0].fbs, ws.sides[0].fbs);
    WHERE_AM_I_UPDATE(aux.w, sides[1].fbs, ws.sides[1].fbs);
    for(auto & t : aux.th) t.w = aux.w;


    /* Currently we assume that we're doing sieving + resieving on
     * both sides, or we're not. In the latter case, we expect to
     * complete the factoring work with batch cofactorization */
    ASSERT_ALWAYS(las.batch || las.batch_print_survivors || (ws.conf.sides[0].lim && ws.conf.sides[1].lim));

    std::shared_ptr<nfs_work_cofac> wc_p;

    {
        wc_p = std::make_shared<nfs_work_cofac>(las, ws);

        rep.total_logI += ws.conf.logI;
        rep.total_J += ws.J;

        std::ostringstream extra;
        if (ws.Q.doing.depth)
            extra << " # within descent, currently at depth " << ws.Q.doing.depth;

        /* should stay consistent with DUPECHECK line printed in
         * sq_finds_relation() */
        std::ostringstream os;
        os << ws.Q;
        verbose_output_vfprint(0, 2, gmp_vfprintf,
                "# "
                "Sieving %s; I=%u; J=%u;%s\n",
                os.str().c_str(),
                1u << ws.conf.logI, ws.J, extra.str().c_str());

        if (!las.allow_composite_q && !ws.Q.doing.is_prime()) {
            verbose_output_vfprint(0, 1, gmp_vfprintf,
                    "# Warning, q=%Zd is not prime\n",
                    (mpz_srcptr) ws.Q.doing.p);
        }
    }

    unsigned int sublat_bound = ws.Q.sublat.m;
    if (sublat_bound == 0)
        sublat_bound = 1;

    for (unsigned int i_cong = 0; i_cong < sublat_bound; ++i_cong) {
        for (unsigned int j_cong = 0; j_cong < sublat_bound; ++j_cong) {
            if (ws.Q.sublat.m) {
                if (i_cong == 0 && j_cong == 0)
                    continue;
                ws.Q.sublat.i0 = i_cong;
                ws.Q.sublat.j0 = j_cong;
                verbose_output_print(0, 1, "# Sublattice (i,j) == (%u, %u) mod %u\n",
                        ws.Q.sublat.i0, ws.Q.sublat.j0, ws.Q.sublat.m);
            }
            do_one_special_q_sublat(ws, wc_p, aux_p, pool);
        }
    }

    /* It's better than before, but still. We're going to keep this data
     * around for longer than we think if we get an exception above. 
     */
    for(auto & wss : ws.sides) {
        wss.precomp_plattice_dense_clear();
        small_sieve_clear(wss.ssd);
    }

    return true;
}/*}}}*/

void prepare_timer_layout_for_multithreaded_tasks(timetree_t & timer)
{
    /* This does nothing. We're just setting up the required
     * empty shell so that all the multithreaded tasks that
     * begin with MARK_TIMER_FOR_SIDE are properly registered
     * under the "sieving on side X" umbrella.
     */
#ifndef DISABLE_TIMINGS
    timetree_t::accounting_child x(timer, tdict_slot_for_threads);
    for (int side = 0; side < 2; ++side) {
        MARK_TIMER_FOR_SIDE(timer, side);
        TIMER_CATEGORY(timer, sieving(side));
    }
#else
    timer.nop();
#endif
}

void las_subjob(las_info & las, int subjob, las_todo_list & todo, las_report & global_report, timetree_t & global_timer)/*{{{*/
{
    where_am_I w MAYBE_UNUSED;
    WHERE_AM_I_UPDATE(w, plas, &las);

    las.set_subjob_binding(subjob);

    std::list<nfs_aux::caller_stuff> aux_good;
    std::list<nfs_aux::caller_stuff> aux_botched;

    double cumulated_wait_time = 0;     /* for this subjob only */
    {
        /* add scoping to control dtor call */
        /* queue 0: main
         * queue 1: ECM
         * queue 2: things that we join almost immediately, but are
         * multithreaded nevertheless: alloc buckets, ...
         */
        thread_pool pool(las.number_of_threads_per_subjob(), 3);
        nfs_work workspaces(las);

        /* {{{ Doc on todo list handling
         * The function las_todo_feed behaves in different
         * ways depending on whether we're in q-range mode or in q-list mode.
         *
         * q-range mode: the special-q's to be handled are specified as a
         * range. Then, whenever the todo list almost runs out, it is
         * refilled if possible, up to the limit q1 (-q0 & -rho just gives a
         * special case of this).
         *
         * q-list mode: the special-q's to be handled are always read from a
         * file. Therefore each new special-q to be handled incurs a
         * _blocking_ read on the file, until EOF. This mode is also used for
         * the descent, which has the implication that the read occurs if and
         * only if the todo list is empty. }}} */

        for(;;) {
            las_output.fflush();
            las_todo_entry * doing_p = todo.feed_and_pop(las.rstate);
            if (!doing_p) break;
            las_todo_entry& doing(*doing_p);

#ifdef DLP_DESCENT
            /* If the next special-q to try is a special marker, it means
             * that we're done with a special-q we started before, including
             * all its spawned sub-special-q's. Indeed, each time we start a
             * special-q from the todo list, we replace it by a special
             * marker. But newer special-q's may enver the todo list in turn
             * (pushed with las_todo_push_withdepth).
             */
            if (todo.is_closing_brace(doing)) {
                las.tree.done_node();
                if (las.tree.depth() == 0) {
                    if (recursive_descent) {
                        /* BEGIN TREE / END TREE are for the python script */
                        fprintf(las_output.output, "# BEGIN TREE\n");
                        las.tree.display_last_tree(las_output.output);
                        fprintf(las_output.output, "# END TREE\n");
                    }
                    las.tree.visited.clear();
                }
                continue;
            }

            /* pick a new entry from the stack, and do a few sanity checks */
            todo.push_closing_brace(doing.depth);
#endif

            /* We set this aside because we may have to rollback our changes
             * if we catch an exception because of full buckets. */
            auto saved_todo = todo.save();

#ifdef DLP_DESCENT
            las.tree.new_node(doing);
#endif

            /* We'll convert that to a shared_ptr later on, because this is
             * to be kept by the cofactoring tasks that will linger on quite
             * late.
             * Note that we must construct this object *outside* the try
             * block, because we want this list to be common to all attempts
             * for this q.
             */
            auto rel_hash_p = std::make_shared<nfs_aux::rel_hash_t>();

            for(;;) {
                /* We're playing a very dangerous game here. rep and timer
                 * below are created here, at the end of a temp list. When
                 * exiting this loop, they will still be live somewhere,
                 * either in aux_good or aux_botched. The nfs_aux data holds
                 * references to this stuff. But this nfs_aux is expected to
                 * live slightly longer, so it will tinker with these structs
                 * at a time where they've already been moved to aux_good (or
                 * even aux_botched, if an exception occurs late during 1l
                 * downsorting). We must be sure that all threads are done
                 * when we finally consume the aux_good and aux_botched
                 * lists!
                 *
                 * Note that the lifetime of the nfs_aux object *is*
                 * thread-safe, because the shared_ptr construct provides
                 * this guarantee. Helgrind, however, is not able to
                 * detect it properly, and sees stuff happening in the
                 * nfs_aux dtor as conflicting with what is done in the
                 * try{} block. This is a false positive. cado-nfs.supp
                 * has an explicit suppression for this.
                 */
                std::list<nfs_aux::caller_stuff> aux_pending;
                aux_pending.push_back({});
                las_report & rep(std::get<0>(aux_pending.back()));
                timetree_t & timer_special_q(std::get<1>(aux_pending.back()));

                /* ready to start over if we encounter an exception */
                try {
                    /* The nfs_aux ctor below starts the special-q timer.
                     * However we must not give it a category right now,
                     * since it is an essential property ot the timer trees
                     * that the root of the trees must not have a nontrivial
                     * category */
                    auto aux_p = std::make_shared<nfs_aux>(las, doing, rel_hash_p, aux_pending.back(), las.number_of_threads_per_subjob());
                    nfs_aux & aux(*aux_p);
                    ACTIVATE_TIMER(timer_special_q);

                    prepare_timer_layout_for_multithreaded_tasks(timer_special_q);

                    bool done = do_one_special_q(las, workspaces, aux_p, pool);

                    if (!done) {
                        /* Then we don't even keep track of the time, it's
                         * totally insignificant.
                         */
                        rep.nr_sq_discarded++;
                        break;
                    }

                    /* for statistics */
                    rep.nr_sq_processed++;

#ifdef DLP_DESCENT
                    postprocess_specialq_descent(las, todo, doing, timer_special_q);
#endif

                    if (!workspaces.cofac_candidates.empty()) {
                        {
                            std::lock_guard<std::mutex> foo(las.L.mutex());
#if 0
                            las.L.reserve(las.L.size() +
                                    workspaces.cofac_candidates.size());
                            /* do we fear that this move could be
                             * expensive ?  We're single-threaded at this
                             * point, so many threads could be waiting
                             * idly while we're moving pointers around.
                             * Worse, we have a mutex on las.L, so other
                             * subjobs might be waiting too. */
                            for(auto & x : workspaces.cofac_candidates)
                                las.L.emplace_back(std::move(x));
                            /* we can release the mutex now */
#else
                            las.L.splice(las.L.end(), workspaces.cofac_candidates);
#endif
                        }
                        workspaces.cofac_candidates.clear();
                    }

                    aux.complete = true;
                    aux_good.splice(aux_good.end(), aux_pending);

                    if (exit_after_rel_found > 1 && rep.reports > 0)
                        break;

                    break;
                } catch (buckets_are_full const & e) {
                    aux_botched.splice(aux_botched.end(), aux_pending);
                    verbose_output_vfprint (2, 1, gmp_vfprintf,
                            "# redoing q=%Zd, rho=%Zd because %s buckets are full\n"
                            "# %s\n",
                            (mpz_srcptr) doing.p, (mpz_srcptr) doing.r,
                            bkmult_specifier::printkey(e.key).c_str(),
                            e.what());

                    /* reason on the bk_multiplier that we used when we
                     * did the allocation ! It is set by
                     * prepare_for_new_q, called from do_one_special_q.
                     */
                    double old_value = workspaces.bk_multiplier.get(e.key);
                    double ratio = (double) e.reached_size / e.theoretical_max_size * 1.05;
                    double new_value = old_value * ratio;
                    double fresh_value = las.get_bk_multiplier().get(e.key);
                    if (fresh_value > new_value) {

                        verbose_output_print(0, 1, "# Global %s bucket multiplier has already grown to %.3f. Not updating, since this will cover %.3f*%d/%d*1.1=%.3f\n",
                                bkmult_specifier::printkey(e.key).c_str(),
                                fresh_value,
                                old_value,
                                e.reached_size,
                                e.theoretical_max_size,
                                new_value
                                );
                    } else {
                        verbose_output_print(0, 1, "# Updating %s bucket multiplier to %.3f*%d/%d*1.1=%.3f\n",
                                bkmult_specifier::printkey(e.key).c_str(),
                                old_value,
                                e.reached_size,
                                e.theoretical_max_size,
                                new_value
                                );
                        las.grow_bk_multiplier(e.key, ratio);
                    }
                    if (las.config_pool.default_config_ptr) {
                        expected_memory_usage(las.config_pool.base, las, true, base_memory);
                    }
                    /* we have to roll back the updates we made to
                     * this structure. */
                    todo.restore(std::move(saved_todo));
                    // las.tree.ditch_node();
                    continue;
                }
                break;
            }

        } // end of loop over special q ideals.

        /* we delete the "pool" and "workspaces" variables at this point. */
        /* The dtor for "pool" is a synchronization point */
        cumulated_wait_time = pool.cumulated_wait_time;
    }

    static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&lock);

    for(auto & P : aux_good) {
        global_report.accumulate_and_clear(std::move(std::get<0>(P)));
        global_timer += std::get<1>(P);
    }
    las_report botched_report;
    timetree_t botched_timer;
    for(auto & P : aux_botched) {
        botched_report.accumulate_and_clear(std::move(std::get<0>(P)));
        botched_timer += std::get<1>(P);
    }

    global_report.nwaste += aux_botched.size();
    global_report.waste += botched_report.tn[0];
    global_report.waste += botched_report.tn[1];
    global_report.waste += botched_report.ttbuckets_fill;
    global_report.waste += botched_report.ttbuckets_apply;
    global_report.waste += botched_report.ttf;
    global_report.waste += botched_report.ttcof;
    global_report.cumulated_wait_time += cumulated_wait_time;

    pthread_mutex_unlock(&lock);

    verbose_output_print(0, 1, "# subjob %d done, now waiting for other jobs\n", subjob);
}/*}}}*/

int main (int argc0, char *argv0[])/*{{{*/
{
    double t0, tts, wct;
    int argc = argc0;
    char **argv = argv0;

#ifdef HAVE_MINGW
    _fmode = _O_BINARY;     /* Binary open for all files */
#endif
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    cxx_param_list pl;
    las_install_sighandlers();

    declare_usage(pl);
    configure_switches(pl);
    configure_aliases(pl);

    argv++, argc--;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        /* Could also be a file */
        FILE * f;
        if ((f = fopen(argv[0], "r")) != NULL) {
            param_list_read_stream(pl, f, 0);
            fclose(f);
            argv++,argc--;
            continue;
        }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        param_list_print_usage(pl, argv0[0], stderr);
        exit(EXIT_FAILURE);
    }

    param_list_parse_int(pl, "trialdiv-first-side", &trialdiv_first_side);
    param_list_parse_int(pl, "exit-early", &exit_after_rel_found);
#if DLP_DESCENT
    param_list_parse_double(pl, "grace-time-ratio", &general_grace_time_ratio);
#endif
    param_list_parse_int(pl, "log-bucket-region", &LOG_BUCKET_REGION);
    set_LOG_BUCKET_REGION();

    las_output.set(pl);

    las_info las(pl);    /* side effects: prints cmdline and flags */

    las_todo_list todo(las.cpoly, pl);

    /* If qmin is not given, use lim on the special-q side by default.
     * This makes sense only if the relevant fields have been filled from
     * the command line.
     *
     * This is a kludge, really.
     */
    if (todo.sqside >= 0 && las.dupqmin[todo.sqside] == ULONG_MAX)
        las.dupqmin[todo.sqside] = las.config_pool.base.sides[todo.sqside].lim;

#ifdef TRACE_K
    init_trace_k(pl);
#endif

    base_memory = Memusage() << 10;

    /* First have a guess at our memory usage in single-threaded mode,
     * and see how many threads must be together. This means more memory
     * per job than our estimate, and then, maybe more threads together
     * as a consequence. This will converge eventually (of course, if
     * there is an explicit number of threads that is supplied, the
     * result is immediate).
     *
     * The significant subtlety is that the las_parallel setting matters
     * more than marginally: all subjobs within a given binding will
     * supposedly share the factor base, so we must take that into
     * account, while the las_parallel ctor will only grok a per-subjob
     * estimate.
     */
    try {
        las.set_parallel(pl);
    } catch (las_parallel_desc::needs_job_ram & e) {
        if (las.config_pool.default_config_ptr) {
            verbose_output_print(0, 2, "# No --job-memory option given, relying on automatic memory estimate\n");
            siever_config const & sc0(las.config_pool.base);
            size_t ram0 = expected_memory_usage_per_binding_zone(sc0, las, false);
            for(int z = 1, s = 1 << 30, n = 1, spin=0 ; ; spin++) {
                if (spin > 10) {
                    fprintf(stderr, "Warning: computation of expected memory does not stabilize after %d attempts, picking the situation as it is\n", spin);
                    break;
                }
                size_t ram1 = expected_memory_usage_per_subjob_worst_logI(sc0, las, n, false);
                size_t jobram = (base_memory / z + ram0) / s + ram1;
                /*
                std::ostringstream os;
                os << z << " " << s << " " << n
                    << " " << (double) ram0 / (1 << 30)
                    << " " << (double) ram1 / (1 << 30)
                    << " " << (double) jobram / (1 << 30);
                fprintf(stderr, "%s\n", os.str().c_str());
                */
                las.set_parallel(pl, (double) jobram / (1 << 30));
                int nz = las.number_of_memory_binding_zones();
                int ns = las.number_of_subjobs_per_memory_binding_zone();
                int nn = las.number_of_threads_per_subjob();
                if (s == ns && n == nn && z == nz)
                    break;
                z = nz; s = ns; n = nn;
            }
        } else {
            verbose_output_print(0, 0, "# No --job-memory option given. Job placement needs either an explicit placement with -t, a complete siever config with a factor base to allow automatic ram estimates, or a --job-memory option\n");
            exit(EXIT_FAILURE);
        }
    }
    las.display_binding_info();
    if (las.config_pool.default_config_ptr)
        expected_memory_usage(las.config_pool.base, las, true, base_memory);

    las.prepare_sieve_shared_data(pl);
    las.load_factor_base(pl);


    /* The global timer and global report structures are not active for
     * now.  Most of the accounting is done per special_q, and this will
     * be summarized once we're done with the multithreaded run.
     */
    las_report global_report;
    timetree_t global_timer;

    t0 = seconds ();
    wct = wct_seconds();

    std::vector<std::thread> subjobs;
#ifdef DLP_DESCENT
    /* In theory we would be able to to multiple descents in parallel, of
     * course, but how we should proceed with the todo list, our brace
     * mechanism, and the descent tree thing is altogether not obvious
     */
    int nsubjobs = 1;
#else
    int nsubjobs = las.number_of_subjobs_total();
#endif
    for(int subjob = 0 ; subjob < nsubjobs ; ++subjob) {
        /* when references are passed through variadic template arguments
         * as for the std::thread ctor, we have automatic decaying unless
         * we use std::ref.
         */
        subjobs.push_back(
                std::thread(las_subjob,
                    std::ref(las),
                    subjob,
                    std::ref(todo),
                    std::ref(global_report),
                    std::ref(global_timer)
                ));
    }
    for(auto & t : subjobs) t.join();

    verbose_output_print(0, 1, "# Cumulated wait time over all threads %.2f\n", global_report.cumulated_wait_time);

#ifdef DLP_DESCENT
    if (recursive_descent) {
        verbose_output_print(0, 1, "# Now displaying again the results of all descents\n");
        las.tree.display_all_trees(las_output.output);
    }
#endif

    las.set_loose_binding();

    if (las.batch)
      {
          timetree_t batch_timer;
          auto z = call_dtor([&]() { global_timer += batch_timer; });
          ACTIVATE_TIMER(batch_timer);

        /* We need to access lim[01] and lpb[01] */
        siever_config const & sc0(las.config_pool.base);
        CHILD_TIMER(batch_timer, "batch cofactorization (time is wrong because of openmp)");
        TIMER_CATEGORY(batch_timer, batch_mixed());
        double extra_time = 0;

        std::array<cxx_mpz, 2> batchP;
        int lpb[2] = { sc0.sides[0].lpb, sc0.sides[1].lpb };

        for(int side = 0 ; side < 2 ; side++) {
            create_batch_file (las.batch_file[side],
                    batchP[side],
                    sc0.sides[side].lim,
                    1UL << las.batchlpb[side],
                    las.cpoly->pols[side],
                    las_output.output,
                    las.number_of_threads_loose(),
                    extra_time);
        }

	double tcof_batch = seconds ();

        /* This one uses openmp, and forks from the current thread (well,
         * I believe so -- it's not entirely clear how openmp deals with
         * cpu binding. At least I presume that it does nothing before
         * the first pragma omp statement.)
         */
	find_smooth (las.L,
                batchP, las.batchlpb, lpb, las.batchmfb,
                las_output.output,
                las.number_of_threads_loose(),
                extra_time);

        /* We may go back to our general thread placement at this point.
         * Currently the code below still uses openmp */

        std::list<relation> rels = factor (las.L,
                las.cpoly,
                las.batchlpb,
                lpb,
		las_output.output,
                las.number_of_threads_loose(),
                extra_time);
        verbose_output_print (0, 1, "# batch reported time for additional threads: %.2f\n", extra_time);
        batch_timer.add_foreign_time(extra_time);

        verbose_output_start_batch();
        nfs_aux::rel_hash_t rel_hash;
        size_t nondup = 0;
        for(auto const & rel : rels) {
            std::ostringstream os;
            nfs_aux::abpair_t ab(rel.a, rel.b);
            bool is_new_rel = rel_hash.insert(ab).second;
            if (!is_new_rel) {
                /* we had this (a,b) pair twice, probably because of a
                 * failed attempt, that was aborted because of an
                 * exception. (occurs only with 2-level sieving) */
                os << "# DUP ";
            } else {
                nondup++;
            }
            os << rel;
            verbose_output_print(0, 1, "%s\n", os.str().c_str());
        }
        verbose_output_end_batch();
        global_report.reports = nondup;

	tcof_batch = seconds () - tcof_batch;
	global_report.ttcof += tcof_batch;
      }

    t0 = seconds () - t0;
    wct = wct_seconds() - wct;

    if (las.adjust_strategy < 2) {
        verbose_output_print (2, 1, "# Average J=%1.0f for %lu special-q's, max bucket fill -bkmult %s\n",
                global_report.total_J / (double) global_report.nr_sq_processed, global_report.nr_sq_processed, las.get_bk_multiplier().print_all().c_str());
    } else {
        verbose_output_print (2, 1, "# Average logI=%1.1f for %lu special-q's, max bucket fill -bkmult %s\n",
                global_report.total_logI / (double) global_report.nr_sq_processed, global_report.nr_sq_processed, las.get_bk_multiplier().print_all().c_str());
    }
    verbose_output_print (2, 1, "# Discarded %lu special-q's out of %u pushed\n",
            global_report.nr_sq_discarded, todo.nq_pushed);

    if (tdict::global_enable >= 1) {
        verbose_output_print (0, 1, "#\n# Hierarchical timings:\n%s", global_timer.display().c_str());

        timetree_t::timer_data_type t = 0;
        auto D = global_timer.filter_by_category();
        for(auto const &c : D)
            t += c.second;
        std::ostringstream os;
        os << std::fixed << std::setprecision(2) << t;
        verbose_output_print (0, 1, "#\n# Categorized timings (total counted time %s):\n", os.str().c_str());
        for(auto const &c : D) {
            std::ostringstream xos;
            xos << std::fixed << std::setprecision(2) << c.second;
            verbose_output_print (0, 1, "# %s: %s\n", 
                    coarse_las_timers::explain(c.first).c_str(),
                    xos.str().c_str());
        }
        verbose_output_print (0, 1, "# total counted time: %s\n#\n", os.str().c_str());
    }
    global_report.display_survivor_counters();


    if (las_output.verbose)
        facul_print_stats (las_output.output);

    /*{{{ Display tally */
#ifdef HAVE_RUSAGE_THREAD
    int dont_print_tally = 0;
#else
    int dont_print_tally = 1;
#endif
    if (bucket_prime_stats) {
        verbose_output_print(2, 1, "# nr_bucket_primes = %lu, nr_div_tests = %lu, nr_composite_tests = %lu, nr_wrap_was_composite = %lu\n",
                 nr_bucket_primes, nr_div_tests, nr_composite_tests, nr_wrap_was_composite);
    }


    verbose_output_print (2, 1, "# Wasted cpu time due to %d bkmult adjustments: %1.2f\n", global_report.nwaste, global_report.waste);

    t0 -= global_report.waste;

    tts = t0;
    tts -= global_report.tn[0];
    tts -= global_report.tn[1];
    tts -= global_report.ttf;
    tts -= global_report.ttcof;

    if (dont_print_tally && las.number_of_threads_total() > 1) 
        verbose_output_print (2, 1, "# Total cpu time %1.2fs [tally available only in mono-thread]\n", t0);
    else
        verbose_output_print (2, 1, "# Total cpu time %1.2fs [norm %1.2f+%1.1f, sieving %1.1f"
                " (%1.1f + %1.1f + %1.1f),"
                " factor %1.1f (%1.1f + %1.1f)] (not incl wasted time)\n",
                t0,
                global_report.tn[0],
                global_report.tn[1],
                tts,
                global_report.ttbuckets_fill,
                global_report.ttbuckets_apply,
                tts-global_report.ttbuckets_fill-global_report.ttbuckets_apply,
		global_report.ttf + global_report.ttcof, global_report.ttf, global_report.ttcof);

    verbose_output_print (2, 1, "# Total elapsed time %1.2fs, per special-q %gs, per relation %gs\n",
                 wct, wct / (double) global_report.nr_sq_processed, wct / (double) global_report.reports);

    /* memory usage */
    if (las_output.verbose >= 1 && las.config_pool.default_config_ptr) {
        expected_memory_usage(las.config_pool.base, las, true, base_memory);
    }
    const long peakmem = PeakMemusage();
    if (peakmem > 0)
        verbose_output_print (2, 1, "# PeakMemusage (MB) = %ld \n",
                peakmem >> 10);
    if (las.suppress_duplicates) {
        verbose_output_print(2, 1, "# Total number of eliminated duplicates: %lu\n", global_report.duplicates);
    }
    verbose_output_print (2, 1, "# Total %lu reports [%1.3gs/r, %1.1fr/sq] in %1.3g elapsed s [%.1f%% CPU]\n",
            global_report.reports, t0 / (double) global_report.reports,
            (double) global_report.reports / (double) global_report.nr_sq_processed,
            wct,
            100*t0/wct);


    print_slice_weight_estimator_stats();

    /*}}}*/

    las.cofac_stats.print();

    /* In essence, almost a dtor, but we want it to be before the pl dtor */
    las_output.release();

    return 0;
}/*}}}*/

