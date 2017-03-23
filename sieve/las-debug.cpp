#include "cado.h"

#include <limits.h>
#include <inttypes.h>
#include <string.h>
#include <stdarg.h>

#include "las-config.h"
#include "las-types.hpp"
#include "las-debug.hpp"
#include "las-coordinates.hpp"
#include "portability.h"

using namespace std;
#if defined(__GLIBC__) && (defined(TRACE_K) || defined(CHECK_UNDERFLOW))
#include <execinfo.h>   /* For backtrace. Since glibc 2.1 */
#endif
#ifdef HAVE_CXXABI_H
/* We use that to demangle C++ names */
#include <cxxabi.h>
#endif

/* The trivial calls for when TRACE_K is *not* defined are inlines in
 * las-debug.h */


/* recall that TRACE_K requires TRACK_CODE_PATH ; so we may safely use
 * all where_am_I types here */

struct trace_Nx_t trace_Nx = { 0, UINT_MAX};
struct trace_ab_t trace_ab = { 0, 0 };
struct trace_ij_t trace_ij = { 0, UINT_MAX, };

/* two norms of the traced (a,b) pair */
mpz_t traced_norms[2];

/* This fills all the trace_* structures from the main one. The main
 * structure is the one for which a non-NULL pointer is passed.
 */
void trace_per_sq_init(sieve_info const & si, const struct trace_Nx_t *Nx,
                       const struct trace_ab_t *ab,
                       const struct trace_ij_t *ij)
{
#ifndef TRACE_K
    if (Nx != NULL || ab != NULL || ij != NULL) {
        fprintf (stderr, "Error, relation tracing requested but this siever "
                 "was compiled without TRACE_K.\n");
        exit(EXIT_FAILURE);
    }
    return;
#endif
    /* At most one of the three coordinates must be specified */
    ASSERT_ALWAYS((Nx != NULL) + (ab != NULL) + (ij != NULL) <= 1);

    if (ab != NULL) {
      trace_ab = *ab;
      /* can possibly fall outside the q-lattice. We have to check for it */
      if (ABToIJ(&trace_ij.i, &trace_ij.j, trace_ab.a, trace_ab.b, si)) {
          IJToNx(&trace_Nx.N, &trace_Nx.x, trace_ij.i, trace_ij.j, si);
      } else {
          verbose_output_print(TRACE_CHANNEL, 0, "# Relation (%" PRId64 ",%" PRIu64 ") to be traced "
                  "is outside of the current q-lattice\n",
                  trace_ab.a, trace_ab.b);
          trace_ij.i=0;
          trace_ij.j=UINT_MAX;
          trace_Nx.N=0;
          trace_Nx.x=UINT_MAX;
      }
    } else if (ij != NULL) {
        trace_ij = *ij;
        IJToAB(&trace_ab.a, &trace_ab.b, trace_ij.i, trace_ij.j, si);
        IJToNx(&trace_Nx.N, &trace_Nx.x, trace_ij.i, trace_ij.j, si);
    } else if (Nx != NULL) {
        trace_Nx = *Nx;
        NxToIJ(&trace_ij.i, &trace_ij.j, trace_Nx.N, trace_Nx.x, si);
        IJToAB(&trace_ab.a, &trace_ab.b, trace_ij.i, trace_ij.j, si);
    }

    if ((trace_ij.j < UINT_MAX && trace_ij.j >= si.J)
         || (trace_ij.i < -(1L << (si.conf.logI_adjusted-1)))
         || (trace_ij.i >= (1L << (si.conf.logI_adjusted-1))))
    {
        verbose_output_print(TRACE_CHANNEL, 0, "# Relation (%" PRId64 ",%" PRIu64 ") to be traced is "
                "outside of the current (i,j)-rectangle (i=%d j=%u)\n",
                trace_ab.a, trace_ab.b, trace_ij.i, trace_ij.j);
        trace_ij.i=0;
        trace_ij.j=UINT_MAX;
        trace_Nx.N=0;
        trace_Nx.x=UINT_MAX;
    }
    if (trace_ij.i || trace_ij.j < UINT_MAX) {
        verbose_output_print(TRACE_CHANNEL, 0, "# Tracing relation (a,b)=(%" PRId64 ",%" PRIu64 ") "
                "(i,j)=(%d,%u), (N,x)=(%u,%u)\n",
                trace_ab.a, trace_ab.b, trace_ij.i, trace_ij.j, trace_Nx.N,
                trace_Nx.x);
    }

    for(int side = 0 ; side < 2 ; side++) {
        mpz_init(traced_norms[side]);
        mpz_poly_homogeneous_eval_sisi(traced_norms[side], 
                si.sides[side].fij, trace_ij.i, trace_ij.j);
    }
}

void trace_per_sq_clear(sieve_info const & si MAYBE_UNUSED)
{
    for(int side = 0 ; side < 2 ; side++)
        mpz_clear(traced_norms[side]);
}

#ifdef TRACE_K
int test_divisible(where_am_I& w)
{
    /* we only check divisibility for the given (N,x) value */
    if (!trace_on_spot_Nx(w.N, w.x))
        return 1;

    /* Note that when we are reaching here through apply_one_bucket, we
     * do not know the prime number. */
    fbprime_t p = w.p;
    if (p==0) return 1;

    const unsigned int logI = w.psi->conf.logI_adjusted;
    const unsigned int I = 1U << logI;

    const unsigned long X = w.x + (w.N << LOG_BUCKET_REGION);
    long i = (long) (X & (I-1)) - (long) (I/2);
    unsigned long j = X >> logI;
    fbprime_t q;

    q = fb_is_power (p, NULL);
    if (q == 0)
        q = p;

    int rc = mpz_divisible_ui_p (traced_norms[w.side], (unsigned long) q);

    if (rc)
        mpz_divexact_ui (traced_norms[w.side], traced_norms[w.side], (unsigned long) q);
    else
        verbose_output_vfprint(TRACE_CHANNEL, 0, gmp_vfprintf, "# FAILED test_divisible(p=%" FBPRIME_FORMAT
                ", N=%d, x=%lu, side %d): i = %ld, j = %lu, norm = %Zd\n",
                w.p, w.N, w.x, w.side, i, j, traced_norms[w.side]);

    return rc;
}

/* {{{ helper: sieve_increase */

string remove_trailing_address_suffix(string const& a, string& suffix)
{
    size_t pos = a.find('+');
    if (pos == a.npos) {
        suffix.clear();
        return a;
    }
    suffix = a.substr(pos);
    return a.substr(0, pos);
}

string get_parenthesized_arg(string const& a, string& prefix, string& suffix)
{
    size_t pos = a.find('(');
    if (pos == a.npos) {
        prefix=a;
        suffix.clear();
        return string();
    }
    size_t pos2 = a.find(')', pos + 1);
    if (pos2 == a.npos) {
        prefix=a;
        suffix.clear();
        return string();
    }
    prefix = a.substr(0, pos);
    suffix = a.substr(pos2 + 1);
    return a.substr(pos+1, pos2-pos-1);
}

/* Do this so that the _real_ caller is always 2 floors up */
void sieve_increase_logging_backend(unsigned char *S, const unsigned char logp, where_am_I& w)
{
    if (!trace_on_spot_Nx(w.N, w.x))
        return;

    ASSERT_ALWAYS(test_divisible(w));

    string caller;

#ifdef __GLIBC__
    {
        void * callers_addresses[3];
        char ** callers = NULL;
        backtrace(callers_addresses, 3);
        callers = backtrace_symbols(callers_addresses, 3);
        caller = callers[2];
        free(callers);
    }

    string xx,yy,zz;
    yy = get_parenthesized_arg(caller, xx, zz);
    if (!yy.empty()) caller = yy;

    if (caller.empty()) {
        caller="<no symbol (static?)>";
    } else {
#ifdef HAVE_CXXABI_H
        string address_suffix;
        caller = remove_trailing_address_suffix(caller, address_suffix);
        int demangle_status;
        {
            char * freeme = abi::__cxa_demangle(caller.c_str(), 0, 0, &demangle_status);
            if (demangle_status == 0) {
                caller = freeme;
                free(freeme);
            }
        }

        /* Get rid of the type signature, it rather useless */
        yy = get_parenthesized_arg(caller, xx, zz);
        caller = xx;

        /* could it be that we have the return type in the name
         * as well ? */

        caller+=address_suffix;
#endif
    }
#endif
    if (w.p) 
        verbose_output_print(TRACE_CHANNEL, 0, "# Add log(%" FBPRIME_FORMAT ",side %d) = %hhu to "
            "S[%u] = %hhu, from BA[%u] -> %hhu [%s]\n",
            w.p, w.side, logp, w.x, *S, w.N, (unsigned char)(*S+logp), caller.c_str());
    else
        verbose_output_print(TRACE_CHANNEL, 0, "# Add log(hint=%lu,side %d) = %hhu to "
            "S[%u] = %hhu, from BA[%u] -> %hhu [%s]\n",
            (unsigned long) w.h, w.side, logp, w.x, *S, w.N, (unsigned char)(*S+logp), caller.c_str());
}

void sieve_increase_logging(unsigned char *S, const unsigned char logp, where_am_I& w)
{
    sieve_increase_logging_backend(S, logp, w);
}

/* Increase the sieve array entry *S by logp, with underflow checking
 * and tracing if desired. w is used only for trace test and output */

void sieve_increase(unsigned char *S, const unsigned char logp, where_am_I& w)
{
    sieve_increase_logging_backend(S, logp, w);
#ifdef CHECK_UNDERFLOW
    sieve_increase_underflow_trap(S, logp, w);
#endif
    *S += logp;
}

#endif  /* TRACE_K */

/* This function is useful both with and without TRACE_K, as the flag
 * controlling it is CHECK_UNDERFLOW
 */
#ifdef CHECK_UNDERFLOW
void sieve_increase_underflow_trap(unsigned char *S, const unsigned char logp, where_am_I& w)
{
    int i;
    unsigned int j;
    int64_t a;
    uint64_t b;
    static unsigned char maxdiff = ~0;

    NxToIJ(&i, &j, w.N, w.x, w.si);
    IJToAB(&a, &b, i, j, w.si);
    if ((unsigned int) logp + *S > maxdiff)
      {
        maxdiff = logp - *S;
        verbose_output_print(TRACE_CHANNEL, 0, "# Error, underflow at (N,x)=(%u, %u), "
                "(i,j)=(%d, %u), (a,b)=(%ld, %lu), S[x] = %hhu, log(%"
                FBPRIME_FORMAT ") = %hhu\n",
                w.N, w.x, i, j, a, b, *S, w.p, logp);
      }
    /* arrange so that the unconditional increase which comes next
     * has the effect of taking the result to maxdiff */
    *S = maxdiff - logp;
}
#endif

/* }}} */

