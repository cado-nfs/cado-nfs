#include "cado.h" // IWYU pragma: keep

/* This compilation units reacts to TRACK_CODE_PATH and uses macros
 * such as WHERE_AM_I_UPDATE.
 * This compilation unit _must_ produce different object files depending
 * on the value of TRACK_CODE_PATH.
 * The WHERE_AM_I_UPDATE macro itself is defined in las-where-am-i.hpp
 */

#include <cstdint>
#include <cinttypes>
#include <climits>
#include <cstdio>
#include <cstdlib>

#include <array>
#include <memory>
#include <string>
#include <vector>

#ifdef HAVE_CXXABI_H
#include <cxxabi.h>
#endif
#ifdef HAVE_EXECINFO
#include <execinfo.h>
#endif

#include <gmp.h>

#include "cxx_mpz.hpp"
#include "fb.hpp"
#include "fb-types.hpp"
#include "las-where-am-i.hpp"
#include "las-where-am-i-debug.hpp"
#include "las-info.hpp"
#include "las-config.hpp"
#include "las-coordinates.hpp"
#include "las-norms.hpp"
#include "las-qlattice.hpp"
#include "las-siever-config.hpp"
#include "las-threads-work-data.hpp"
#include "las-where-am-i-proxy.hpp"
#include "macros.h"
#include "verbose.h"
#include "params.h"

#ifndef TRACE_K
#error "This file *must* be compiled with TRACE_K defined"
#define TRACE_K 1
#endif

int extern_trace_on_spot_ab(cxx_mpz const & a, cxx_mpz const & b) {
    return trace_on_spot_ab(a, b);
}

int extern_trace_on_spot_ab(int64_t a, uint64_t b) {
    return trace_on_spot_ab(a, b);
}


#ifdef TRACK_CODE_PATH
where_am_I::where_am_I() : pimpl{ new impl{} } { }
where_am_I::~where_am_I() { delete pimpl; }
where_am_I::where_am_I(where_am_I const & x) : pimpl(new impl(*x.pimpl)) {
}
where_am_I & where_am_I::operator=(where_am_I const & x) {
    *pimpl = *x.pimpl;
    return *this;
}
#endif

void where_am_I::decl_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "traceab", "Relation to trace, in a,b format");
    param_list_decl_usage(pl, "traceij", "Relation to trace, in i,j format");
    param_list_decl_usage(pl, "traceNx", "Relation to trace, in N,x format");
}

/* The trivial calls for when TRACE_K is *not* defined are inlines in
 * las-debug.h */


/* recall that TRACE_K requires TRACK_CODE_PATH ; so we may safely use
 * all where_am_I types here */

trace_Nx_t trace_Nx { 0, UINT_MAX};
trace_ab_t trace_ab { 0, 0 };
trace_ij_t trace_ij { 0, UINT_MAX, };

/* Those are from the parameter list. */
static std::unique_ptr<trace_ab_t> pl_ab;
static std::unique_ptr<trace_ij_t> pl_ij;
static std::unique_ptr<trace_Nx_t> pl_Nx;

/* two norms of the traced (a,b) pair */
std::vector<cxx_mpz> traced_norms;

void where_am_I::interpret_parameters(cxx_param_list & pl)
{
    struct trace_ab_t ab;
    struct trace_ij_t ij;
    struct trace_Nx_t Nx;
    int have_trace_ab = 0, have_trace_ij = 0, have_trace_Nx = 0;

#ifdef SUPPORT_LARGE_Q
    std::pair<cxx_mpz, cxx_mpz> r;
    have_trace_ab = param_list_parse(pl, "traceab", r);
    if (have_trace_ab) {
        ab.a = r.first;
        ab.b = r.second;
    }
#else
    const char *abstr = param_list_lookup_string(pl, "traceab");
    if (abstr != NULL) {
        if (sscanf(abstr, "%" SCNd64",%" SCNu64, &ab.a, &ab.b) == 2)
            have_trace_ab = 1;
        else {
            fprintf (stderr, "Invalid value for parameter: -traceab %s\n",
                     abstr);
            exit (EXIT_FAILURE);
        }
    }
#endif

    const char *ijstr = param_list_lookup_string(pl, "traceij");
    if (ijstr != NULL) {
        if (sscanf(ijstr, "%d,%u", &ij.i, &ij.j) == 2) {
            have_trace_ij = 1;
        } else {
            fprintf (stderr, "Invalid value for parameter: -traceij %s\n",
                     ijstr);
            exit (EXIT_FAILURE);
        }
    }

    const char *Nxstr = param_list_lookup_string(pl, "traceNx");
    if (Nxstr != NULL) {
        if (sscanf(Nxstr, "%u,%u", &Nx.N, &Nx.x) == 2)
            have_trace_Nx = 1;
        else {
            fprintf (stderr, "Invalid value for parameter: -traceNx %s\n",
                     Nxstr);
            exit (EXIT_FAILURE);
        }
    }
    if (have_trace_ab) pl_ab = std::unique_ptr<trace_ab_t>(new trace_ab_t(ab));
    if (have_trace_ij) pl_ij = std::unique_ptr<trace_ij_t>(new trace_ij_t(ij));
    if (have_trace_Nx) pl_Nx = std::unique_ptr<trace_Nx_t>(new trace_Nx_t(Nx));
}

/* This fills all the trace_* structures from the main one. The main
 * structure is the one for which a non-NULL pointer is passed.
 */
void where_am_I::begin_special_q(
        nfs_work const & ws,
        special_q_data_class auto const & Q)
{
    int const logI = ws.conf.logI;
    unsigned int const J = ws.J;

    /* At most one of the three coordinates must be specified */
    ASSERT_ALWAYS((pl_Nx != NULL) + (pl_ab != NULL) + (pl_ij != NULL) <= 1);

    if (pl_ab) {
      trace_ab = *pl_ab;
      /* can possibly fall outside the q-lattice. We have to check for it */
      if (Q.convert_ab_to_ij(trace_ij.i, trace_ij.j, trace_ab.a, trace_ab.b)) {
          convert_ij_to_Nx(trace_Nx.N, trace_Nx.x, trace_ij.i, trace_ij.j, logI);
      } else {
          verbose_fmt_print(3 /* TRACE_CHANNEL */, 0, "# Relation ({},{}) to be traced "
                  "is outside of the current q-lattice\n",
                  trace_ab.a, trace_ab.b);
          trace_ij.i=0;
          trace_ij.j=UINT_MAX;
          trace_Nx.N=0;
          trace_Nx.x=UINT_MAX;
          return;
      }
    } else if (pl_ij) {
        trace_ij = *pl_ij;
        Q.convert_ij_to_ab(trace_ab.a, trace_ab.b, trace_ij.i, trace_ij.j);
        convert_ij_to_Nx(trace_Nx.N, trace_Nx.x, trace_ij.i, trace_ij.j, logI);
    } else if (pl_Nx) {
        trace_Nx = *pl_Nx;
        if (trace_Nx.x < ((size_t) 1 << LOG_BUCKET_REGION)) {
            convert_Nx_to_ij(trace_ij.i, trace_ij.j, trace_Nx.N, trace_Nx.x, logI);
            Q.convert_ij_to_ab(trace_ab.a, trace_ab.b, trace_ij.i, trace_ij.j);
        } else {
            fprintf(stderr, "Error, tracing requested for x=%u but"
                    " this siever was compiled with LOG_BUCKET_REGION=%d\n",
                    trace_Nx.x, LOG_BUCKET_REGION);
            exit(EXIT_FAILURE);
        }
    }

    if ((trace_ij.j < UINT_MAX && trace_ij.j >= J)
         || (trace_ij.i < -(1L << (logI-1)))
         || (trace_ij.i >= (1L << (logI-1))))
    {
        verbose_fmt_print(3 /* TRACE_CHANNEL */, 0, "# Relation ({},{}) to be traced is "
                "outside of the current (i,j)-rectangle (i=%d j=%u)\n",
                trace_ab.a, trace_ab.b, trace_ij.i, trace_ij.j);
        trace_ij.i=0;
        trace_ij.j=UINT_MAX;
        trace_Nx.N=0;
        trace_Nx.x=UINT_MAX;
        return;
    }
    if (trace_ij.i || trace_ij.j < UINT_MAX) {
        verbose_fmt_print(3 /* TRACE_CHANNEL */, 0, "# Tracing relation (a,b)=({},{}) "
                "(i,j)=({},{}), (N,x)=({},{})\n",
                trace_ab.a, trace_ab.b, trace_ij.i, trace_ij.j, trace_Nx.N,
                trace_Nx.x);
    }

    traced_norms.resize(ws.las.cpoly->nb_polys);
    for(int side = 0 ; side < ws.las.cpoly->nb_polys ; side++) {
        int i = trace_ij.i;
        unsigned j = trace_ij.j;
        Q.sublat.adjustIJ(i, j);
        ws.sides[side].lognorms.norm(traced_norms[side], i, j, Q);
    }
}

template void where_am_I::begin_special_q(nfs_work const &, qlattice_basis const &);
template void where_am_I::begin_special_q(nfs_work const &, siqs_special_q_data const &);

int test_divisible(where_am_I const & w)
{
    /* we only check divisibility for the given (N,x) value */
    if (!trace_on_spot_Nx(w->N, w->x))
        return 1;

    /* Note that when we are reaching here through apply_one_bucket, we
     * do not know the prime number. */
    fbprime_t const p = w->p;
    if (p==0) return 1;

    const unsigned int logI = w->logI;
    const unsigned int I = 1U << logI;

    const unsigned long X = w->x + (w->N << LOG_BUCKET_REGION);
    long const i = (long) (X & (I-1)) - (long) (I/2);
    unsigned long const j = X >> logI;
    fbprime_t q;

    q = fb_is_power (p, NULL);
    if (q == 0)
        q = p;

    int const rc = mpz_divisible_ui_p (traced_norms[w->side], (unsigned long) q);

    if (rc)
        mpz_divexact_ui (traced_norms[w->side], traced_norms[w->side], (unsigned long) q);
    else
        verbose_fmt_print(3 /* TRACE_CHANNEL */, 0, "# FAILED test_divisible(p={},N={}, x={}, side {}): i = {}, j = {}, norm = {}\n",
                w->p, w->N, w->x, w->side, (long) i, j, traced_norms[w->side]);

    return rc;
}

/* {{{ helper: sieve_increase */

#if defined(HAVE_CXXABI_H) && defined(HAVE_EXECINFO)
static std::string remove_trailing_address_suffix(std::string const& a, std::string& suffix)
{
    size_t const pos = a.find('+');
    if (pos == a.npos) {
        suffix.clear();
        return a;
    }
    suffix = a.substr(pos);
    return a.substr(0, pos);
}

static std::string get_parenthesized_arg(std::string const& a, std::string& prefix, std::string& suffix)
{
    size_t const pos = a.find('(');
    if (pos == a.npos) {
        prefix=a;
        suffix.clear();
        return std::string();
    }
    size_t const pos2 = a.find(')', pos + 1);
    if (pos2 == a.npos) {
        prefix=a;
        suffix.clear();
        return std::string();
    }
    prefix = a.substr(0, pos);
    suffix = a.substr(pos2 + 1);
    return a.substr(pos+1, pos2-pos-1);
}
#endif

/* Do this so that the _real_ caller is always 2 floors up. Must *NOT* be
 * a static function, for this very reason ! */
static void sieve_increase_logging_backend(unsigned char *S, const unsigned char logp, where_am_I const & w)
{
    if (!trace_on_spot_Nx(w->N, w->x))
        return;

    ASSERT_ALWAYS(test_divisible(w));

    std::string caller;

#ifdef HAVE_EXECINFO
    {
        void * callers_addresses[3];
        char ** callers = NULL;
        backtrace(callers_addresses, 3);
        callers = backtrace_symbols(callers_addresses, 3);
        caller = callers[2];
        free(callers);
    }

    std::string xx,yy,zz;
    yy = get_parenthesized_arg(caller, xx, zz);
    if (!yy.empty()) caller = yy;

    if (caller.empty()) {
        caller="<no symbol (static?)>";
    } else {
#ifdef HAVE_CXXABI_H
        std::string address_suffix;
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
    if (w->p) 
        verbose_fmt_print(3 /* TRACE_CHANNEL */, 0, "# Add log({},side {}) = {} to "
            "S[{}] = {}, from BA[{}] -> {} [{}]\n",
            w->p, w->side, logp, w->x, *S, w->N, (unsigned char)(*S+logp), caller.c_str());
    else
        verbose_fmt_print(3 /* TRACE_CHANNEL */, 0, "# Add log(hint={},side {}) = {} to "
            "S[{}] = {}, from BA[{}] -> {} [{}]\n",
            (unsigned long) w->h, w->side, logp, w->x, *S, w->N, (unsigned char)(*S+logp), caller.c_str());
}

/* Produce logging as sieve_increase() does, but don't actually update
   the sieve array. */
void sieve_increase_logging(unsigned char *S, const unsigned char logp, where_am_I const & w)
{
    sieve_increase_logging_backend(S, logp, w);
}

/* Increase the sieve array entry *S by logp, with underflow checking
 * and tracing if desired. w is used only for trace test and output */

void sieve_increase(unsigned char *S, const unsigned char logp, where_am_I const & w)
{
    sieve_increase_logging_backend(S, logp, w);
#ifdef CHECK_UNDERFLOW
    sieve_increase_underflow_trap(S, logp, w);
#endif
    *S += logp;
}


/* This function is useful both with and without TRACE_K, as the flag
 * controlling it is CHECK_UNDERFLOW
 */
#ifdef CHECK_UNDERFLOW
void sieve_increase_underflow_trap(unsigned char *S, const unsigned char logp, where_am_I const & w)
{
    int i;
    unsigned int j;
    int64_t a;
    uint64_t b;
    static unsigned char maxdiff = ~0;

    convert_Nx_to_ij(&i, &j, w->N, w->x, w->logI);
    w->Q->convert_ij_to_ab(&a, &b, i, j);
    if ((unsigned int) logp + *S > maxdiff)
      {
        maxdiff = logp - *S;
        verbose_fmt_print(3 /* TRACE_CHANNEL */, 0, "# Error, underflow at (N,x)=({}, {}), "
                "(i,j)=({}, {}), (a,b)=({}, {}), S[x] = {}, log(%"
                FBPRIME_FORMAT ") = {}\n",
                w->N, w->x, i, j, a, b, *S, w->p, logp);
      }
    /* arrange so that the unconditional increase which comes next
     * has the effect of taking the result to maxdiff */
    *S = maxdiff - logp;
}
#endif


/* }}} */
