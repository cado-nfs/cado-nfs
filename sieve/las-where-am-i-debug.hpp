#ifndef CADO_LAS_WHERE_AM_I_DEBUG_HPP
#define CADO_LAS_WHERE_AM_I_DEBUG_HPP

// IWYU pragma: private, include "las-where-am-i.hpp"

#include <array>
#include <climits>
#include <cstdint>
#include <cstddef>
#include <vector>

#include "las-where-am-i-proxy.hpp"
#include "fb-types.hpp"
#include "fb.hpp"
#include "las-config.hpp"

struct cxx_mpz;
struct las_info;
struct qlattice_basis;

/* Define CHECK_UNDERFLOW to check for underflow when subtracting
   the rounded log(p) from sieve array locations */
//#define CHECK_UNDERFLOW

/* Define TRACK_CODE_PATH in order to have the where_am_I structures
 * propagate info on the current situation of the data being handled.
 * This more or less makes the variables global, in that every function
 * can then access the totality of the variables. But it's for debug and
 * inspection purposes only.
 *
 * Note that WANT_ASSERT_EXPENSIVE, a flag which exists in broader
 * context, augments the scope of the tracking here by performing a
 * divisibility test on each sieve update. This is obviously very
 * expensive, but provides nice checking.
 *
 * Another useful tool for debugging is the sieve-area checksums that get
 * printed with verbose output (-v) enabled.
 */
#define xxxTRACK_CODE_PATH
#define xxxWANT_ASSERT_EXPENSIVE

/* TRACE_K *requires* TRACK_CODE_PATH -- or it displays rubbish */
#if defined(TRACE_K) && !defined(TRACK_CODE_PATH)
#define TRACK_CODE_PATH
#endif

/* idem for CHECK_UNDERFLOW */
#if defined(CHECK_UNDERFLOW) && !defined(TRACK_CODE_PATH)
#define TRACK_CODE_PATH
#endif

/*  where_am_I (debug) */
struct where_am_I::impl {
    int logI = 0;
    struct side_data {
        const fb_factorbase::slicing * fbs = NULL;
    };
    std::vector<side_data> sides;
    fbprime_t p = 0;        /* current prime or prime power, when applicable */
    fbroot_t r = 0;         /* current root */
    slice_index_t i = 0;    /* Slice index, if applicable */
    slice_offset_t h = 0;   /* Prime hint, if not decoded yet */
    unsigned int j = 0;     /* row number in bucket */
    unsigned int x = 0;     /* value in bucket */
    unsigned int N = 0;     /* bucket number */
    int side = 0;
    const las_info * plas = NULL;
};

#define WHERE_AM_I_UPDATE(w, field, value) (w)->field = (value)

extern int test_divisible(where_am_I const & w);

struct trace_Nx_t { unsigned int N; unsigned int x; };
#ifdef SUPPORT_LARGE_Q
struct trace_ab_t { cxx_mpz a; cxx_mpz b; };
#else
struct trace_ab_t { int64_t a; uint64_t b; };
#endif
struct trace_ij_t { int i; unsigned int j; };

extern struct trace_Nx_t trace_Nx;
extern struct trace_ab_t trace_ab;
extern struct trace_ij_t trace_ij;

extern std::vector<cxx_mpz> traced_norms;

static inline int trace_on_spot_N(unsigned int N) {
    if (trace_Nx.x == UINT_MAX) return 0;
    return N == trace_Nx.N;
}

static inline int trace_on_spot_Nx(unsigned int N, unsigned int x) {
    if (trace_Nx.x == UINT_MAX) return 0;
    return N == trace_Nx.N && x == trace_Nx.x;
}

static inline int trace_on_range_Nx(unsigned int N, unsigned int x0, unsigned int x1) {
    if (trace_Nx.x == UINT_MAX) return 0;
    return N == trace_Nx.N && x0 <= trace_Nx.x && trace_Nx.x < x1;
}

static inline int trace_on_spot_x(uint64_t x) {
    return x == (((uint64_t)trace_Nx.N) << LOG_BUCKET_REGION)
        + (uint64_t)trace_Nx.x;
}

static inline int trace_on_spot_ab(cxx_mpz const & a, cxx_mpz const & b) {
    return a == trace_ab.a && b == trace_ab.b;
}

static inline int trace_on_spot_ab(int64_t a, uint64_t b) {
    return a == trace_ab.a && b == trace_ab.b;
}

static inline int trace_on_spot_ij(int i, unsigned int j) {
    return i == trace_ij.i && j == trace_ij.j;
}

void sieve_increase_logging(unsigned char *S, const unsigned char logp, where_am_I const & w);
void sieve_increase(unsigned char *S, const unsigned char logp, where_am_I const & w);

#endif	/* CADO_LAS_WHERE_AM_I_DEBUG_HPP */
