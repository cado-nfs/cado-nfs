#include "cado.h" // IWYU pragma: keep

#define xxxLOG_BUCKET_REGION_IS_A_CONSTANT

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>

// IWYU pragma: no_include <ext/alloc_traits.h>
#include <vector>
#include <sstream>
#include <algorithm>
#include <memory>                             // for allocator_traits<>::val...
#include <string>                             // for string
#include <type_traits>                        // for is_same

#include <sys/time.h>

#include <gmp.h>

#include "fb-types.hpp"
#include "las-smallsieve.hpp"           // for ssp_simple_t, ssp_t
#include "las-where-am-i-proxy.hpp"           // for where_am_I
#include "las-where-am-i.hpp"
#include "macros.h"
#include "params.h"
#include "sieve/las-smallsieve-glue.hpp"
#include "sieve/las-smallsieve-lowlevel.hpp"
#include "test-smallsieve-mock.hpp"


static int consistency_check_mode = 0;
static int quiet = 0;
static int abort_on_fail = 0;
static int only_complete_functions = 0;

/* this is really a mock structure just for the fun of it. */
static sublat_t sl;




/******************************************************************/
/* now provides routines for testing */

template<typename T, typename U> struct list_car {};

/* we create lists of functions, because from all these tidbits we need
 * to generate the more complete sieving functions below.
 */
struct all_generic_candidates {
    typedef
            list_car<manual_oldloop,
            list_car<manual_oldloop_nounroll,
            list_nil>> type0;
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    typedef list_car<assembly_generic_loop8,
            list_car<assembly_generic_loop12,
            list_car<assembly_generic_loop16,
            list_car<assembly_generic_loop8p0,
            list_car<assembly_generic_loop12p0,
            list_car<assembly_generic_loop16p0,
            list_car<assembly_generic_loop8p1,
            list_car<assembly_generic_loop12p1,
            list_car<assembly_generic_loop16p1,
            list_car<assembly_generic_loop12p2,
            list_car<assembly_generic_oldloop,
            type0> >>> >>> >>> > type;
#else
    using type = type0;
#endif
};

template<int bit> struct all_candidates_for_evenline {
    using type = all_generic_candidates::type;
};
template<int bit> struct all_candidates_for_oddline {
    using type = all_generic_candidates::type;
};
template<> struct all_candidates_for_evenline<0> {
    using type0 = list_car<manual0, all_generic_candidates::type>;
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    using type = list_car<assembly0, type0>;
#else
    using type = type0;
#endif
};
template<> struct all_candidates_for_evenline<1> {
    using type0 = list_car<manual0, all_generic_candidates::type>;
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    using type = list_car<assembly0, type0>;
#else
    using type = type0;
#endif
};
template<> struct all_candidates_for_evenline<2> {
    using type0 = list_car<manual1, all_generic_candidates::type>;
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    using type = list_car<assembly1, type0>;
#else
    using type = type0;
#endif
};
template<> struct all_candidates_for_evenline<3> {
    using type0 = list_car<manual2, all_generic_candidates::type>;
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    using type = list_car<assembly2, type0>;
#else
    using type = type0;
#endif
};
template<> struct all_candidates_for_evenline<4> {
    using type0 = list_car<manual3, all_generic_candidates::type>;
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    using type = list_car<assembly3, type0>;
#else
    using type = type0;
#endif
};
template<> struct all_candidates_for_evenline<5> {
    using type0 = list_car<manual4, all_generic_candidates::type>;
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    using type = list_car<assembly4, type0>;
#else
    using type = type0;
#endif
};
template<> struct all_candidates_for_evenline<6> {
    using type0 = all_generic_candidates::type;
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    using type = list_car<assembly5, type0>;
#else
    using type = type0;
#endif
};
template<> struct all_candidates_for_evenline<7> {
    using type0 = all_generic_candidates::type;
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    using type = list_car<assembly6, type0>;
#else
    using type = type0;
#endif
};
template<> struct all_candidates_for_oddline<0> {
    using type0 = list_car<manual0, all_generic_candidates::type>;
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    using type = list_car<assembly0, type0>;
#else
    using type = type0;
#endif
};
template<> struct all_candidates_for_oddline<1> {
    using type0 = list_car<manual1, all_generic_candidates::type>;
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    using type = list_car<assembly1, type0>;
#else
    using type = type0;
#endif
};
template<> struct all_candidates_for_oddline<2> {
    using type0 = list_car<manual2, all_generic_candidates::type>;
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    using type = list_car<assembly2, type0>;
#else
    using type = type0;
#endif
};
template<> struct all_candidates_for_oddline<3> {
    using type0 = list_car<manual3, all_generic_candidates::type>;
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    using type = list_car<assembly3, type0>;
#else
    using type = type0;
#endif
};
template<> struct all_candidates_for_oddline<4> {
    using type0 = list_car<manual4, all_generic_candidates::type>;
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    using type = list_car<assembly4, type0>;
#else
    using type = type0;
#endif
};
template<> struct all_candidates_for_oddline<5> {
    using type0 = all_generic_candidates::type;
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    using type = list_car<assembly5, type0>;
#else
    using type = type0;
#endif
};
template<> struct all_candidates_for_oddline<6> {
    using type0 = all_generic_candidates::type;
#ifdef HAVE_GCC_STYLE_AMD64_INLINE_ASM
    using type = list_car<assembly6, type0>;
#else
    using type = type0;
#endif
};

/* For 2^(I-k) <= p < 2^(I-k+1), we have: 2^k >= 2^I/p > 2^(k-1), and
 *
 * so that the number of intervals of width p that fit within 2^I,
 * which is (2^I-1)/p, is between:
 *
 * 2^k > (2^I-1)/p >= 2^(k-1)
 *
 * The number of hits (endpoints of the intervals) is thus <= 2^k.
 *
 * For a lower bound, 2^(k-1) is attained because of the offset.
 *
 * Note that on even lines, p behaves as if it was doubled.
 */

static inline size_t sieve_full_line(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_OLD(pi, p_or_2p, S1, logp);
    return pi - S1;
}
static inline size_t sieve_full_line_new_half(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_NEW_LOOP12(pi, p_or_2p, S1, logp);
    return pi - S1;
}
static inline size_t sieve_full_line_new(unsigned char * S0, unsigned char * S1, size_t x0 MAYBE_UNUSED, size_t pos, size_t p_or_2p, unsigned char logp, where_am_I w MAYBE_UNUSED)
{
    unsigned char * pi = S0 + pos;
    SMALLSIEVE_ASSEMBLY_NEW_LOOP12(pi, p_or_2p, S1, logp);
    return pi - S1;
}

static void current_I18_branch(std::vector<int> const & positions, std::vector<ssp_simple_t> const& primes, unsigned char * S, int logI, unsigned int N, where_am_I & w MAYBE_UNUSED) /* {{{ */
{
    SMALLSIEVE_COMMON_DEFS();
    ASSERT_ALWAYS(positions.size() == primes.size());
    for(auto const & ssp : primes) {
        int pos = positions[&ssp - primes.data()];

        const fbprime_t p = ssp.get_p();
        const fbprime_t r = ssp.get_r();
        const unsigned char logp = ssp.logp;
        unsigned char * S0 = S;

        unsigned int const i_compens_sublat = sl.i0 & 1;
        unsigned int j = j0;


        /* we sieve over the area [S0..S0+(i1-i0)], which may
         * actually be just a fragment of a line. After that, if
         * (i1-i0) is different from I, we'll break anyway. So
         * whether we add I or (i1-i0) to S0 does not matter much.
         */
        if (j & 1) {
            WHERE_AM_I_UPDATE(w, j, j - j0);
            sieve_full_line(S0, S0 + (i1 - i0), S0 - S,
                    pos, p, logp, w);
            S0 += I;
            j++;
            pos += r; if (pos >= (int) p) pos -= p;
        }
        for( ; j < j1 ; ) {
            /* for j even, we sieve only odd pi, so step = 2p. */
            int const xpos = ((i_compens_sublat + pos) & 1) ? pos : (pos+p);
            sieve_full_line(S0, S0 + (i1 - i0), S0 - S,
                    xpos, p+p, logp, w);
            S0 += I;
            pos += r; if (pos >= (int) p) pos -= p;
            if (++j >= j1) break;

            /* now j odd again */
            WHERE_AM_I_UPDATE(w, j, j - j0);
            sieve_full_line(S0, S0 + (i1 - i0), S0 - S,
                    pos, p, logp, w);
            S0 += I;
            pos += r; if (pos >= (int) p) pos -= p;
            ++j;
        }
    }
}/*}}}*/
#if 0
void modified_I18_branch_C(std::vector<int> const & positions, std::vector<ssp_simple_t> const& primes, unsigned char * S, int logI, unsigned int N, where_am_I & w MAYBE_UNUSED) /* {{{ */
{
    SMALLSIEVE_COMMON_DEFS();
    ASSERT_ALWAYS(positions.size() == primes.size());
    for(auto const & ssp : primes) {
        int pos = positions[&ssp - primes.data()];

        const fbprime_t p = ssp.get_p();
        const fbprime_t r = ssp.get_r();
        const unsigned char logp = ssp.logp;
        unsigned char * S0 = S;

        /* we sieve over the area [S0..S0+(i1-i0)], which may
         * actually be just a fragment of a line. After that, if
         * (i1-i0) is different from I, we'll break anyway. So
         * whether we add I or (i1-i0) to S0 does not matter much.
         */
        fbprime_t const spx = (p+p)^p;
        fbprime_t px = (j0&1)?p:(p+p);
        /* The increment is i1-i0, not I, to cover the case where
         * B <= I; we then have i1-i0 = B
         * because there's a min. (and in that case, there's only a
         * single line).
         */
        for(unsigned int j = j0 ; j < j1 ; ++j) {
            /* for j even, we sieve only odd pi, so step = 2p. */
            WHERE_AM_I_UPDATE(w, j, j);
            sieve_full_line(S0, S0 + (i1 - i0), S0 - S,
                    pos + (((sl.i0^pos^1)&1&~j)?p:0), px, logp, w);
            S0 += I;
            pos += r; if (pos >= (int) p) pos -= p;
            px ^= spx;
        }
    }
}/*}}}*/
#endif
static void legacy_branch(std::vector<int> const & positions, std::vector<ssp_simple_t> const& primes, unsigned char * S, int logI, unsigned int N, where_am_I & w MAYBE_UNUSED) /*{{{*/
{
    SMALLSIEVE_COMMON_DEFS();
    const unsigned long bucket_region = (1UL << LOG_BUCKET_REGION);
    ASSERT_ALWAYS(positions.size() == primes.size());
    ASSERT_ALWAYS(logI <= LOG_BUCKET_REGION);

    for(auto const & ssp : primes) {
        int pos = positions[&ssp - primes.data()];

        const fbprime_t p = ssp.get_p();
        const fbprime_t r = ssp.get_r();
        const unsigned char logp = ssp.logp;

        unsigned long j;
        const int test_divisibility MAYBE_UNUSED = 0; /* very slow, but nice for debugging */
        const unsigned long nj = bucket_region >> logI; /* Nr. of lines per bucket region */

        WHERE_AM_I_UPDATE(w, p, p);

        //XXX XXX XXX const unsigned char logp = ssd->logp[k];
        unsigned char *S_ptr = S;
        size_t p_or_2p = p;
        ASSERT((fbprime_t) pos < p);
        unsigned int i_compens_sublat = 0;
        if (sl.m != 0) {
            i_compens_sublat = sl.i0 & 1;
        }
        j = 0;
        if (row0_is_oddj) goto j_odd;
        // it is possible to make this code work with sieve_full_line by
        // uncommenting all lines that match S0. Performance is
        // apparently identical, but it's not fair to say that this *is*
        // the reference code, as far as performance is concerned.
        // unsigned char * S0 MAYBE_UNUSED;
        do {
            unsigned char *pi;
            WHERE_AM_I_UPDATE(w, j, j);
            // S0 = S_ptr;
            pi = S_ptr + pos;
            S_ptr += I;
            /* for j even, we sieve only odd pi, so step = 2p. */
            p_or_2p += p_or_2p;
            if (!((i_compens_sublat + pos) & 1)) {
                pi += p;
            }
            SMALLSIEVE_ASSEMBLY_OLD(pi, p_or_2p, S_ptr, logp);
            // sieve_full_line(S0, S_ptr, S0 - S, pi - S0, p_or_2p, logp, w);
            pos += r; if (pos >= (int) p) pos -= p;
            /* Next line */
            if (++j >= nj) break;
            p_or_2p >>= 1;
j_odd:
            WHERE_AM_I_UPDATE(w, j, j);
            // S0 = S_ptr;
            pi = S_ptr + pos;
            S_ptr += I;
            SMALLSIEVE_ASSEMBLY_OLD(pi, p_or_2p, S_ptr, logp);
            // sieve_full_line(S0, S_ptr, S0 - S, pi - S0, p_or_2p, logp, w);
            pos += r; if (pos >= (int) p) pos -= p;
        } while (++j < nj);
    }
}
/*}}}*/
static void devel_branch(std::vector<int> const & positions, std::vector<ssp_simple_t> const& primes, unsigned char * S, int logI, unsigned int N, where_am_I & w MAYBE_UNUSED) /*{{{*/
{
    SMALLSIEVE_COMMON_DEFS();
    ASSERT_ALWAYS(positions.size() == primes.size());
    for(auto const & ssp : primes) {
        int pos = positions[&ssp - primes.data()];

        const fbprime_t p = ssp.get_p();
        const fbprime_t r = ssp.get_r();
        const unsigned char logp = ssp.logp;
        unsigned char * S0 = S;

        unsigned int j = j0;

        /* we sieve over the area [S0..S0+(i1-i0)], which may
         * actually be just a fragment of a line. After that, if
         * (i1-i0) is different from I, we'll break anyway. So
         * whether we add I or (i1-i0) to S0 does not matter much.
         */
        /* TODO: sublat !!!! j actually corresponds to row S(j), with
         * S(j)=sublatm*j + sublatj0. Therefore:
         *  - j0&1 == 0 is not sufficient to tell the parity of S(j)
         *  - if sublatm is even, we not every other line is even (resp.
         *    odd).
         *  - if sublatm is even, we should actually take advantage of
         *    it, and run specific code, probably at the outer loop
         *    level.
         */
        if (j0 & 1)
            goto j_odd_devel;

        for( ; j < j1 ; ) {
            /* for S(j) even, we sieve only odd pi, so step = 2p. */
            {
            int const xpos = ((sl.i0 + pos) & 1) ? pos : (pos+p);
            sieve_full_line_new_half(S0, S0 + (i1 - i0), S0 - S,
                    xpos, p+p, logp, w);
            }
            S0 += I;
            pos += r; if (pos >= (int) p) pos -= p;
            if (++j >= j1) break;

j_odd_devel:
            /* now S(j) odd again */
            WHERE_AM_I_UPDATE(w, j, j - j0);
            sieve_full_line_new(S0, S0 + (i1 - i0), S0 - S,
                    pos, p, logp, w);
            S0 += I;
            pos += r; if (pos >= (int) p) pos -= p;
            ++j;
        }
    }
}
/*}}}*/
static void legacy_mod_branch(std::vector<int> const & positions, std::vector<ssp_simple_t> const& primes, unsigned char * S, int logI, unsigned int N, where_am_I & w MAYBE_UNUSED) /*{{{*/
{
    SMALLSIEVE_COMMON_DEFS();
    const unsigned long bucket_region = (1UL << LOG_BUCKET_REGION);
    ASSERT_ALWAYS(positions.size() == primes.size());

    if (logI > LOG_BUCKET_REGION) {
        for(auto const & ssp : primes) {
            int pos = positions[&ssp - primes.data()];

            const fbprime_t p = ssp.get_p();
            const fbprime_t r = ssp.get_r();
            const unsigned char logp = ssp.logp;

            unsigned long j;
            const int test_divisibility MAYBE_UNUSED = 0; /* very slow, but nice for debugging */
            const unsigned long nj = bucket_region >> logI; /* Nr. of lines per bucket region */

            WHERE_AM_I_UPDATE(w, p, p);

            //XXX XXX XXX const unsigned char logp = ssd->logp[k];
            unsigned char *S_ptr = S;
            unsigned char *S1;
            size_t p_or_2p = p;
            ASSERT(pos < (int) p);
            unsigned int i_compens_sublat = 0;
            if (sl.m != 0) {
                i_compens_sublat = sl.i0 & 1;
            }
            j = 0;
            if (row0_is_oddj) goto j_odd0;
            // it is possible to make this code work with sieve_full_line by
            // uncommenting all lines that match S0. Performance is
            // apparently identical, but it's not fair to say that this *is*
            // the reference code, as far as performance is concerned.
            // unsigned char * S0 MAYBE_UNUSED;
            do {
                unsigned char *pi;
                WHERE_AM_I_UPDATE(w, j, j);
                // S0 = S_ptr;
                pi = S_ptr + pos;
                S1 = S_ptr + (i1 - i0);
                S_ptr += I;
                /* for j even, we sieve only odd pi, so step = 2p. */
                p_or_2p += p_or_2p;
                if (!((i_compens_sublat + pos) & 1)) {
                    pi += p;
                }
                SMALLSIEVE_ASSEMBLY_OLD(pi, p_or_2p, S1, logp);
                // sieve_full_line(S0, S_ptr, S0 - S, pi - S0, p_or_2p, logp, w);
                pos += r; if (pos >= (int) p) pos -= p;
                /* Next line */
                if (++j >= nj) break;
                p_or_2p >>= 1;
j_odd0:
                WHERE_AM_I_UPDATE(w, j, j);
                // S0 = S_ptr;
                pi = S_ptr + pos;
                S1 = S_ptr + (i1 - i0);
                S_ptr += I;
                SMALLSIEVE_ASSEMBLY_OLD(pi, p_or_2p, S1, logp);
                pos += r; if (pos >= (int) p) pos -= p;
            } while (++j < nj);
        }
    } else {
        for(auto const & ssp : primes) {
            int pos = positions[&ssp - primes.data()];

            const fbprime_t p = ssp.get_p();
            const fbprime_t r = ssp.get_r();
            const unsigned char logp = ssp.logp;

            unsigned long j;
            const int test_divisibility MAYBE_UNUSED = 0; /* very slow, but nice for debugging */
            const unsigned long nj = bucket_region >> logI; /* Nr. of lines per bucket region */

            WHERE_AM_I_UPDATE(w, p, p);

            //XXX XXX XXX const unsigned char logp = ssd->logp[k];
            unsigned char *S_ptr = S;
            size_t p_or_2p = p;
            ASSERT(pos < (int) p);
            unsigned int i_compens_sublat = 0;
            if (sl.m != 0) {
                i_compens_sublat = sl.i0 & 1;
            }
            j = 0;
            if (row0_is_oddj) goto j_odd;
            // it is possible to make this code work with sieve_full_line by
            // uncommenting all lines that match S0. Performance is
            // apparently identical, but it's not fair to say that this *is*
            // the reference code, as far as performance is concerned.
            // unsigned char * S0 MAYBE_UNUSED;
            do {
                unsigned char *pi;
                WHERE_AM_I_UPDATE(w, j, j);
                // S0 = S_ptr;
                pi = S_ptr + pos;
                S_ptr += I;
                /* for j even, we sieve only odd pi, so step = 2p. */
                p_or_2p += p_or_2p;
                if (!((i_compens_sublat + pos) & 1)) {
                    pi += p;
                }
                SMALLSIEVE_ASSEMBLY_OLD(pi, p_or_2p, S_ptr, logp);
                pos += r; if (pos >= (int) p) pos -= p;
                /* Next line */
                if (++j >= nj) break;
                p_or_2p >>= 1;
j_odd:
                WHERE_AM_I_UPDATE(w, j, j);
                // S0 = S_ptr;
                pi = S_ptr + pos;
                S_ptr += I;
                SMALLSIEVE_ASSEMBLY_OLD(pi, p_or_2p, S_ptr, logp);
                pos += r; if (pos >= (int) p) pos -= p;
            } while (++j < nj);
        }
    }
}
/*}}}*/

template<typename even_code, typename odd_code, bool fragment>
static void devel_branch_meta(std::vector<int> const & positions, std::vector<ssp_simple_t> const& primes, unsigned char * S, int logI, unsigned int N, where_am_I & w MAYBE_UNUSED) /*{{{*/
{
    SMALLSIEVE_COMMON_DEFS();
    ASSERT_ALWAYS(positions.size() == primes.size());
    for(auto const & ssp : primes) {
        int pos = positions[&ssp - primes.data()];

        const fbprime_t p = ssp.get_p();
        const fbprime_t r = ssp.get_r();
        const unsigned char logp = ssp.logp;
        unsigned char * S0 = S;

        unsigned int j = j0;

        /* we sieve over the area [S0..S0+(i1-i0)], which may
         * actually be just a fragment of a line. After that, if
         * (i1-i0) is different from I, we'll break anyway. So
         * whether we add I or (i1-i0) to S0 does not matter much.
         */
        /* TODO: sublat !!!! j actually corresponds to row S(j), with
         * S(j)=sublatm*j + sublatj0. Therefore:
         *  - j0&1 == 0 is not sufficient to tell the parity of S(j)
         *  - if sublatm is even, we not every other line is even (resp.
         *    odd).
         *  - if sublatm is even, we should actually take advantage of
         *    it, and run specific code, probably at the outer loop
         *    level.
         */
        if (j0 & 1)
                goto j_odd_devel0;

        for( ; j < j1 ; ) {
            /* for j even, we sieve only odd pi, so step = 2p. */
            {
            int const xpos = ((sl.i0 + pos) & 1) ? pos : (pos+p);
            even_code()(S0, S0 + (i1 - i0), S0 - S, xpos, p+p, logp, w);
            }
            S0 += I;
            pos += r; if (pos >= (int) p) pos -= p;
            if (++j >= j1) break;

j_odd_devel0:
            /* now j odd again */
            WHERE_AM_I_UPDATE(w, j, j - j0);
            odd_code()(S0, S0 + (i1 - i0), S0 - S, pos, p, logp, w);
            S0 += I;
            pos += r; if (pos >= (int) p) pos -= p;
            ++j;
        }
    }
}/*}}}*/


static void generated(std::vector<int> const & positions, std::vector<ssp_simple_t> const & primes, unsigned char * S, int logI, unsigned int N, where_am_I & w MAYBE_UNUSED) /*{{{*/
{
    std::vector<ssp_t> const not_nice_primes;
    small_sieve SS(positions, primes, not_nice_primes, S, logI, N, sl);
    // SS.pattern_sieve2(w);
    // SS.pattern_sieve3(w);
    SS.normal_sieve(w);
}
/*}}}*/


struct cand_func {
    typedef void (*ss_func)(std::vector<int> const & positions, std::vector<ssp_simple_t> const& primes, unsigned char * S, int logI, unsigned int N, where_am_I & w MAYBE_UNUSED);
    bool sel;
    ss_func f;
    const char * name;
    cand_func(bool sel, ss_func f, const char * name) : sel(sel), f(f), name(name) {}
};

using candidate_list = std::vector<cand_func>;

/* code factory: first we create a list of functions that have the
 * currently recorded best code for even lines, and then we iterate on
 * the various options for odd lines
 */
template<int bit, bool fragment> struct factory_for_bit_round1 {
    using Be = typename best_evenline<bit>::type;
    using Bo = typename best_oddline<bit>::type;
    template<typename T> struct iterator {
        void operator()(candidate_list &) {}
    };
    template<typename T, typename U> struct iterator<list_car<T, U>> {
        void operator()(candidate_list & res) const {
            bool const sel = std::is_same<Bo, T>::value;
            res.push_back(cand_func {sel, &devel_branch_meta<Be, T, fragment>, T::name});
            iterator<U>()(res);
        }
    };
    void operator()(candidate_list & res) const {
        iterator<typename all_candidates_for_oddline<bit>::type>()(res);
    }
};
/* and now the other way around. */
template<int bit, bool fragment> struct factory_for_bit_round2 {
    using Be = typename best_evenline<bit>::type;
    using Bo = typename best_oddline<bit>::type;
    template<typename T> struct iterator {
        void operator()(candidate_list &) {}
    };
    template<typename T, typename U> struct iterator<list_car<T, U>> {
        void operator()(candidate_list & res) const {
            bool const sel = std::is_same<Be, T>::value;
            res.push_back(cand_func {sel, &devel_branch_meta<T, Bo, fragment>, T::name});
            iterator<U>()(res);
        }
    };
    void operator()(candidate_list & res) const {
        iterator<typename all_candidates_for_evenline<bit>::type>()(res);
    }
};

struct bench_base {
#define BBLACK(X) "\e[01;30m" X "\e[00m"
#define BRED(X)   "\e[01;31m" X "\e[00m"
#define BGREEN(X) "\e[01;32m" X "\e[00m"
#define BYELLOW(X)"\e[01;33m" X "\e[00m"
#define BBLUE(X)  "\e[01;34m" X "\e[00m"
#define BVIOLET(X)"\e[01;35m" X "\e[00m"
#define BNAVY(X)  "\e[01;36m" X "\e[00m"
#define BLACK(X) "\e[00;30m" X "\e[00m"
#define RED(X)   "\e[00;31m" X "\e[00m"
#define GREEN(X) "\e[00;32m" X "\e[00m"
#define YELLOW(X)"\e[00;33m" X "\e[00m"
#define BLUE(X)  "\e[00;34m" X "\e[00m"
#define VIOLET(X)"\e[00;35m" X "\e[00m"
#define NAVY(X)  "\e[00;36m" X "\e[00m"

    std::vector<ssp_simple_t> allprimes;
    std::vector<int> positions;
    unsigned char * S;
    size_t B;
    int logI;
    int logA;

    bench_base(int logB, int logI, int logA) : logI(logI), logA(logA) {
        B = 1UL << logB;
        int const rc = posix_memalign((void**)&S, B, B);
        ASSERT_ALWAYS(rc == 0);
    }
    bench_base(bench_base const &) = delete;
    ~bench_base() { free(S); }
    bool test(candidate_list const & cand, const char * pfx="", std::string const& goal="current selection") {
        if (consistency_check_mode) {
            return test_correctness(cand);
        }
        if (cand.empty()) return true;
        std::vector<unsigned char *> refS;
        size_t const I = 1UL << logI;
        int const Nmax = 1 << (logA - LOG_BUCKET_REGION);

        std::vector<double> timings;
        int sel_index = -1;

        where_am_I w;
        for(auto const& bf : cand) {
            if (bf.sel) sel_index = &bf-cand.data();

            /* Do it once just to warm up */
            positions.clear();
            for(auto const & ssp : allprimes)
                positions.push_back((I/2)%ssp.get_p());
            memset(S, 0, B);
            for(int N = 0 ; N < Nmax ; N++) {
                (*bf.f)(positions, allprimes, S, logI, N, w);
            }

            /* and now with the intent of computing timings.  */
            positions.clear();
            for(auto const & ssp : allprimes)
                positions.push_back((I/2)%ssp.get_p());
            memset(S, 0, B);
            clock_t const tt = clock();
            for(int N = 0 ; N < Nmax ; N++) {
                (*bf.f)(positions, allprimes, S, logI, N, w);
            }

            timings.push_back((double) (clock()-tt) / CLOCKS_PER_SEC);
            printf(".");
            fflush(stdout);
            /* FIXME: we're only checking the very last region, here ! */
            if (!refS.empty()) {
                if (memcmp(S, refS.front(), B) != 0) {
                    fprintf(stderr, "inconsistency between %s and %s\n",
                            bf.name, cand.front().name);
                    for(size_t i = 0, n = 0 ; i < B && n < 16 ; i++) {
                        if (S[i] != refS.front()[i]) {
                            fprintf(stderr, "%04x: %02x != %02x\n",
                                    (unsigned int) i,
                                    (unsigned int) S[i],
                                    (unsigned int) refS.front()[i]);
                            n++;
                        }
                    }
                    if (abort_on_fail) abort();
                }
            }
            unsigned char * Scopy;
            int const rc = posix_memalign((void**)&Scopy, B, B);
            ASSERT_ALWAYS(rc == 0);
            memcpy(Scopy, S, B);
            refS.push_back(Scopy);
        }
        printf("\n");

        size_t const best_index = std::ranges::min_element(timings) - timings.begin();
        double const best_time = timings[best_index];
        if (sel_index >= 0) {
            bool const sel_best = (size_t) sel_index == best_index;
            bool const sel_tied = timings[sel_index] <= 1.05 * best_time;;
            for(auto const& bf : cand) {
                size_t const index = &bf-cand.data();
                double const tt = timings[index];
                bool const is_best = index == best_index;
                bool const near_best = tt <= 1.05 * best_time;
                printf("%s%-26s:\t%.3f\t", pfx, bf.name, tt);
                if (bf.sel) {
                    if (sel_best) {
                        printf(BGREEN("%s") "\n", goal.c_str());
                    } else if (sel_tied) {
                        printf(GREEN("%s, tied") "\n", goal.c_str());
                    } else {
                        printf(BRED("%s, wrong (+%.1f%%)?") "\n",
                                goal.c_str(),
                                100.0*(timings[sel_index]/best_time-1)
                                );
                    }
                } else if (is_best && !sel_best) {
                    if (sel_tied) {
                         printf(GREEN("current best, tied") "\n");
                    } else {
                        printf(BRED("current best") "\n");
                    }
                } else if (near_best) {
                    printf(GREEN("tied")"\n");
                } else {
                    printf("\n");
                }
            }
        } else {
            for(auto const& bf : cand) {
                size_t const index = &bf-cand.data();
                double const tt = timings[index];
                bool const is_best = index == best_index;
                bool const near_best = tt <= 1.05 * best_time;
                printf("%s%-24s:\t%.3f\t", pfx, bf.name, tt);
                if (is_best) {
                    printf(GREEN("current best") "\n");
                } else if (near_best) {
                    printf(GREEN("tied")"\n");
                } else {
                    printf("\n");
                }
            }
        }

        for(unsigned char * Sp : refS)
            free(Sp);
        return true;
    }

    /* This one has the loop order reversed. We're no longer testing
     * speed, we're checking correctness.
     */
    bool test_correctness(candidate_list const & cand) {
        if (cand.empty()) return true;
        if (cand.size()==1) {
            fprintf(stderr, "warning, %s is the only function we have, there's nothing to check againts...\n", cand.front().name);
            return true;
        }
        size_t const I = 1UL << logI;
        int const Nmax = 1 << (logA - LOG_BUCKET_REGION);

        positions.clear();
        for(auto const & ssp : allprimes)
            positions.push_back((I/2)%ssp.get_p());

        std::vector<unsigned char *> refS;
        std::vector<std::vector<int>> refpos;
        for(size_t i = 0 ; i < cand.size() ; i++) {
            unsigned char * Scopy;
            int const rc = posix_memalign((void**)&Scopy, B, B);
            ASSERT_ALWAYS(rc == 0);
            refS.push_back(Scopy);
            /* make a copy */
            refpos.push_back(positions);
        }
        printf("\n");
        if (!quiet) printf("Testing all %d buckets (logA=%d, logB=%d, logI=%d) for functions:\n",
                Nmax, logA, LOG_BUCKET_REGION, logI);
        for(auto const& bf : cand)
            if (!quiet) printf("  %s\n", bf.name);
        bool ok=true;
        std::vector<bool> ok_perfunc (cand.size(), true);
        int ndisp = 1;
        where_am_I w;
        for(int N = 0 ; N < Nmax ; N++) {
            for(auto const& bf : cand) {
                size_t const index = &bf-cand.data();
                memset(refS[index], 0, B);
                (*bf.f)(refpos[index], allprimes, refS[index], logI, N, w);
            }
            for(auto const& bf : cand) {
                size_t const index = &bf-cand.data();
                if (!ok_perfunc[index]) continue;
                if (index > 0) {
                    if (memcmp(refS[index], refS.front(), B) != 0) {
                        ok = ok_perfunc[index] = false;
                        fprintf(stderr, "bucket %d: inconsistency between %s and %s\n",
                                N, bf.name, cand.front().name);
                        for(size_t i = 0, n = 0 ; i < B && n < 16 ; i++) {
                            if (refS[index][i] != refS.front()[i]) {
                                fprintf(stderr, "%04x: %02x != %02x\n",
                                        (unsigned int) i,
                                        (unsigned int) refS[index][i],
                                        (unsigned int) refS.front()[i]);
                                n++;
                            }
                        }
                    }
                    if (refpos[index] != refpos.front()) {
                        ok = ok_perfunc[index] = false;
                        fprintf(stderr, "bucket %d: inconsistency in post-sieve positions between %s and %s\n",
                                N, bf.name, cand.front().name);
                        for(size_t i = 0, n = 0 ; i < allprimes.size() && n < 16 ; i++) {
                            if (refpos[index][i] != refpos.front()[i]) {
                                fprintf(stderr, "p=%u r=%u ; %d!=%d\n",
                                        allprimes[i].get_p(),
                                        allprimes[i].get_r(),
                                        refpos[index][i],
                                        refpos.front()[i]);
                                n++;
                            }
                        }
                    }
                }
            }
            if (abort_on_fail && !ok) abort();
            if (!quiet) for( ; (N+1)*16 >= (ndisp * Nmax) ; ndisp++) {
                printf(".");
                fflush(stdout);
            }
        }
        if (!quiet) printf("\n");
        for(unsigned char * Sp : refS)
            free(Sp);
        return ok;
    }
};

static void store_primes(std::vector<ssp_simple_t>& allprimes, int bmin, int bmax, gmp_randstate_t rstate)
{
    mpz_t pz;
    mpz_init(pz);
    mpz_set_ui(pz, 1u << bmin);
    for(;;) {
        mpz_nextprime(pz, pz);
        if (mpz_cmp_ui(pz, 1u << bmax) >= 0) break;
        unsigned long const p = mpz_get_ui(pz);
        unsigned long const r = gmp_urandomm_ui(rstate, p);
        /* not a problem if we put garbage here */
        unsigned char const logp = log2(p)/log2(1.6);
        allprimes.emplace_back(p, r, logp);
    }
    if (!quiet) printf("created a list of %zu primes\n", allprimes.size());
    mpz_clear(pz);
}

static void declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "q",  "quiet mode (for tests, mostly)");
    param_list_decl_usage(pl, "C",  "run tests, not timings");
    param_list_decl_usage(pl, "F",  "abort on test failure");
    param_list_decl_usage(pl, "I",  "set logI (default = 16)");
    param_list_decl_usage(pl, "A",  "set logA (default = 2*logI-1)");
#ifndef LOG_BUCKET_REGION_IS_A_CONSTANT
    param_list_decl_usage(pl, "B",  "set LOG_BUCKET_REGION (default = 16)");
#endif
    param_list_decl_usage(pl, "bmin",  "restrict test or timings to primes >= 2^bmin (i.e. (bmin+1)-bit primes)");
    param_list_decl_usage(pl, "bmax",  "restrict test or timings to primes < 2^bmax (i.e. up to bmax-bit primes)");
    param_list_decl_usage(pl, "only-complete-functions",  "restrict to testing the complete small sieve functions");
}

int main(int argc0, char const * argv0[])
{
    int logI = 16;
    int logA = 0;
    int bmin = 0;
    int bmax = 0;
    int argc = argc0;
    char const **argv = argv0;
    cxx_param_list pl;

    declare_usage(pl);
    param_list_configure_switch(pl, "-q", &quiet);
    param_list_configure_switch(pl, "-C", &consistency_check_mode);
    param_list_configure_switch(pl, "-F", &abort_on_fail);
    param_list_configure_switch(pl, "--only-complete-functions", &only_complete_functions);
    argv++, argc--;
    for( ; argc ; ) {
        if (param_list_update_cmdline(pl, &argc, &argv)) { continue; }
        fprintf(stderr, "Unhandled parameter %s\n", argv[0]);
        param_list_print_usage(pl, argv0[0], stderr);
        exit(EXIT_FAILURE);
    }


    param_list_parse_int(pl, "I", &logI);
    param_list_parse_int(pl, "A", &logA);
#ifndef LOG_BUCKET_REGION_IS_A_CONSTANT
    param_list_parse_int(pl, "B", &LOG_BUCKET_REGION);
#endif
    param_list_parse_int(pl, "bmin", &bmin);
    param_list_parse_int(pl, "bmax", &bmax);

    if (!bmax) bmax = logI;
    if (!logA) logA = 2*logI-1;

    gmp_randstate_t rstate;
    gmp_randinit_default(rstate);

    int errors = 0;

    if (!only_complete_functions) {
        // std::vector<std::pair<size_t, size_t>> bounds_perbit(logI + 1, {0,0});
        for(int b = bmin ; b < bmax ; b++) {
            if (b == 0) continue;
#if 0
            /* primes within [1<<b, 1<<(b+1)] are (b+1)-bit primes. With respect to the
             * semantics that we use in the template code, those are counted as
             * (logI-bit) off the bound. So they will use the code in
             * all_candidates_for_evenline<logI-bit> and
             * all_candidates_for_oddline<logI-bit>.
             */
            size_t i0 = 0;
            for( ; i0 < allprimes.size() && !(allprimes[i0].get_p() >> b) ; i0++);
            size_t i1 = i0;
            for( ; i1 < allprimes.size() && (allprimes[i1].get_p() >> b) ; i1++);
            bounds_perbit[b] = { i0, i1 };
#endif
            int const logfill = std::min(logI, LOG_BUCKET_REGION);

            printf("===== now doing specific tests for %d-bit primes using candidates<%d> =====\n", b+1, logfill-b);
            bench_base bbase(LOG_BUCKET_REGION, logI, logA);
            // std::vector<ssp_simple_t>& allprimes(bbase.allprimes);
            // std::vector<int>& positions(bbase.positions);
            // unsigned char * & S (bbase.S);
            store_primes(bbase.allprimes, b, b + 1, rstate);


            // std::vector<ssp_simple_t> lp(allprimes.begin() + i0, allprimes.begin() + i1);
            // std::vector<int>& lx(positions.begin() + i0, positions.begin() + i1);
            candidate_list cand1, cand2;

            if (logI > LOG_BUCKET_REGION) {
                // (true!) bool fragment = logI > LOG_BUCKET_REGION;

                switch(logfill-b) {
#define CASE(xxx)						\
                    case xxx:					\
                            factory_for_bit_round1<xxx, true>()(cand1);	\
                            factory_for_bit_round2<xxx, true>()(cand2);	\
                    break
                    CASE(1); CASE(2); CASE(3); CASE(4);
                    CASE(5); CASE(6); CASE(7); CASE(8);
                    CASE(9); CASE(10); CASE(11); CASE(12);
                    CASE(13); CASE(14); CASE(15);
#undef CASE
                    default:
                    if (logfill <= b) {
                        factory_for_bit_round1<0, true>()(cand1);
                        factory_for_bit_round2<0, true>()(cand2);
                    } else {
                        fprintf(stderr, "sorry, not handled\n");
                    }
                }
            } else {
            switch(logfill-b) {
#define CASE(xxx)						\
                case xxx:					\
                            factory_for_bit_round1<xxx, false>()(cand1);	\
                            factory_for_bit_round2<xxx, false>()(cand2);	\
                break
                CASE(1); CASE(2); CASE(3); CASE(4);
                CASE(5); CASE(6); CASE(7); CASE(8);
                CASE(9); CASE(10); CASE(11); CASE(12);
                CASE(13); CASE(14); CASE(15);
#undef CASE
                default:
                if (logfill <= b) {
                        factory_for_bit_round1<0, false>()(cand1);
                        factory_for_bit_round2<0, false>()(cand2);
                } else {
                    fprintf(stderr, "sorry, not handled\n");
                }
            }
            }

            std::vector<unsigned char *> const refS;

            if (!cand1.empty()) {
                std::ostringstream goal_name;
                goal_name << "best_oddline<" << logfill-b << ">";
                printf("  testing odd lines ");
                errors += !bbase.test(cand1, "    ", goal_name.str());
            }
            if (!cand2.empty()) {
                std::ostringstream goal_name;
                goal_name << "best_evenline<" << logfill-b << ">";
                printf("  testing even lines ");
                errors += !bbase.test(cand2, "    ", goal_name.str());
            }
        }
    }

    printf("===== now doing tests for complete functions =====\n");
    {
        bench_base bbase(LOG_BUCKET_REGION, logI, logA);
        // std::vector<ssp_simple_t>& allprimes(bbase.allprimes);
        // std::vector<int>& positions(bbase.positions);
        // unsigned char * & S (bbase.S);
        store_primes(bbase.allprimes, bmin, bmax, rstate);

        candidate_list funcs {
            // modified_I18_branch_C,
                { false, current_I18_branch, "I18" },
                { false, devel_branch, "devel" },
                { false, generated, "devel-meta" },
        };

        if (logI <= LOG_BUCKET_REGION) funcs.emplace_back(false, legacy_branch, "legacy");
        funcs.emplace_back(false, legacy_mod_branch, "legacy_mod");

        errors += !bbase.test(funcs);
    }

    gmp_randclear(rstate);

    return errors ? EXIT_FAILURE : EXIT_SUCCESS;
}

