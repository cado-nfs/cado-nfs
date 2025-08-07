#include "cado.h" // IWYU pragma: keep
#include <climits>

#include <algorithm>

#ifndef LINGEN_BINARY
#include <gmp.h>
#endif
#include "lingen_bw_dimensions.hpp"
#include "lingen_expected_pi_length.hpp"
#include "lingen_matpoly_select.hpp"
#include "macros.h"


std::tuple<unsigned int, unsigned int> get_minmax_delta(std::vector<unsigned int> const & delta)/*{{{*/
{
    unsigned int maxdelta = 0;
    unsigned int mindelta = UINT_MAX;
    for(auto x : delta) {
        x = std::max(x, maxdelta);
        x = std::min(x, maxdelta);
    }
    return std::make_tuple(mindelta, maxdelta);
}/*}}}*/

#if 0
static unsigned int get_min_delta(std::vector<unsigned int> const & delta)/*{{{*/
{
    unsigned int mindelta, maxdelta;
    std::tie(mindelta, maxdelta) = get_minmax_delta(delta);
    return mindelta;
}/*}}}*/
static unsigned int get_max_delta(std::vector<unsigned int> const & delta)/*{{{*/
{
    unsigned int mindelta, maxdelta;
    std::tie(mindelta, maxdelta) = get_minmax_delta(delta);
    return maxdelta;
}/*}}}*/
#endif

template<bool is_binary>
static unsigned int
sizeinbase(typename matpoly<is_binary>::arith_hard * ab);

#ifndef LINGEN_BINARY
template<>
unsigned int
sizeinbase<false>(matpoly<false>::arith_hard * ab)
{
    mpz_srcptr p = ab->characteristic();
    if (mpz_cmp_ui(p, 1024) >= 0) {
        return mpz_sizeinbase(p, 2);
        // l *= ab->degree();    /* roughly log_2(#K) */
    } else {
        // mpz_pow_ui(p, p, ab->degree());
        return mpz_sizeinbase(p, 2);
    }
}
#endif

#ifdef LINGEN_BINARY
template<>
unsigned int
sizeinbase<true>(matpoly<true>::arith_hard *) {
    // K is GF(2), period.
    return 1;
}
#endif

template<bool is_binary>
unsigned int expected_pi_length(bw_dimensions<is_binary> & d, unsigned int len)/*{{{*/
{
    /* The idea is that we want something which may account for something
     * exceptional, bounded by probability 2^-64. This corresponds to a
     * column in e (matrix of size m*b) to be spontaneously equal to
     * zero. This happens with probability (#K)^-m.
     * The solution to
     * (#K)^(-m*x) > 2^-64
     * is m*x*log_2(#K) < 64
     *
     * We thus need to get an idea of the value of log_2(#K).
     *
     * (Note that we know that #K^abgroupsize(ab) < 2^64, but that bound
     * might be very gross).
     *
     * The same esitmate can be used to appreciate what we mean by
     * ``luck'' in the end. If a column happens to be zero more than
     * expected_pi_length(d,0) times in a row, then the cause must be
     * more than sheer luck, and we use it to detect generating rows.
     */

    unsigned int const m = d.m;
    unsigned int const n = d.n;
    unsigned int const b = m + n;
    auto * ab MAYBE_UNUSED = & d.ab;
    unsigned int const res = 1 + iceildiv(len * m, b);
    unsigned int const l = sizeinbase<is_binary>(ab);
    // unsigned int safety = iceildiv(abgroupsize(ab), m * sizeof(abelt));
    unsigned int safety = iceildiv(64, m * l);
    // this +1 is here because I've sometimes seen the c30 fail in the
    // lingen step because of the pi size check. I wonder whether I
    // should have a +1 here, or maybe a +t0. Anyway. Quite sure that
    // it's in the whereabouts of adding a small offset.
    safety++;
    return res + safety;
}/*}}}*/

template<bool is_binary>
unsigned int expected_pi_length(bw_dimensions<is_binary> & d, std::vector<unsigned int> const & delta, unsigned int len)/*{{{*/
{
    // see comment above.

    unsigned int mi, ma;
    std::tie(mi, ma) = get_minmax_delta(delta);

    return expected_pi_length(d, len) + ma - mi;
}/*}}}*/

template<bool is_binary>
unsigned int expected_pi_length_lowerbound(bw_dimensions<is_binary> & d, unsigned int len)/*{{{*/
{
    /* generically we expect that len*m % (m+n) columns have length
     * 1+\lfloor(len*m/(m+n))\rfloor, and the others have length one more.
     * For one column to have a length less than \lfloor(len*m/(m+n))\rfloor,
     * it takes probability 2^-(m*l) using the notations above. Therefore
     * we can simply count 2^(64-m*l) accidental zero cancellations at
     * most below the bound.
     * In particular, it is sufficient to derive from the code above!
     */
    unsigned int const m = d.m;
    unsigned int const n = d.n;
    unsigned int const b = m + n;
    auto * ab MAYBE_UNUSED = & d.ab;
    unsigned int const res = 1 + (len * m) / b;
    unsigned int const l = sizeinbase<is_binary>(ab);
    unsigned int const safety = iceildiv(64, m * l);
    return safety < res ? res - safety : 0;
}/*}}}*/

#ifdef LINGEN_BINARY
template
unsigned int expected_pi_length(bw_dimensions<true> & d, unsigned int len);
template
unsigned int expected_pi_length(bw_dimensions<true> & d, std::vector<unsigned int> const & delta, unsigned int len);
template
unsigned int expected_pi_length_lowerbound<true>(bw_dimensions<true> & d, unsigned int len);
#else
template
unsigned int expected_pi_length(bw_dimensions<false> & d, unsigned int len);
template
unsigned int expected_pi_length(bw_dimensions<false> & d, std::vector<unsigned int> const & delta, unsigned int len);
template
unsigned int expected_pi_length_lowerbound<false>(bw_dimensions<false> & d, unsigned int len);
#endif /* LINGEN_BINARY */
