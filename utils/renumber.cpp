#include "cado.h" // IWYU pragma: keep

#include <cstdint>
#include <cstdio>
#include <climits>
#include <cmath>

#include <algorithm>
#include <array>
#include <bit>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <sys/stat.h>

#include <gmp.h>
#include "fmt/format.h"

#include "badideals.hpp"
#include "cxx_mpz.hpp"
#include "misc.h"
#include "mmap_allocator.hpp"
#include "getprime.h"
#include "gmp_aux.h"
#include "fstream_maybe_compressed.hpp"
#include "arith/mod_ul.h"
#include "mpz_poly.h"
#include "omp_proxy.h"
#include "params.hpp"
#include "renumber.hpp"
#include "rootfinder.h"
#include "stats.h"
#include "macros.h"
#include "typedefs.h"

#if defined(__GLIBCXX__) && defined(_GLIBCXX_DEBUG) && defined(_GLIBCXX_DEBUG_DISABLE_CHECK_PARTITIONED)
namespace __gnu_debug {
template<>
inline bool __check_partitioned_lower<std::vector<std::array<p_r_values_t, 2>>::const_iterator>(
        std::vector<std::array<p_r_values_t, 2>>::const_iterator,
        std::vector<std::array<p_r_values_t, 2>>::const_iterator,
        std::vector<std::array<p_r_values_t, 2>>::value_type const &)
{
    /* This is really ugly, but unfortunately when _GLIBCXX_DEBUG is
     * enabled, we're ending up with O(n) checks in
     * __check_partitioned_lower whenever we do a dichotomic search
     * through the renumber table. Which is of course horribly expensive,
     * especially so when our test code does O(n) dichotomic searches by
     * itself...
     *
     * The cure can be one of:
     *  - do not use _GLIBCXX_DEBUG (it's currently _not_ used in our
     *  tests)
     *  - use it, but define _GLIBCXX_DEBUG_DISABLE_CHECK_PARTITIONED
     *  along with it.
     */
    return true;
}
}
#endif

/* Some documentation on the internal encoding of the renumber table...
 *
 * XXX This format has changed incompatibly in 2020 XXX
 *
 * The goal is to have a table that converts to/from two formats:
 *
 *  - an integer index
 *
 *  - a triple (side, prime number p, root in [0..p]) that represents a
 *    prime ideal.
 *    (yes, p included -- that means a projective root, i.e. an ideal
 *    above p that divides J=<1,alpha>^-1).
 *
 * We want to minimize the storage, and guarantee cheap lookups in both
 * cases.
 */

/* {{{ As of 20200515, the known uses of this table in cado-nfs, taking out
 * the uses that pertain only to the building and handling of the
 * database itself, are the following.
 *
 * index_from_p_r :
 *      filter/dup2.cpp
 *      sieve/las-dlog-base.cpp
 *      sieve/fake_rels.cpp
 * 
 * indices_from_p_a_b :
 *      filter/dup2.cpp
 *      sieve/las-dlog-base.cpp  -- not yet, but maybe at some point.
 * 
 * p_r_from_index but as a forward iterator only
 *      utils/galois_action.cpp
 *      filter/reconstructlog.cpp
 *      sieve/fake_rels.cpp
 *      filter/reconstructlog.cpp
 *
 * p_r_from_index, random access:
 *      misc/convert_rels.cpp
 *
 * is_bad :
 *      filter/dup2.cpp
 *      sieve/las-dlog-base.cpp
 *      sieve/fake_rels.cpp
 *      (there, we do something very simple, and only pick a random index for a bad ideal, in [index, index+nb) )
 *
 * is_additional_column :
 *      filter/reconstructlog.cpp
 *      sieve/fake_rels.cpp
 *
 * }}} */

/* In the above uses, the lookup (side,p,r) -> index is critical in the
 * dup2 and descent steps, while the lookup index -> (side,p,r) is most
 * often used only for forward iterators, where the context can be used
 * to carry some extra information.
 */

/* There are various strategies that can be used to have an in-memory
 * structure that makes these lookups possible.
 */

/* My personal preference would be to bite the bullet and allow for two
 * integers per ideal, always storing both p and r explicitly:
 *  - The _from_index lookups would be trivial (we could keep the encoding
 *    of r as vr = side'*(p+1)+r).
 *  - The _from_p_r lookups would also be doable with very simple-minded
 *    dichotomy (even the standard library can be put to some use here)
 *  - The database in itself would be human readable.
 */
/*}}}*/

/* {{{ exceptions */
renumber_t::corrupted_table::corrupted_table(std::string const & x)
    : std::runtime_error(std::string("Renumber table is corrupt: ") + x)
{}

static renumber_t::corrupted_table wrong_entry(p_r_values_t p, p_r_values_t vr) {/*{{{*/
    auto msg = fmt::format("above p=0x{:x}, the index 0x{:x} makes no sense",
            p, vr);
    return renumber_t::corrupted_table(msg);
}/*}}}*/
static renumber_t::corrupted_table prime_is_too_large(p_r_values_t p) {/*{{{*/
    auto msg = fmt::format("prime 0x{:x} is too large!", p);
#if SIZEOF_INDEX != 8
    msg += "\n"
        "Note: above 2^32 ideals or relations,"
        " add FLAGS_SIZE=\"-DSIZEOF_P_R_VALUES=8 -DSIZEOF_INDEX=8\""
        " to local.sh\n";
#endif
    return renumber_t::corrupted_table(msg);
}/*}}}*/
static renumber_t::corrupted_table prime_maps_to_garbage(int format, p_r_values_t p, index_t i, p_r_values_t q)/*{{{*/
{
    auto msg = fmt::format("cached index for prime p=0x{:x} is 0x{:x},"
            " which points to q=0x{:x} (should be {}) ;"
            " note: isprime(p)=={}",
            p, i, q,
            format == renumber_t::format_flat ? "p" : "vp",
            ulong_isprime(p));
    return renumber_t::corrupted_table(msg);
}/*}}}*/
static renumber_t::corrupted_table prime_maps_to_garbage(int format, p_r_values_t p, index_t i)/*{{{*/
{
    auto msg = fmt::format("cached index for prime p=0x{:x} is 0x{:x},"
            " which points nowhere (should be {}) ;"
            " note: isprime(p)=={}",
            p, i,
            format == renumber_t::format_flat ? "p" : "vp",
            ulong_isprime(p));
    return renumber_t::corrupted_table(msg);
}/*}}}*/
static renumber_t::corrupted_table cannot_find_pr(renumber_t::p_r_side x)/*{{{*/
{
    auto msg = fmt::format("cannot find p=0x{:x}, r=0x{:x} on side {}",
            x.p, x.r, x.side);
    /* #30012 and #30048 are cases where we see composites. #30048 in
     * particular sees the code crashing here, which is really cryptic.
     * We can alleviate that with a simple test.
     */
    if (!ulong_isprime(x.p)) {
        msg += " ; NOTE that this p is not prime"
            ", which is very unexpected. See bug #30048";
    }
    return renumber_t::corrupted_table(msg);
}/*}}}*/
static renumber_t::corrupted_table cannot_lookup_p_a_b_in_bad_ideals(renumber_t::p_r_side x, int64_t a, uint64_t b)/*{{{*/
{
    auto msg = fmt::format(
            "failed bad ideal lookup for (a,b)=({},{})"
            " at p=0x{:x} on side {}",
            a, b, x.p, x.side);
    return renumber_t::corrupted_table(msg);
}/*}}}*/
static renumber_t::corrupted_table cannot_lookup_p_a_b_in_bad_ideals(renumber_t::p_r_side x, mpz_srcptr a, mpz_srcptr b)/*{{{*/
{
    auto msg = fmt::format(
            "failed bad ideal lookup for (a,b)=({},{})"
            " at p=0x{:x} on side {}",
            cxx_mpz(a), cxx_mpz(b), x.p, x.side);
    return renumber_t::corrupted_table(msg);
}/*}}}*/
static renumber_t::corrupted_table parse_error(std::string const & what)/*{{{*/
{
    return renumber_t::corrupted_table(fmt::format("parse error ({})", what));
}/*}}}*/
/* }}} */

/* {{{ helper functions relative to the "traditional" and "variant"
 * formats
 */

/* This one is a helper only, because we use it for both T ==
 * p_r_values_t and T == uint64_t */
template<typename T>
static inline T vp_from_p(T p, int n, int c, int format MAYBE_UNUSED)
{
    /* The final "+d" is not necessary, but we keep it for compatibility */
    int const d = 0;
    return (n - c) * (p + 1) + d;
}

/* in the "usual" case, we have vp==p+1, but that is no longer true when
 * we have two algebraic sides, or when we have multiple sides.
 */
p_r_values_t renumber_t::compute_vp_from_p (p_r_values_t p) const/*{{{*/
{
    int const n = get_nb_polys();
    int const c = (get_rational_side() >= 0);
    return vp_from_p(p, n, c, format);
}/*}}}*/


p_r_values_t renumber_t::compute_vr_from_p_r_side (renumber_t::p_r_side x) const/*{{{*/
{
    if (x.side == get_rational_side()) {
        /* the rational root _always_ has r encoded implicitly. In truth,
         * we rarely need to look it up, except that we have to decide on
         * something to store in the table...
         */
        return compute_vp_from_p(x.p);
    }
    p_r_values_t vr = x.side * (x.p + 1) + x.r;
    if (get_rational_side() >= 0 && x.side > get_rational_side())
        vr -= x.p + 1;
    return vr;
}/*}}}*/

renumber_t::p_r_side renumber_t::compute_p_r_side_from_p_vr (p_r_values_t p, p_r_values_t vr) const/*{{{*/
{
    /* Note that vr is only used for the encoding of non-rational ideals.
     */
    p_r_side res { .p=p, .r=0, .side=0 };

    res.r = vr;
    for(res.side = 0 ; res.side < get_nb_polys() ; res.side++) {
        if (res.side == get_rational_side())
            continue;
        if (res.r <= p)
            return res;
        res.r -= p + 1;
    }
    if (res.r == 0) {
        if (get_rational_side() >= 0) {
            res.side = get_rational_side();
            return res;
        }
    }
    throw wrong_entry(p, vr);
}/*}}}*/
/* }}} */

/* sort in decreasing order. Faster than qsort for ~ < 15 values in r[] */

void renumber_t::cook_into(unsigned long p,
        std::vector<std::vector<unsigned long>> & roots,
        fragment & F) const
{
    size_t total_nroots = 0;

    /* Note that roots always has a root on the rational side, even
     * though it's only a zero -- the root itself isn't computed.
     */
    for (int i = 0; i < get_nb_polys() ; i++)
        total_nroots += roots[i].size();

    /* primes with no ideal above them leave no trace at all */
    if (total_nroots == 0) return;

    for (int i = 0; i < get_nb_polys() ; i++) {
        ASSERT_ALWAYS(roots[i].size() <= UCHAR_MAX);
        F.nroots.push_back((uint8_t) roots[i].size());
    }

    for (int side = 0 ; side < get_nb_polys(); side++) {
        /* reverse the ordering of the *ROOTS* (not of the sides), because
         * our goal is to remain compatible with the old-format indexing
         */
        for (auto it = roots[side].rbegin() ; it != roots[side].rend() ; ++it) {
            p_r_side const x { (p_r_values_t) p, (p_r_values_t) *it, side };
            F.flat.emplace_back(
                    std::array<p_r_values_t, 2> {{
                    (p_r_values_t) p,
                    compute_vr_from_p_r_side (x)
                    }});
        }
    }
}

/* return the number of bad ideals above x (and therefore zero if
 * x is not bad). If the ideal is bad, put in the reference [first] the
 * first index that corresponds to the bad ideals.
 */
int renumber_t::is_bad (index_t & first, p_r_side x) const
{
    if (x.p > bad_ideals_max_p) return 0;

    /* bad ideals start after the additional columns in the table. */
    first = above_add;

    for (auto const & I : bad_ideals) {
        if (x == I.first)
            return I.second.nbad;
        first += I.second.nbad;
    }
    return 0;
}

int renumber_t::is_bad (p_r_side x) const
{
    index_t first;
    return is_bad(first, x);
}

int renumber_t::is_bad (index_t & first, index_t i) const
{
    if (i >= above_bad)
        return 0;
    i -= above_add;
    first = above_add;
    for (auto const & I : bad_ideals) {
        index_t const n = I.second.nbad;
        if (i < n) return n;
        i -= n;
        first += n;
    }
    throw corrupted_table("bad bad ideals");
}

/* XXX This is an important part of the index_from lookup. This one needs
 * to do a dichotomy, one way or another.
 *
 * Note that we return the index relative to the internal table, which is
 * shifted by [[above_bad]] compared to the indices that we return in the
 * public interface.
 *
 * if p is a non-prime, we return a lower bound for the first entry that
 * is >= p.
 */
auto renumber_t::get_first_iterator_from_p(p_r_values_t p) const -> const_iterator
{
    // p_r_values_t side = x.side;
    if (p < index_from_p_cache.size()) {
        index_t i;
        /* Do this loop so that we can safely return something that is
         * valid even if p is not a prime.
         */
        for( ; ; p++) {
            if (p >= index_from_p_cache.size())
                return get_first_iterator_from_p(p);
            i = index_from_p_cache[p];
            if (i != std::numeric_limits<index_t>::max())
                break;
        }
        if (UNLIKELY((above_bad + i) >= size()))
            throw prime_maps_to_garbage(format, p, i);
        if (format == format_flat) {
            if (UNLIKELY(flat_data[i][0] != p))
                throw prime_maps_to_garbage(format, p, i, flat_data[i][0]);
        }
        index_t const outer_idx = i + above_bad;
        return { *this, outer_idx };
    }

    /* Note that if lpbmax is > the width of the type, then
     * check_needed_bits() should have complained already.
     *
     * If lpbmax is == the width of the type, then the table fits (if
     * there's only one non-rational side). But we can't simply check p
     * >> lpbmax, since the check overflows...
     *
     */
    unsigned int const lpbmax = *std::max_element(lpb.begin(), lpb.end());
    if (lpbmax < (8 * sizeof(p_r_values_t)) && UNLIKELY(p >> lpbmax))
        throw prime_is_too_large(p);

    std::array<p_r_values_t, 2> const p0 {{ p, 0 }};
    {
        auto it = std::lower_bound(flat_data.begin(), flat_data.end(), p0);
        if (it == flat_data.end())
            throw prime_is_too_large(p);
        if ((*it)[0] != p) {
            /* we can reach here if p is not prime and we simply want a
             * good lower bound on the entries that are >= p. For this,
             * std::lower_bound really is the thing that we want.
             *
             * Here, (*it)[0] is a prime larger than p. But by the
             * definition of lower_bound, (*--it)[0] is a prime that
             * is less than p, and we definitely don't want to return
             * that. So the good answer is really the thing that is
             * pointed to by the iterator it.
             */
            //throw prime_maps_to_garbage(format, p, i, (*it)[0]);
        }
        index_t const outer_idx = it - flat_data.begin() + above_bad;
        return { *this, outer_idx };
    }
}

index_t renumber_t::index_from_p_r (p_r_side x) const
{
    index_t i;
    if (is_bad(i, x))
        return i;
    if (x.p == 0) {
        // additional columns. By convention they're attached to p==0.
        // Note that also by convention/tradition we use only a single
        // additional column when we have two non-monic sides, in which
        // case it's only counted on side 0.
        if (has_merged_additional_column()) {
            if (x.side == 0)
                return 0;
            throw cannot_find_pr(x);
        }
        i = 0;
        for(auto side : get_sides_of_additional_columns()) {
            if (side == x.side)
                return i;
            i++;
        }
        throw cannot_find_pr(x);
    }
    p_r_values_t const vr = compute_vr_from_p_r_side(x);

    for(auto it = get_first_iterator_from_p(x.p) ; it.raw()[0] == x.p ; ++it) {
        if (it.raw()[1] == vr)
            return it.i;
    }
    throw cannot_find_pr(x);
}

/* This returns the smallest outer index i0 that stores a prime ideal above a
 * prime >=p0 ; note that p0 does not have to be prime.
 * This function may return get_max_index() if no prime >= p0 is in the
 * table.
 */
index_t renumber_t::index_from_p(p_r_values_t p0) const
{
    return iterator_from_p(p0).i;
}
index_t renumber_t::index_from_p(p_r_values_t p0, int side) const
{
    return iterator_from_p(p0, side).i;
}
    
renumber_t::const_iterator renumber_t::iterator_from_p(p_r_values_t p0) const
{
    if (p0 >> get_max_lpb())
        return end();
    const_iterator it = begin();
    if (p0 >= bad_ideals_max_p)
        it = get_first_iterator_from_p(p0);
    for( ; it != end() && (*it).p < p0 ; ++it);
    return it;
}

renumber_t::const_iterator renumber_t::iterator_from_p(p_r_values_t p0, int side) const
{
    const_iterator it = iterator_from_p(p0);
    for( ; it != end() && (*it).side != side ; ++it);
    return it;
}
    
/* This used to be handle_bad_ideals in filter/filter_badideals.cpp ; in
 * fact, this really belongs here.
 */
std::pair<index_t, std::vector<int>> renumber_t::indices_from_p_a_b(p_r_side x, int e, int64_t a, uint64_t b) const
{
    /* We want this interface to be called for bad ideals only.
     *
     * A priori we don't need x.r since it can be reocmputed. However,
     * since in all cases the caller just computed this r a moment ago,
     * and it's used as a lookup key in the badideals data, let's use it.
     * x.p and x.side are of course necessary.
     *
     */

    /* bad ideals start after the additional columns in the table. */
    index_t first = above_add;
    std::vector<int> const exps;
    for (auto const & I : bad_ideals) {
        if (x == I.first) {
            std::vector<int> exps;

            /* work here, and return */
            for(auto J : I.second.branches) {
                int k = J.k;
                p_r_values_t pk = x.p;
                ASSERT_ALWAYS(pk);
                for( ; --k ; ) {
                    p_r_values_t const pk1 = pk * x.p;
                    /* we want to be sure that we don't wrap around ! */
                    ASSERT_ALWAYS(pk1 > pk);
                    pk = pk1;
                }
                p_r_values_t const rk = mpz_get_uint64(J.r);
                /* rk represents an element in P^1(Z/p^k). Ir rk < pk, it's
                 * (rk:1), otherwise it's (1:p*(rk-pk)) (I think)
                 */
                p_r_values_t uk = rk;
                p_r_values_t vk = 1;
                if (rk >= pk) {
                    uk = 1;
                    vk = rk - pk;
                }
                modulusul_t m;
                residueul_t ma, mb, muk, mvk;
                modul_initmod_ul (m, pk);
                modul_init (ma, m);
                modul_init (mb, m);
                modul_init (muk, m);
                modul_init (mvk, m);
                modul_set_int64 (ma, a, m);
                modul_set_uint64 (mb, b, m);
                modul_set_int64  (muk, uk, m);
                modul_set_uint64 (mvk, vk, m);
                modul_mul(ma, ma, mvk, m);
                modul_mul(mb, mb, muk, m);
                modul_sub(ma, ma, mb, m);
                if (modul_intcmp_ul(ma, 0) == 0) {
                    /* we found the right ideal ! */
                    for(int const v : J.v) {
                        if (v >= 0) {
                            exps.push_back(v);
                        } else {
                            ASSERT_ALWAYS(e >= -v);
                            exps.push_back(e + v);
                        }
                    }
                    return { first, exps };
                }
            }
            /* then it's a fatal error ! */
            break;
        }
        first += I.second.nbad;
    }
    throw cannot_lookup_p_a_b_in_bad_ideals(x, a, b);
}

/* Use the fact that this interface is for bad ideals only and that bad ideals
 * are small.
 * Reduce a and b module p^m where m is maximal power of p appearing in the
 * branches.
 */
std::pair<index_t, std::vector<int>> renumber_t::indices_from_p_a_b(
            p_r_side x,
            int e,
            mpz_srcptr a,
            mpz_srcptr b) const
{
    /* bad ideals start after the additional columns in the table. */
    int m = 1;
    p_r_values_t pm = x.p;
    ASSERT_ALWAYS(pm);

    for (auto const & I : bad_ideals) {
        if (x == I.first) {
            for(auto J : I.second.branches) {
                for( ; m < J.k; ++m) {
                    p_r_values_t const pm1 = pm * x.p;
                    /* we want to be sure that we don't wrap around ! */
                    ASSERT_ALWAYS(pm1 > pm);
                    pm = pm1;
                }
            }
            int64_t sa = (int64_t) mpz_tdiv_uint64(a, pm);
            if (mpz_sgn(a) < 0) {
                sa = -sa;
            }
            return indices_from_p_a_b(x, e, sa, mpz_tdiv_uint64(b, pm));
        }
    }
    throw cannot_lookup_p_a_b_in_bad_ideals(x, a, b);
}

/* additional columns _must_ be handled differently at this point */
renumber_t::p_r_side renumber_t::p_r_from_index (index_t i) const
{
    if (i < above_add) {
        /* In the special case where we have two non-monic polynomials,
         * we have a single additional column, hence above_add=1 and i=0.
         * We return {0,0,(arbitrarily the first of the two sides)} in
         * that case.
         */
        for(auto side : get_sides_of_additional_columns()) {
            if (i-- == 0)
                return { 0, 0, side };
        }
        throw corrupted_table("bad additional columns");
    }
    if (i < above_bad) {
        i -= above_add;
        for (auto const & I : bad_ideals) {
            if (i < (index_t) I.second.nbad)
                return I.first;
            i -= I.second.nbad;
        }
        throw corrupted_table("bad bad ideals");
    }
    i -= above_bad;

    p_r_values_t const p = flat_data[i][0];
    p_r_values_t const vr = flat_data[i][1];
    return compute_p_r_side_from_p_vr(p, vr);
}

static uint64_t previous_prime_of_powers_of_2[65] = { 0x0, 0x0, 0x3, 0x7, 0xd,
        0x1f, 0x3d, 0x7f, 0xfb, 0x1fd, 0x3fd, 0x7f7, 0xffd, 0x1fff, 0x3ffd,
        0x7fed, 0xfff1, 0x1ffff, 0x3fffb, 0x7ffff, 0xffffd, 0x1ffff7, 0x3ffffd,
        0x7ffff1, 0xfffffd, 0x1ffffd9, 0x3fffffb, 0x7ffffd9, 0xfffffc7,
        0x1ffffffd, 0x3fffffdd, 0x7fffffff, 0xfffffffb, 0x1fffffff7,
        0x3ffffffd7, 0x7ffffffe1, 0xffffffffb, 0x1fffffffe7, 0x3fffffffd3,
        0x7ffffffff9, 0xffffffffa9, 0x1ffffffffeb, 0x3fffffffff5, 0x7ffffffffc7,
        0xfffffffffef, 0x1fffffffffc9, 0x3fffffffffeb, 0x7fffffffff8d,
        0xffffffffffc5, 0x1ffffffffffaf, 0x3ffffffffffe5, 0x7ffffffffff7f,
        0xfffffffffffd1, 0x1fffffffffff91, 0x3fffffffffffdf, 0x7fffffffffffc9,
        0xfffffffffffffb, 0x1fffffffffffff3, 0x3ffffffffffffe5, 0x7ffffffffffffc9,
        0xfffffffffffffa3, 0x1fffffffffffffff, 0x3fffffffffffffc7,
        0x7fffffffffffffe7, 0xffffffffffffffc5 };


static void check_needed_bits(unsigned int nbits)
{
    if (nbits > 8 * sizeof(p_r_values_t))
        throw std::overflow_error("p_r_values_t is too small to store ideals, "
                "recompile with FLAGS_SIZE=\"-DSIZEOF_P_R_VALUES=8\"\n");
}

unsigned int renumber_t::needed_bits() const
{
    /* based on the lpbs and the number of sides, return 32 or 64 */
    uint64_t const p = previous_prime_of_powers_of_2[get_max_lpb()];
    uint64_t const vp = vp_from_p(p, get_nb_polys(), get_rational_side() >= 0, format);
    if (nbits (vp) <= 32)
        return 32;
    else
        return 64;
}

void renumber_t::set_format(int f)
{
    ASSERT_ALWAYS (above_all == above_bad);
    ASSERT_ALWAYS (above_cache == above_bad);
    ASSERT_ALWAYS(f == format_flat || f == format_binary);
    format = f;
}

void renumber_t::read_header(std::istream& is)
{
    std::ios_base::fmtflags const ff = is.flags();
    is >> std::dec;

    for(std::string s; std::ws(is).peek() == '#' ; getline(is, s) ) ;

    std::string s;
    getline(is, s);
    std::istringstream iss(s);
    int f = 0;
    if (iss >> f && (f == format_flat || f == format_binary)) {
        format = f;
    } else {
        throw cado::error(
                        "Renumber format error. Got {}, expected {} or {} instead. You must regenerate the renumber table with the freerel tool.",
                    f, format_flat, format_binary);
    }


    ASSERT_ALWAYS(above_all == above_add);
    {
        // we only have to parse the large prime bounds
        for(std::string s; std::ws(is).peek() == '#' ; getline(is, s) ) ;
        for(auto & x : lpb) is >> x;
        if (!is) throw parse_error("header");
        read_bad_ideals(is);
    }

    if (format == format_binary) {
        /* The header ends with one more line, which says where the data
         * begins. Consuming it to its end -- it is padded with spaces --
         * leaves us exactly on the data, which is the only way there
         * for a reader that cannot seek.
         */
        for(std::string t; std::ws(is).peek() == '#' ; getline(is, t) ) ;
        size_t off = 0;
        unsigned int es = 0;
        char bo = 0;
        is >> off >> es >> bo;
        if (!is)
            throw parse_error("header, binary data placement");
        if (es != 4 && es != 8)
            throw corrupted_table(fmt::format(
                        "binary renumber table with {}-byte entries", es));
        char const native =
            std::endian::native == std::endian::little ? 'L' : 'B';
        if (bo != native)
            throw corrupted_table(fmt::format(
                        "binary renumber table was written on a {}-endian "
                        "machine, this one is {}-endian",
                        bo == 'L' ? "little" : "big",
                        native == 'L' ? "little" : "big"));
        /* We write that offset as a multiple of 4096, so this can only
         * be a corrupt or foreign file. Reading it would mean a
         * misaligned pointer, which is not something to do quietly.
         */
        constexpr size_t a = alignof(decltype(flat_data)::value_type);
        if (off % a)
            throw corrupted_table(fmt::format(
                        "binary renumber table announces its data at offset"
                        " {}, which is not a multiple of {}", off, a));
        binary_data_offset = off;
        binary_element_size = es;

        std::string padding;
        getline(is, padding);
        if (!is)
            throw parse_error("header, binary data placement");
    }
    is.flags(ff);

    check_needed_bits(needed_bits());
}

/* This reads the bad ideals section of the new-format renumber file
 */
void renumber_t::read_bad_ideals(std::istream& is)
{
    std::ios_base::fmtflags const ff = is.flags();

    ASSERT_ALWAYS (above_all == above_bad);
    ASSERT_ALWAYS (above_cache == above_bad);
    above_bad = above_add;
    bad_ideals_max_p = 0;
    for(int side = 0; is && side < get_nb_polys(); side++) {
        for(std::string s; std::ws(is).peek() == '#' ; getline(is, s) ) ;
        int x;
        int n;
        /* The "side, number of bad ideals" is in decimal */
        is >> std::dec;
        is >> x >> n;
        ASSERT_ALWAYS(x == side);
        for( ; n-- ; ) {
            badideal b(is);
            p_r_side const x { (p_r_values_t) mpz_get_ui(b.p), (p_r_values_t) mpz_get_ui(b.r), side };
            bad_ideals.emplace_back(x, b);
            above_bad += b.nbad;
            if (x.p >= bad_ideals_max_p)
                bad_ideals_max_p = x.p;
        }
    }
    above_all = above_cache = above_bad;
    if (!is.good())
        throw parse_error("header, stream failed after reading bad ideals");
    is.flags(ff);
}

void renumber_t::write_header(std::ostream& os) const
{
    std::ios_base::fmtflags const ff = os.flags();
    os << std::dec;
    /* The traditional format doesn't even accept comments in the very
     * first line !
     */
    // the comment is never significant. However only the newer formats
    // give us the opportunity to stick in a marker for the file format.
    // So if we want (for the moment) to provide old-format renumber
    // files that can be parsed by old code, we can only convey the
    // format information as an (unparsed) side note.
    os << "# Renumber file using format " << format << "\n";

    /* Write the polynomials as comments */
    for (int i = 0; i < get_nb_polys() ; i++) {
        os << "# pol" << i << ": "
            << cxx_mpz_poly(cpoly[i]).print_poly("x")
            << "\n";
    }

    {
        os << format << "\n";
        // number of additional columns is implicit anyway.
        os << "# large prime bounds:\n";
        for (int i = 0; i < get_nb_polys() ; i++) {
            if (i) os << " ";
            os << lpb[i];
        }
        os << "\n";
    }

    os << "# " << above_add << " additional columns";
    if (has_merged_additional_column())
        os << " (combined for both sides)";
    os << "\n";
    os.flags(ff);
}

void renumber_t::write_bad_ideals(std::ostream& os) const
{
    std::ios_base::fmtflags const ff = os.flags();
    /* Write first the bad ideal information at the beginning of file */
    for(int side = 0; os && side < get_nb_polys(); side++) {
        os << "# bad ideals on side " << side << ": ";
        {
            unsigned int n = 0;
            for(auto const & b : bad_ideals) {
                if (b.first.side == side) {
                    if (n++) os << "+";
                    n++;
                    os << b.second.nbad;
                }
            }
            if (n == 0) os << "not used";
        }
        os << "\n";
        unsigned int n = 0;
        for(auto const & b : bad_ideals)
            if (b.first.side == side) n++;
        os << side << " " << n << "\n";
        for(auto const & b : bad_ideals) {
            if (b.first.side == side) {
                os << b.second;
            }
        }
    }
    os.flags(ff);
    os << "# renumber table for all indices above " << above_bad << ":\n";
}

/* The binary blob begins at a well-defined offset in the file, so that
 * it can be mmapped. Getting there means padding the header, and the
 * padding is written as one last comment line, so that a reader that
 * cannot seek (the file went through a compressor) reaches the blob by
 * reading lines until it has read that one.
 */
std::string renumber_t::header_string_with_padding() const
{
    ASSERT_ALWAYS(format == format_binary);

    /* A fixed number, not this machine's page size: the file must not
     * depend on where it was written. Readers align down to their own
     * page size anyway. All it really has to be is a multiple of the
     * alignment of an entry, which it amply is -- but a change to the
     * entry type should not get there unnoticed.
     */
    constexpr size_t blob_alignment = 4096;
    static_assert(
            blob_alignment % alignof(decltype(flat_data)::value_type) == 0,
            "the renumber table entries want an alignment that the "
            "binary format does not provide");

    std::ostringstream os;
    write_header(os);
    write_bad_ideals(os);
    std::string s = os.str();

    /* The last line of the header says where the data begins, and is
     * padded with spaces so that it ends exactly there. That is what
     * lets a reader that cannot seek land on the data: it reads that
     * line to its end, and it is there. The offset is written with a
     * fixed width, so the length of the line does not depend on it and
     * we can compute the one from the other.
     */
    auto placement = [this](size_t off) {
        return fmt::format("{:016d} {} {}", off, sizeof(p_r_values_t),
                std::endian::native == std::endian::little ? 'L' : 'B');
    };

    size_t const off = ((s.size() + placement(0).size() + 1
                + blob_alignment - 1) / blob_alignment) * blob_alignment;

    s += placement(off);
    s.append(off - s.size() - 1, ' ');
    s += '\n';
    ASSERT_ALWAYS(s.size() == off);

    return s;
}

std::vector<int> renumber_t::get_sides_of_additional_columns() const
{
    std::vector<int> res;
    for(int side = 0 ; side < get_nb_polys() ; side++) {
        mpz_poly_srcptr f = cpoly[side];
        /*
         * We used to not count an additional ideal J for degree 1
         * polynomials.  IMHO this is wrong.  We do have a J ideal,
         * and it must be counted so. It's cosmetic, but still.
         */
        if (/* f->deg > 1 && */ !mpz_poly_is_monic(f))
            res.push_back(side);
    }
    return res;
}

void renumber_t::use_additional_columns_for_dl()
{
    ASSERT_ALWAYS(above_all == 0);
    above_add = get_sides_of_additional_columns().size();
    if (has_merged_additional_column()) {
        /* This is a minor optimization. When we have two non-monic
         * sides, a single additional column is sufficient.
         */
        above_add = 1;
    }
    above_bad = above_add;
    above_cache = above_add;
    above_all = above_add;
}

void renumber_t::compute_bad_ideals()
{
    /* There are two use cases. Under normal circumstances, we're reading
     * the header, at this point. Which means that above_all =
     * above_cache = above_add = above_add. I.e. there are no bad ideals
     * stored, no cached ideals, and the table is empty.
     *
     * There's also the use case of
     * recompute_debug_number_theoretic_stuff(). There, we're recomputing
     * something that we mostly already know. In that case, we don't want
     * to destroy the old values of above_cache and above_all, of course.
     */
    // ASSERT_ALWAYS (above_all == above_bad);
    // ASSERT_ALWAYS (above_cache == above_bad);

    unsigned int const old_nbad = above_bad - above_add;

    /* most useful for traditional format, where we need to do this on
     * every open. Well, it's super cheap anyway
     *
     * Note that in all cases, we also use this function to compute bad
     * ideals when we create the table in the first place, from
     * freerel.cpp
     */
    above_bad = above_add;
    bad_ideals.clear();
    bad_ideals_max_p = 0;
    for(int side = 0 ; side < get_nb_polys() ; side++) {
        cxx_mpz_poly f(cpoly[side]);
        if (f->deg == 1) continue;
        for(auto const & b : badideals_for_polynomial(f, side)) {
            p_r_values_t const p = mpz_get_ui(b.p);
            p_r_values_t const r = mpz_get_ui(b.r);
            p_r_side const x { p, r, side };
            bad_ideals.emplace_back(x, b);
            above_bad += b.nbad;
            if (p >= bad_ideals_max_p)
                bad_ideals_max_p = p;
        }
    }
    if (old_nbad || above_all > above_add) {
        ASSERT_ALWAYS(above_bad - above_add == old_nbad);
    } else {
        above_all = above_cache = above_bad;
    }
}

void renumber_t::compute_ramified_primes()
{
    for(int side = 0 ; side < get_nb_polys() ; side++) {
        cxx_mpz disc;
        cxx_mpz_poly f(cpoly[side]);
        mpz_poly_discriminant(disc, f);
        mpz_mul(disc, disc, mpz_poly_lc(f));
        small_primes.push_back(trial_division(disc, 10000000, disc));
    }
}

void renumber_t::read_table(std::istream& is)
{
    stats_data_t stats;
    uint64_t nprimes = 0; // sigh... *must* be ulong for stats().
    stats_init(stats, stdout, &nprimes, 23, "Read", "primes", "", "p");
    for(std::string s; std::ws(is).peek() == '#' ; getline(is, s) ) ;
    for(p_r_values_t p, r ; is >> p >> r ; ) {
        flat_data.emplace_back(std::array<p_r_values_t, 2> {{ p, r }});
        above_all++;
        nprimes++;
    }
    stats_print_progress(stats, nprimes, 0, 0, 1);

    compute_index_from_p_cache();
}

void renumber_t::compute_index_from_p_cache()
{
    index_t i = 0;
    for( ; i < flat_data.size() ; i++)  {
        auto pvr = flat_data[i];
        p_r_values_t const p = pvr[0];
        if (p >> RENUMBER_MAX_LOG_CACHED)
            break;
        if (p < index_from_p_cache.size())
            continue;
        index_from_p_cache.insert(index_from_p_cache.end(),
                p - index_from_p_cache.size(),
                std::numeric_limits<index_t>::max());
        ASSERT_ALWAYS(index_from_p_cache.size() == p);
        index_from_p_cache.push_back(i);
    }
    above_cache = above_bad + i;
}

/* Read the binary blob. The fast path, which is the whole point of the
 * format, is the one where the entries have the width that this binary
 * was compiled for, and the file is a real file: then the table is
 * mmapped and nothing is read at all. Otherwise we fall back to reading
 * and converting, which is still vastly cheaper than parsing text.
 */
void renumber_t::read_table_binary(std::istream & is,
        std::string const & filename, bool may_mmap)
{
    using entry = decltype(flat_data)::value_type;

    size_t const row_bytes = 2 * size_t(binary_element_size);

    /* read_header() left us exactly on the data. If the stream can tell
     * us where that is, check that the header was telling the truth:
     * the mmap path trusts that offset, so it had better be right.
     */
    {
        auto const pos = is.tellg();
        if (pos >= 0 && size_t(pos) != binary_data_offset)
            throw corrupted_table(fmt::format(
                        "binary renumber table announces its data at offset"
                        " {}, but the header ends at offset {}",
                        binary_data_offset, size_t(pos)));
    }

    size_t nrows = 0;
    bool have_nrows = false;
    if (may_mmap) {
        struct stat sbuf[1];
        if (stat(filename.c_str(), sbuf) == 0) {
            size_t const total = size_t(sbuf->st_size);
            if (total < binary_data_offset
                    || (total - binary_data_offset) % row_bytes)
                throw corrupted_table("binary table has a truncated entry");
            nrows = (total - binary_data_offset) / row_bytes;
            have_nrows = true;
        }
    }

    /* read_header() has checked that binary_data_offset suits the
     * entries; mmapped_file::mapping takes care of rounding it down to
     * a page boundary and compensating.
     */
    if (nrows && have_nrows && binary_element_size == sizeof(p_r_values_t)) {
        using namespace mmap_allocator_details;
        mmapped_file source(filename, READ_ONLY,
                off_t(binary_data_offset), nrows * row_bytes);
        decltype(flat_data) y(mmap_allocator<entry>(
                    source, off_t(binary_data_offset), nrows));
        y.mmap(nrows);
        flat_data.swap(y);
        std::cout << fmt::format(
                "# INFO: {} entries mmapped from {} at offset {}\n",
                nrows, filename, binary_data_offset);
    } else {
        /* When we can't stat the file (it went through a compressor),
         * we don't know how many entries there are. The number of
         * ideals is close to the number of primes below the large prime
         * bounds, so this is a good enough guess to avoid reallocating.
         */
        if (have_nrows) {
            flat_data.reserve(nrows);
        } else {
            double guess = 0;
            for(auto l : lpb) {
                double const x = ldexp(1.0, int(l));
                guess += x / log(x);
            }
            flat_data.reserve(size_t(guess * 1.05) + 1024);
        }

        constexpr size_t batch = 1 << 16;
        std::vector<char> buf(batch * row_bytes);
        for(size_t done = 0 ; !have_nrows || done < nrows ; ) {
            size_t want = batch;
            if (have_nrows && nrows - done < want)
                want = nrows - done;
            is.read(buf.data(), std::streamsize(want * row_bytes));
            size_t const got = size_t(is.gcount()) / row_bytes;
            if (size_t(is.gcount()) % row_bytes)
                throw corrupted_table("binary table has a truncated entry");
            if (!got) {
                if (have_nrows)
                    throw corrupted_table("binary table is too short");
                break;
            }
            if (binary_element_size == sizeof(p_r_values_t)) {
                auto const * p = reinterpret_cast<entry const *>(buf.data());
                flat_data.append(p, p + got);
            } else if (binary_element_size == 4) {
                auto const * p = reinterpret_cast<uint32_t const *>(buf.data());
                for(size_t i = 0 ; i < got ; i++)
                    flat_data.push_back(entry {{
                            p_r_values_t(p[2*i]), p_r_values_t(p[2*i+1]) }});
            } else {
                auto const * p = reinterpret_cast<uint64_t const *>(buf.data());
                for(size_t i = 0 ; i < got ; i++) {
                    if (p[2*i] > std::numeric_limits<p_r_values_t>::max()
                            || p[2*i+1] > std::numeric_limits<p_r_values_t>::max())
                        throw prime_is_too_large(p[2*i]);
                    flat_data.push_back(entry {{
                            p_r_values_t(p[2*i]), p_r_values_t(p[2*i+1]) }});
                }
            }
            done += got;
        }
        std::cout << fmt::format(
                "# INFO: {} entries read from {} ({}-byte entries{})\n",
                flat_data.size(), filename, binary_element_size,
                binary_element_size == sizeof(p_r_values_t)
                ? "" : ", converted");
    }

    above_all = above_bad + flat_data.size();

    compute_index_from_p_cache();
}

void renumber_t::read_from_file(std::string const & filename, bool for_dl)
{
    if (cpoly.nsides() == 0)
        throw std::runtime_error("you must first call renumber_t::renumber_t(cxx_cado_poly const &)");
    ifstream_maybe_compressed is(filename);
    if (for_dl)
        use_additional_columns_for_dl();
    read_header(is);
    info(std::cout);
    if (format == format_binary)
        read_table_binary(is, filename, !is.is_pipe());
    else
        read_table(is);
    more_info(std::cout);
    
    /* It's used by inertia_from_p_r */
    compute_ramified_primes();
}

void renumber_t::recompute_debug_number_theoretic_stuff()
{
    /* explain_indexed_relations really insists on having the ramified
     * primes and the bad ideals computed anew, because that fills some
     * fields which are _not_ retrieved from the renumber output file
     */
    compute_ramified_primes();
    compute_bad_ideals();
}

std::string renumber_t::debug_data(index_t i) const
{
    p_r_side const x = p_r_from_index (i);

    std::string s = fmt::format("i=0x{:x} tab[i]=", i);

    if (is_additional_column (i)) {
        s += fmt::format("# added column for {}",
                has_merged_additional_column()
                ? "both sides combined"
                : fmt::format("side {}", x.side));
    } else if (is_bad(i)) {
        s += "# bad ideal";
        index_t j = i - above_add;
        for(auto const & b : bad_ideals) {
            if (j < (index_t) b.second.nbad) {
                s += fmt::format(" (number {}/{})", 1+j, b.second.nbad);
                break;
            }
            j -= b.second.nbad;
        }
        s += fmt::format(" above ({},{}) on side {}", x.p, x.r, x.side);
    } else {
        i -= above_bad;
        s += fmt::format(" (0x{:x},0x{:x})", flat_data[i][0], flat_data[i][1]);
        if (x.side == get_rational_side()) {
            s += fmt::format(" p=0x{:x} rat side {}", x.p, x.side);
        } else {
            s += fmt::format(" p=0x{:x} r=0x{:x} side {}", x.p, x.r, x.side);
            if (x.r == x.p)
                s += " proj";
        }
    }
    return s;
}

/* return valid sagemath code that describes the ideal */
std::string renumber_t::debug_data_sagemath(index_t i) const
{
    p_r_side x = p_r_from_index (i);
    if (is_additional_column (i)) {
        if (has_merged_additional_column()) {
            return "J0J1";
        } else {
            return fmt::format("J{0}", x.side);
        }
    } else if (is_bad(i)) {
        index_t j = i - above_add;
        for(auto const & b : bad_ideals) {
            if (j < (index_t) b.second.nbad) {
                if (b.second.sagemath_string.empty()) {
                    throw std::runtime_error("call compute_bad_ideals() first!\n");
                }
                return b.second.sagemath_string[j];
            }
            j -= b.second.nbad;
        }
    } else {
        i -= above_bad;
        if (x.side == get_rational_side()) {
            return fmt::format("OK{0}.fractional_ideal({1})", x.side, x.p);
        } else {
            /* XXX if p divides the discriminant, make sure that we treat
             * ramified primes correctly. Otherwise a simple and stupid
             * approach can work.
             */
            cxx_mpz_poly const f(cpoly[x.side]);
            for(auto const & y : small_primes[x.side]) {
                if (x.p == y.first)
                    return generic_sagemath_string(f, x.side, x.p, x.r);
            }

            /* back to the easy case */
            if (x.r == x.p) {
                return fmt::format("(OK{0}.fractional_ideal({1})+J{0})",
                        x.side, x.p);
            } else {
                return fmt::format("(OK{0}.fractional_ideal({1},alpha{0}-{2})*J{0})",
                        x.side, x.p, x.r);
            }
        }
    }
    throw corrupted_table("bad bad ideals");
}

/* return machine parseable code that describes the ideal */
/* The first token can be one of:
 *  J
 *  rat
 *  proj
 *  easy
 *  generic
 * which all lead to different ways of parsing the data.
 *
 * J is for additional columns. It is followed by the side info, or maybe
 * the sideS info in the odd case where we have a merged additional
 * column.
 *
 * rat is for rational primes, followed by the side index and the prime
 *
 * proj is for projective primes, followed by the side index and the prime
 *
 * easy is for easy algebrac primes, followed by the side index, the prime, and the root of alpha
 *
 * generic is for the rest. The following data is the side index, the
 * prime p, a denominator, and the coefficients of a polynomial which, when
 * divided by the denominator, yield an algebraic integer beta such that
 * (p,beta) is a prime ideal
 */
std::string renumber_t::debug_data_machine_description(index_t i) const
{
    p_r_side x = p_r_from_index (i);
    if (is_additional_column (i)) {
        if (has_merged_additional_column()) {
            return "J 0 1";
        } else {
            return fmt::format("J {0}", x.side);
        }
    } else if (is_bad(i)) {
        index_t j = i - above_add;
        for(auto const & b : bad_ideals) {
            if (j < (index_t) b.second.nbad) {
                if (b.second.machine_description.empty()) {
                    throw std::runtime_error("call compute_bad_ideals() first!\n");
                }
                return fmt::format("generic {} {}",
                        x.side,
                        b.second.machine_description[j]);
            }
            j -= b.second.nbad;
        }
    } else {
        // i -= above_bad;
        if (x.side == get_rational_side()) {
            return fmt::format("rat {} {}", x.side, x.p);
        } else {
            /* XXX if p divides the discriminant, make sure that we treat
             * ramified primes correctly. Otherwise a simple and stupid
             * approach can work.
             */
            cxx_mpz_poly const f(cpoly[x.side]);
            for(auto const & y : small_primes[x.side]) {
                if (x.p == y.first) {
                    return fmt::format("generic {} {}",
                            x.side,
                            generic_machine_description(f, x.side, x.p, x.r));
                }
            }

            /* back to the easy case */
            if (x.r == x.p) {
                return fmt::format("proj {} {}", x.side, x.p);
            } else {
                return fmt::format("easy {} {} {}", x.side, x.p, x.r);
            }
        }
    }
    throw corrupted_table("bad bad ideals");
}
/* can be called rather early, in fact (after the bad ideals computation) */
void renumber_t::info(std::ostream & os) const
{
    const char * P = "# INFO: ";
    os << "# Information on renumber table:\n";

    std::string const format_string =
        format == format_binary ? "binary" : "flat";

    os << P << "format = " << format_string << " (" << format << ")\n";
    os << P << "sizeof(p_r_values_t) = " << sizeof(p_r_values_t) << "\n";
    os << P << "nb_bits = " << needed_bits() << "\n";
    os << P << "number of polynomials = " << get_nb_polys() << "\n";
    if (get_rational_side() < 0)
        os << P << "There is no rational side\n";
    else
        os << P << "Polynomial on side " << get_rational_side() << " is rational\n";
    os << P << "#additional columns = " << above_add;
    if (above_add) {
        os << ", on side(s)";
        for(auto s : get_sides_of_additional_columns())
            os << " " << s;
        if (has_merged_additional_column()) {
            os << " (one column for both sides combined)";
        }
    }
    os << "\n";
    os << P << "#badideals = " << above_bad - above_add
        << ", above " << bad_ideals.size() << " (p,r) pairs"
        << " [max_p = " << bad_ideals_max_p << "]\n";
    for(int side = 0 ; side < get_nb_polys() ; side++) {
        os << P << "pol" << side << ": ";
        if (mpz_poly_is_monic(get_poly(side)))
            os << "monic\n";
        else
            os << "not monic\n";
    }
    for(int side = 0 ; side < get_nb_polys() ; side++) {
        os << P << "lpb" << side << " = " << lpb[side] << "\n";
    }
}
void renumber_t::more_info(std::ostream & os) const
{
    const char * P = "# INFO: ";
    os << "# Extra information on renumber table:\n";
    os << P << "size = " << size() << "\n";
}

void renumber_t::builder_declare_usage(cxx_param_list & pl)
{
    pl.declare_usage("renumber", "output file for renumbering table");
    pl.declare_usage("renumber_format", "format of the renumbering table (\"binary\" or \"flat\")");
    pl.declare_usage("renumber_rounds", "number of rounds used to build the renumbering table (the table is built by all threads at once, and one round's worth of fragments is held in memory)");
}

void renumber_t::builder_lookup_parameters(cxx_param_list & pl)
{
    pl.lookup("renumber");
    pl.lookup("renumber_format");
    pl.lookup("renumber_rounds");
}

/* This is the core of the renumber table building routine. Part of this
 * code used to exist in freerel.cpp file.
 */
struct renumber_t::builder{/*{{{*/
    renumber_t & R;
    std::ostream * os_p;
    renumber_t::hook * hook;
    unsigned int nrounds;
    stats_data_t stats;
    uint64_t nprimes = 0; // sigh... *must* be ulong for stats().
    index_t R_max_index; // we *MUST* follow it externally, since we're not storing the table in memory.

    builder(renumber_t & R, std::ostream * os_p, renumber_t::hook * hook,
            unsigned int nrounds)
        : R(R)
        , os_p(os_p)
        , hook(hook)
        , nrounds(nrounds)
        , R_max_index(R.get_max_index())
    {
        /* will print report at 2^10, 2^11, ... 2^23 computed primes
         * and every 2^23 primes after that */
        stats_init(stats, stdout, &nprimes, 23, "Processed", "primes", "", "p");
    }
    void progress() {
        if (stats_test_progress(stats))
            stats_print_progress(stats, nprimes, 0, 0, 0);
    }
    ~builder() {
        stats_print_progress(stats, nprimes, 0, 0, 1);
    }
    builder(builder const &) = delete;
    builder(builder &&) = delete;
    builder& operator=(builder const &) = delete;
    builder& operator=(builder &&) = delete;

    index_t operator()();
    void produce(fragment & F, unsigned long p0, unsigned long p1,
            gmp_randstate_ptr rstate) const;
    void finish(fragment & F) const;
    void emit(fragment & F);
};/*}}}*/

/* Compute the part of the table that lies above the primes in [p0, p1).
 * This is where all the time goes, and the only thing it touches is the
 * fragment it is given.
 */
void renumber_t::builder::produce(fragment & F,
        unsigned long p0, unsigned long p1,
        gmp_randstate_ptr rstate) const/*{{{*/
{
    int const nsides = R.get_nb_polys();

    F.clear();

    {
        /* One ideal per prime and per side, on average, for the sides
         * whose large prime bound is not exceeded (Chebotarev). Doing
         * this keeps the memory that the fragments hold close to the
         * size of the data they hold.
         */
        double ideals = 0;
        for(int side = 0 ; side < nsides ; side++) {
            double const b = ldexp(1.0, int(R.get_lpb(side)));
            if (double(p0) < b)
                ideals += nprimes_interval(double(p0), std::min(double(p1), b));
        }
        F.flat.reserve(size_t(ideals * 1.05) + 64);
        F.nroots.reserve(size_t(ideals * 1.05) + 64);
    }

    std::vector<std::vector<unsigned long>> all_roots(nsides);

    for(unsigned long const p : prime_range(p0, p1)) {
        F.nprimes_seen++;
        for (int side = 0; side < nsides; side++) {
            std::vector<unsigned long> & roots = all_roots[side];
            mpz_poly_srcptr f = R.get_poly(side);

            roots.clear();

            if (UNLIKELY(p >> R.get_lpb(side))) {
                continue;
            } else if (f->deg == 1) {
                roots.assign(1, 0);
            } else {
                /* Note that roots are sorted */
                roots = mpz_poly_roots(f, p, rstate);
            }

            /* Check for a projective root ; append it (so that the list of roots
             * is still sorted)
             */

            if ((int) roots.size() != R.get_poly_deg(side)
                    && mpz_divisible_ui_p(mpz_poly_lc(f), p))
                roots.push_back(p);

            /* take off bad ideals from the list, if any. */
            if (p <= R.bad_ideals_max_p) { /* can it be a bad ideal ? */
                for (size_t i = 0; i < roots.size() ; i++) {
                    unsigned long const r = roots[i];
                    if (!R.is_bad(p, r, side))
                        continue;
                    /* bad ideal -> remove this root from the list */
                    roots.erase(roots.begin() + i);
                    i--;
                }
            }
        }

        R.cook_into(p, all_roots, F);
    }
}/*}}}*/

/* Everything that needs to know where the fragment sits in the whole
 * table, but is still independent from one fragment to the next.
 */
void renumber_t::builder::finish(fragment & F) const/*{{{*/
{
    if (hook) {
        int const nsides = R.get_nb_polys();
        index_t idx = F.base;
        size_t cursor = 0;
        for(size_t g = 0 ; g < F.nroots.size() ; g += nsides) {
            /* cook_into() only records primes that have at least one
             * ideal above them, so this entry does exist, and its first
             * coordinate is the prime itself.
             */
            p_r_values_t const p = F.flat[cursor][0];
            (*hook)(R, p, idx, &F.nroots[g], F.hook_text);
            size_t n = 0;
            for(int side = 0 ; side < nsides ; side++)
                n += F.nroots[g + side];
            cursor += n;
            idx += n;
        }
        ASSERT_ALWAYS(cursor == F.flat.size());
    }

    if (os_p && R.get_format() != format_binary) {
        for(auto x : F.flat)
            F.text += fmt::format("{} {}\n", x[0], x[1]);
    }
}/*}}}*/

/* The only part that is done in sequence. */
void renumber_t::builder::emit(fragment & F)/*{{{*/
{
    if (os_p) {
        if (R.get_format() == format_binary) {
            using entry = decltype(F.flat)::value_type;
            os_p->write(
                    reinterpret_cast<char const *>(F.flat.data()),
                    std::streamsize(F.flat.size() * sizeof(entry)));
        } else {
            (*os_p) << F.text;
        }
    } else {
        R.flat_data.append(F.flat.begin(), F.flat.end());
    }

    if (hook)
        hook->flush(F.hook_text);

    nprimes += F.nprimes_seen;

    F.clear();

    progress();
}/*}}}*/

index_t renumber_t::builder::operator()()/*{{{*/
{
    /* Generate the renumbering table. The prime range is cut into
     * intervals that hold about the same number of primes, and these
     * are dealt with in rounds of one interval per thread. Rounds are
     * what bounds the memory that the fragments take: with the default
     * of 16 rounds, at most one sixteenth of the table is in flight.
     */
    unsigned long const lpbmax = 1UL << R.get_max_lpb();
    auto const nthreads = size_t(omp_get_max_threads());

    /* Intervals hold the same number of primes, but not the same
     * amount of work: finding the roots modulo a large prime costs
     * more than modulo a small one, and within one round the largest
     * prime is far above the smallest. Cutting each round in more
     * pieces than there are threads, and handing them out
     * dynamically, is what evens that out.
     */
    constexpr size_t granularity = 8;

    size_t per_round = nthreads * granularity;
    size_t nintervals = per_round * nrounds;
    {
        /* Cutting the range in pieces that are too small is pointless,
         * and the test suite goes as low as lpb=10.
         */
        double const np = double(lpbmax) / log(double(lpbmax));
        auto const most = size_t(std::max(1.0, np / 1024));
        if (nintervals > most) {
            nintervals = most;
            per_round = std::max(size_t(1), nintervals / nrounds);
        }
    }
    auto const splits = subdivide_primes_interval(2, lpbmax, nintervals);

    if (!os_p) {
        /* The number of ideals is close to the number of primes below
         * the large prime bounds, so this is a good enough guess to
         * never have to reallocate. Over-reserving costs address space
         * only.
         */
        double guess = 0;
        for(int side = 0 ; side < R.get_nb_polys() ; side++) {
            double const x = ldexp(1.0, int(R.get_lpb(side)));
            guess += x / log(x);
        }
        R.flat_data.reserve(size_t(guess * 1.05) + 1024);
    }

    std::vector<fragment> frags(per_round);
    std::vector<cxx_gmp_randstate> rstate_per_thread(nthreads);

    for(size_t i0 = 0 ; i0 < nintervals ; i0 += per_round) {
        size_t const n = std::min(per_round, nintervals - i0);

#pragma omp parallel for schedule(dynamic, 1)
        for(size_t j = 0 ; j < n ; j++)
            produce(frags[j], splits[i0 + j], splits[i0 + j + 1],
                    rstate_per_thread[omp_get_thread_num()]);

        for(size_t j = 0 ; j < n ; j++) {
            frags[j].base = R_max_index;
            R_max_index += frags[j].flat.size();
        }

#pragma omp parallel for schedule(dynamic, 1)
        for(size_t j = 0 ; j < n ; j++)
            finish(frags[j]);

        for(size_t j = 0 ; j < n ; j++)
            emit(frags[j]);

        if (!os_p) {
            R.above_all = R_max_index;
            ASSERT_ALWAYS(R.flat_data.size() == R.above_all - R.above_bad);
        }
    }

    if (!os_p)
        R.compute_index_from_p_cache();

    return R_max_index;
}/*}}}*/

index_t renumber_t::build(bool for_dl, hook * f)
{
    cxx_param_list pl;
    return build(pl, for_dl, f);
}

index_t renumber_t::build(cxx_param_list & pl, bool for_dl, hook * f)
{
    const char * renumberfilename = pl.lookup_old("renumber");
    const char * format_string = pl.lookup_old("renumber_format");

    unsigned int nrounds = 16;
    pl.parse("renumber_rounds", nrounds);
    if (!nrounds)
        throw std::runtime_error("renumber_rounds must be positive");

    if (format_string == nullptr) {
        set_format(format_binary);
    } else if (std::string(format_string) == "binary") {
        set_format(format_binary);
    } else if (std::string(format_string) == "flat") {
        set_format(format_flat);
    } else {
        throw std::runtime_error("cannot use this renumber format");
    }

    if (for_dl)
        use_additional_columns_for_dl();
    /* We can pass some hint primes as well, in a vector */
    compute_bad_ideals();

    info(std::cout);

    compute_ramified_primes();

    check_needed_bits(needed_bits());

    std::unique_ptr<std::ostream> out;

    if (renumberfilename) {
        out.reset(new ofstream_maybe_compressed(renumberfilename));

        if (format == format_binary) {
            *out << header_string_with_padding();
        } else {
            write_header(*out);
            write_bad_ideals(*out);
        }
    }

    index_t const ret = builder(*this, out.get(), f, nrounds)();

    more_info(std::cout);

    return ret;
}

renumber_t::const_iterator renumber_t::begin() const
{
    return const_iterator(*this, 0);
}
renumber_t::const_iterator renumber_t::end() const
{
    return const_iterator(*this, above_bad + flat_data.size());
}

renumber_t::p_r_side renumber_t::const_iterator::operator*() const {
    if (i < table->above_add) {
        /* See comment in p_r from index about the special case with 2
         * non-monic sides: we return {0,0,0} in that case, since
         * above_add == 1 */
        index_t j = 0;
        for(auto side : table->get_sides_of_additional_columns()) {
            if (j++ == i)
                return { 0, 0, side };
        }
    }
    if (i < table->above_bad) {
        /* annoying. we don't exactly have the pointer to the bad
         * ideal, we have to recover it. */
        index_t ii0 = i;
        for(auto const & I : table->bad_ideals) {
            if (ii0 < (index_t) I.second.nbad)
                return I.first;
            ii0 -= I.second.nbad;
        }
    }
    if (i == table->get_max_index()) {
        return p_r_side {
            std::numeric_limits<p_r_values_t>::max(),
            std::numeric_limits<p_r_values_t>::max(),
            0 };
    }
    auto const & [ p, vr ] = raw();
    return table->compute_p_r_side_from_p_vr(p, vr);
}

std::array<p_r_values_t, 2> renumber_t::const_iterator::raw() const
{
    ASSERT_ALWAYS(i >= table->above_bad);
    ASSERT_ALWAYS(i < table->get_max_index());
        return table->flat_data[i - table->above_bad];
}

renumber_t::const_iterator renumber_t::const_iterator::operator++(int)
{
    const_iterator ret = *this;
    ++*this;
    return ret;
}
renumber_t::const_iterator& renumber_t::const_iterator::operator++()
{
    /* We used to have quite fancy tests. It's now a lot simpler. The
     * only thing that we used to check, and that we could still choose
     * to check if we so wish (perhaps as an expensive assert) is the
     * fact that if an ideal was a bad one AND not the last bad ideal
     * with respect to the above_bad bound, then the resulting ideal
     * should still be bad.
     */
    i++;
    return *this;
}

int renumber_t::inertia_from_p_r(p_r_side x) const
{
    cxx_mpz_poly const f(cpoly[x.side]);
    if (f.degree() == 1) return 1;

    for(auto const & y : small_primes[x.side]) {
        if (x.p == y.first) {
            static std::mutex m;
            {
                std::lock_guard<std::mutex> const dummy(m);
                auto it = exceptional_inertia.find(x);
                if (it != exceptional_inertia.end())
                    return it->second;
            }
            cxx_gmp_randstate const state;
            int inertia = get_inertia_of_prime_ideal(f, x.p, x.r);
            {
                std::lock_guard<std::mutex> const dummy(m);
                exceptional_inertia[x] = inertia;
            }
            if (inertia != 1) {
                std::cerr << fmt::format("# computed exceptional inertia {} for ideal {},{} on side {}\n", inertia, x.p, x.r, x.side);
            }
            return inertia;
        }
    }
    return 1;
}
