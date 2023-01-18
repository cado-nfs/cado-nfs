#include "cado.h" // IWYU pragma: keep
// IWYU pragma: no_include <ext/alloc_traits.h>
#include <algorithm>
#include <fstream>      // std::ifstream // IWYU pragma: keep
#include <iomanip>      // std::hex // IWYU pragma: keep
#include <iostream>     // std::cout
#include <limits>
#include <list>
#include <memory>              // for allocator_traits<>::value_type, unique...
#include <sstream>      // std::ostringstream // IWYU pragma: keep
#include <stdexcept>
#include <string>
#include <vector>
#include <cstring>             // for strcmp
#include <cstdio> // stdout // IWYU pragma: keep
#include <climits> // UINT_MAX // IWYU pragma: keep
#include <gmp.h>               // for mpz_get_ui, mpz_divisible_ui_p, mpz_t
#include "badideals.hpp"
#include "cxx_mpz.hpp"         // for cxx_mpz
#include "mpz_mat.h"         // for mpz_mat
#include "numbertheory.hpp"  // for factorization_of_prime
#include "getprime.h"          // for getprime_mt, prime_info_clear, prime_i...
#include "gmp_aux.h"           // for ulong_isprime, nbits, mpz_get_uint64
#include "gzip.h"       // ifstream_maybe_compressed
#include "mod_ul.h"     // modulusul_t
#include "mpz_poly.h"   // mpz_poly
#include "omp_proxy.h" // IWYU pragma: keep
#include "params.h"            // for cxx_param_list, param_list_lookup_string
#include "renumber.hpp"
#include "rootfinder.h" // mpz_poly_roots
#include "stats.h"      // for the builder process
#include "macros.h"
#include "fmt/format.h"

#if defined(_GLIBCXX_DEBUG) && defined(_GLIBCXX_DEBUG_DISABLE_CHECK_PARTITIONED)
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
 *      filter/filter_galois.cpp
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

constexpr const int renumber_t::format_flat; /*{{{*/
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

#if 0
static renumber_t::corrupted_table cannot_find_p(p_r_values_t p) {/*{{{*/
    std::ostringstream os;
    os << std::hex;
    os << "cannot find data for prime 0x" << p;
    os << " ; note: isprime(p)==" << ulong_isprime(p);
#if SIZEOF_INDEX != 8
    os << "\n"
        << "Note: above 2^32 ideals or relations,"
        << " add FLAGS_SIZE=\"-DSIZEOF_P_R_VALUES=8 -DSIZEOF_INDEX=8\""
        << " to local.sh\n";
#endif
    return renumber_t::corrupted_table(os.str());
}/*}}}*/
#endif

static renumber_t::corrupted_table wrong_entry(p_r_values_t p, p_r_values_t vr) {/*{{{*/
    std::ostringstream os;
    os << std::hex;
    os << "above prime 0x" << p
        << ", the index 0x" << vr << " makes no sense";
    return renumber_t::corrupted_table(os.str());
}/*}}}*/
static renumber_t::corrupted_table prime_is_too_large(p_r_values_t p) {/*{{{*/
    std::ostringstream os;
    os << std::hex;
    os << "prime 0x" << p << " is too large!";
#if SIZEOF_INDEX != 8
    os << "\n"
        << "Note: above 2^32 ideals or relations,"
        << " add FLAGS_SIZE=\"-DSIZEOF_P_R_VALUES=8 -DSIZEOF_INDEX=8\""
        << " to local.sh\n";
#endif
    return renumber_t::corrupted_table(os.str());
}/*}}}*/
static renumber_t::corrupted_table prime_maps_to_garbage(int format, p_r_values_t p, index_t i, p_r_values_t q)/*{{{*/
{
    std::ostringstream os;
    os << std::hex;
    os << "cached index for prime p=0x" << p;
    os << " is 0x" << i << ", which points to q=0x" << q;
    if (format == renumber_t::format_flat)
        os << " (should be p)";
    else
        os << " (should be vp)";
    os << " ; note: isprime(p)==" << ulong_isprime(p);
    return renumber_t::corrupted_table(os.str());
}/*}}}*/
static renumber_t::corrupted_table prime_maps_to_garbage(int format, p_r_values_t p, index_t i)/*{{{*/
{
    std::ostringstream os;
    os << std::hex;
    os << "cached index for prime p=0x" << p;
    os << " is 0x" << i << ", which points nowhere";
    if (format == renumber_t::format_flat)
        os << " (should be p)";
    else
        os << " (should be vp)";
    os << " ; note: isprime(p)==" << ulong_isprime(p);
    return renumber_t::corrupted_table(os.str());
}/*}}}*/
static renumber_t::corrupted_table cannot_find_pr(renumber_t::p_r_side x)/*{{{*/
{
    std::ostringstream os;
    os << std::hex;
    os << "cannot find p=0x" << x.p << ", r=0x" << x.r << " on side " << x.side;
    /* #30012 and #30048 are cases where we see composites. #30048 in
     * particular sees the code crashing here, which is really cryptic.
     * We can alleviate that with a simple test.
     */
    if (!ulong_isprime(x.p)) {
        os << " ; NOTE that this p is not prime"
            << ", which is very unexpected. See bug #30048";
    }
    return renumber_t::corrupted_table(os.str());
}/*}}}*/
static renumber_t::corrupted_table cannot_lookup_p_a_b_in_bad_ideals(renumber_t::p_r_side x, int64_t a, uint64_t b)/*{{{*/
{
    std::ostringstream os;
    os << "failed bad ideal lookup for"
        << " (a,b)=(" << a << "," << b << ")"
        << " at p=0x" << std::hex << x.p
        << " on side " << x.side;
    return renumber_t::corrupted_table(os.str());
}/*}}}*/
static renumber_t::corrupted_table parse_error(std::string const & what)/*{{{*/
{
    std::ostringstream os;
    os << "parse error (" << what << ")";
    return renumber_t::corrupted_table(os.str());
}/*}}}*/
/* }}} */

/* {{{ helper functions relative to the "traditional" and "variant"
 * formats
 */

/* This one is a helper only, because we use it for both T ==
 * p_r_values_t and T == uint64_t */
template<typename T> inline T vp_from_p(T p, int n, int c, int format MAYBE_UNUSED)
{
    /* The final "+d" is not necessary, but we keep it for compatibility */
    int d = 0;
    return (n - c) * (p + 1) + d;
}

/* only used for format == format_{traditional,variant} */
p_r_values_t renumber_t::compute_vp_from_p (p_r_values_t p) const/*{{{*/
{
    int n = get_nb_polys();
    int c = (get_rational_side() >= 0);
    return vp_from_p(p, n, c, format);
}/*}}}*/


/* only used for format == format_{traditional,variant} */
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

/* only used for format == format_{traditional,variant} */
renumber_t::p_r_side renumber_t::compute_p_r_side_from_p_vr (p_r_values_t p, p_r_values_t vr) const/*{{{*/
{
    /* Note that vr is only used for the encoding of non-rational ideals.
     */
    p_r_side res { p, 0, 0 };

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

renumber_t::cooked renumber_t::cook(unsigned long p, std::vector<std::vector<unsigned long>> & roots) const
{
    cooked C;

    size_t total_nroots = 0;

    /* Note that all_roots always a root on the rational side, even
     * though it's only a zero -- the root itself isn't computed.
     */
    for (int i = 0; i < get_nb_polys() ; i++) {
        C.nroots.push_back(roots[i].size());
        total_nroots += roots[i].size();
    }

    if (total_nroots == 0) return C;

    for (int side = 0 ; side < get_nb_polys(); side++) {
        /* reverse the ordering of the *ROOTS* (not of the sides), because
         * our goal is to remain compatible with the old-format indexing
         */
        for (auto it = roots[side].rbegin() ; it != roots[side].rend() ; ++it) {
            p_r_side x { (p_r_values_t) p, (p_r_values_t) *it, side };
            C.flat.emplace_back(
                    std::array<p_r_values_t, 2> {{
                    (p_r_values_t) p,
                    compute_vr_from_p_r_side (x)
                    }});
        }
    }
    std::ostringstream os;
    for(auto x : C.flat)
        os << x[0] << " " << x[1] << "\n";
    C.text = os.str();
    return C;
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
        index_t n = I.second.nbad;
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
index_t renumber_t::get_first_index_from_p(p_r_values_t p) const
{
    // p_r_values_t side = x.side;
    if (p < index_from_p_cache.size()) {
        index_t i;
        /* Do this loop so that we can safely return something that is
         * valid even if p is not a prime.
         */
        for( ; ; p++) {
            if (p >= index_from_p_cache.size())
                return get_first_index_from_p(p);
            i = index_from_p_cache[p];
            if (i != std::numeric_limits<index_t>::max())
                break;
        }
        if (UNLIKELY(i >= flat_data.size()))
            throw prime_maps_to_garbage(format, p, i);
        if (UNLIKELY(flat_data[i][0] != p))
            throw prime_maps_to_garbage(format, p, i, flat_data[i][0]);
        return i;
    }

    /* Note that if lpbmax is > the width of the type, then
     * check_needed_bits() should have complained already.
     *
     * If lpbmax is == the width of the type, then the table fits (if
     * there's only one non-rational side). But we can't simply check p
     * >> lpbmax, since the check overflows...
     *
     */
    unsigned int lpbmax = *std::max_element(lpb.begin(), lpb.end());
    if (lpbmax < (8 * sizeof(p_r_values_t)) && UNLIKELY(p >> lpbmax))
        throw prime_is_too_large(p);

    std::array<p_r_values_t, 2> p0 {{ p, 0 }};
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
    index_t i = it - flat_data.begin();
    return i;
}

index_t renumber_t::index_from_p_r (p_r_side x) const
{
    index_t i;
    if (is_bad(i, x))
        return i;
    if (x.p == 0) {
        // additional columns. By convention they're attached to p==0.
        // Note that also by convention/tradition we use only a single
        // additional column where we have two non-monic sides, in which
        // case it's only counted on side 0.
        if (get_nb_polys() == 2 && get_sides_of_additional_columns().size() == 2) {
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
    i = get_first_index_from_p(x.p);
    p_r_values_t vr = compute_vr_from_p_r_side(x);

        /* The "flat" format has really simple lookups */
        for( ; flat_data[i][0] == x.p ; i++) {
            if (flat_data[i][1] == vr)
                return above_bad + i;
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
    if (p0 >> get_max_lpb())
        return get_max_index();
    index_t i;
    if (p0 < bad_ideals_max_p)
        i = 0;
    else
        i = above_bad + get_first_index_from_p(p0);
    for(const_iterator it(*this, i) ; it != end() && (*it).p < p0 ; ++it, ++i);
    return i;
}
index_t renumber_t::index_from_p(p_r_values_t p0, int side) const
{
    index_t i = index_from_p(p0);
    for(const_iterator it(*this, i) ; it != end() && (*it).side != side ; ++it, ++i);
    return i;
}
    
renumber_t::const_iterator renumber_t::iterator_from_p(p_r_values_t p0) const
{
    if (p0 >> get_max_lpb())
        return end();
    const_iterator it = begin();
    if (p0 >= bad_ideals_max_p)
        it.i = (above_bad + get_first_index_from_p(p0));
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
    std::vector<int> exps;
    for (auto const & I : bad_ideals) {
        if (x == I.first) {
            std::vector<int> exps;

            /* work here, and return */
            for(auto J : I.second.branches) {
                int k = J.k;
                p_r_values_t pk = x.p;
                ASSERT_ALWAYS(pk);
                for( ; --k ; ) {
                    p_r_values_t pk1 = pk * x.p;
                    /* we want to be sure that we don't wrap around ! */
                    ASSERT_ALWAYS(pk1 > pk);
                    pk = pk1;
                }
                p_r_values_t rk = mpz_get_uint64(J.r);
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
                    for(int v : J.v) {
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

/* additional columns _must_ be handled differently at this point */
renumber_t::p_r_side renumber_t::p_r_from_index (index_t i) const
{
    if (i < above_add) {
        /* In the special case where we have two non-monic polynomials,
         * we have a single additional column, hence above_add=1 and i=0.
         * We return {0,0,(aribitrarily the first of the two sides)} in
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

    p_r_values_t p = flat_data[i][0];
    p_r_values_t vr = flat_data[i][1];
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


void check_needed_bits(unsigned int nbits)
{
    if (nbits > 8 * sizeof(p_r_values_t))
        throw std::overflow_error("p_r_values_t is too small to store ideals, "
                "recompile with FLAGS_SIZE=\"-DSIZEOF_P_R_VALUES=8\"\n");
}

unsigned int renumber_t::needed_bits() const
{
    /* based on the lpbs and the number of sides, return 32 or 64 */
    uint64_t p = previous_prime_of_powers_of_2[get_max_lpb()];
    uint64_t vp = vp_from_p(p, get_nb_polys(), get_rational_side() >= 0, format);
    if (nbits (vp) <= 32)
        return 32;
    else
        return 64;
}

void renumber_t::set_format(int f)
{
    ASSERT_ALWAYS (above_all == above_bad);
    ASSERT_ALWAYS (above_cache == above_bad);
    ASSERT_ALWAYS(f == format_flat);
    format = f;
}

void renumber_t::read_header(std::istream& is)
{
    std::ios_base::fmtflags ff = is.flags();
    is >> std::dec;

    for(std::string s; std::ws(is).peek() == '#' ; getline(is, s) ) ;

    std::string s;
    getline(is, s);
    std::istringstream iss(s);
    int f = 0;
    if (iss >> f && (f == format_flat)) {
        format = f;
    } else {
        throw std::runtime_error(fmt::format(FMT_STRING(
                        "Renumber format error. Got {}, expected {} instead. You must regenerate the renumber table with the freerel tool."),
                    f, format_flat));
    }

    ASSERT_ALWAYS(above_all == above_add);
    {
        // we only have to parse the large prime bounds
        for(std::string s; std::ws(is).peek() == '#' ; getline(is, s) ) ;
        for(auto & x : lpb) is >> x;
        if (!is) throw parse_error("header");
        read_bad_ideals(is);
    }
    is.flags(ff);

    check_needed_bits(needed_bits());
}

/* This reads the bad ideals section of the new-format renumber file
 */
void renumber_t::read_bad_ideals(std::istream& is)
{
    std::ios_base::fmtflags ff = is.flags();

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
            p_r_side x { (p_r_values_t) mpz_get_ui(b.p), (p_r_values_t) mpz_get_ui(b.r), side };
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
    std::ios_base::fmtflags ff = os.flags();
    os << std::dec;
    /* The traditional format doesn't even accept comments in the very
     * first line !
     */
    // the comment is never significant. However only the newer formats
    // give us the opportunity to stick in a marker for the file format.
    // So if we want (for the moment) to provide old-format renumber
    // files that can be parsed by old code, we can only convey the
    // format information as an (unparsed) side note.
    os << "# Renumber file using format " << format << std::endl;

    /* Write the polynomials as comments */
    for (int i = 0; i < get_nb_polys() ; i++) {
        os << "# pol" << i << ": "
            << cxx_mpz_poly(cpoly->pols[i]).print_poly("x")
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
    if (get_nb_polys() == 2 && get_sides_of_additional_columns().size() == 2)
        os << " (combined for both sides)";
    os << "\n";
    os.flags(ff);
}

void renumber_t::write_bad_ideals(std::ostream& os) const
{
    std::ios_base::fmtflags ff = os.flags();
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
        os << std::endl;
        unsigned int n = 0;
        for(auto const & b : bad_ideals)
            if (b.first.side == side) n++;
        os << side << " " << n << std::endl;
        for(auto const & b : bad_ideals) {
            if (b.first.side == side) {
                os << b.second;
            }
        }
    }
    os.flags(ff);
    os << "# renumber table for all indices above " << above_bad << ":\n";
}

std::vector<int> renumber_t::get_sides_of_additional_columns() const
{
    std::vector<int> res;
    for(int side = 0 ; side < get_nb_polys() ; side++) {
        mpz_poly_srcptr f = cpoly->pols[side];
        if (f->deg > 1 && !mpz_poly_is_monic(f))
            res.push_back(side);
    }
    return res;
}

void renumber_t::use_additional_columns_for_dl()
{
    ASSERT_ALWAYS(above_all == 0);
    above_add = get_sides_of_additional_columns().size();
    if (get_nb_polys() == 2 && get_sides_of_additional_columns().size() == 2) {
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

    unsigned int old_nbad = above_bad - above_add;

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
        cxx_mpz_poly f(cpoly->pols[side]);
        if (f->deg == 1) continue;
        for(auto const & b : badideals_for_polynomial(f, side)) {
            p_r_values_t p = mpz_get_ui(b.p);
            p_r_values_t r = mpz_get_ui(b.r);
            p_r_side x { p, r, side };
            bad_ideals.emplace_back(x, b);
            above_bad += b.nbad;
            if (p >= bad_ideals_max_p)
                bad_ideals_max_p = p;
        }
    }
    if (old_nbad) {
        ASSERT_ALWAYS(above_bad - above_add == old_nbad);
    } else {
        above_all = above_cache = above_bad;
    }
}

void renumber_t::compute_ramified_primes()
{
    for(int side = 0 ; side < get_nb_polys() ; side++) {
        cxx_mpz disc;
        cxx_mpz_poly f(cpoly->pols[side]);
        mpz_poly_discriminant(disc, f);
        mpz_mul(disc, disc, f->coeff[f->deg]);
        small_primes.push_back(trial_division(disc, 10000000, disc));
    }
}

void renumber_t::use_cooked(p_r_values_t p, cooked const & C)
{
    if (C.empty()) return;
    /* In the current format, we have
     * above_all - above_bad == flat_daat.size().
     * Note that this used to not be the case with the old formats.
     */
    index_t pos_hard = flat_data.size();
    above_all = use_cooked_nostore(above_all, p, C);
    flat_data.insert(flat_data.end(), C.flat.begin(), C.flat.end());
    if (!(p >> RENUMBER_MAX_LOG_CACHED) && p >= index_from_p_cache.size()) {
        index_from_p_cache.insert(index_from_p_cache.end(),
                p - index_from_p_cache.size(),
                std::numeric_limits<index_t>::max());
        ASSERT_ALWAYS(index_from_p_cache.size() == p);
        index_from_p_cache.push_back(pos_hard);
        above_cache = above_all;
    }
}
index_t renumber_t::use_cooked_nostore(index_t n0, p_r_values_t p MAYBE_UNUSED, cooked const & C)
{
    if (C.empty()) return n0;
    for(auto n : C.nroots) n0 += n;
    return n0;
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

    {
        index_t i = 0;
        for( ; i < flat_data.size() ; i++)  {
            auto pvr = flat_data[i];
            p_r_values_t p = pvr[0];
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
}

void renumber_t::read_from_file(const char * filename, int for_dl)
{
    ifstream_maybe_compressed is(filename);
    if (for_dl)
        use_additional_columns_for_dl();
    read_header(is);
    info(std::cout);
    read_table(is);
    more_info(std::cout);
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
    p_r_side x = p_r_from_index (i);
    std::ostringstream os;

    os << "i=0x" << std::hex << i;

    if (is_additional_column (i)) {
        if (get_nb_polys() == 2 && get_sides_of_additional_columns().size() == 2) {
            os << " tab[i]=#"
                << " added column for both sides combined";
        } else {
            os << " tab[i]=#"
                << " added column for side " << std::dec << x.side;
        }
    } else if (is_bad(i)) {
        os << " tab[i]=#"
            << " bad ideal";
        index_t j = i - above_add;
        for(auto const & b : bad_ideals) {
            if (j < (index_t) b.second.nbad) {
                os << std::dec
                    << " (number " << (1+j) << "/" << b.second.nbad << ")";
                break;
            }
            j -= b.second.nbad;
        }
        os  << " above"
            << " (" << x.p << "," << x.r << ")"
            << " on side " << x.side;
    } else {
        i -= above_bad;
        os << std::hex;
        {
            os << " tab[i]=";
            os << " (0x" << flat_data[i][0] << ",0x" << flat_data[i][1] << ")";
        }
        os << " p=0x" << x.p;
        if (x.side == get_rational_side()) {
            os << " rat";
            os << " side " << std::dec << x.side;
        } else {
            os << " r=0x" << x.r;
            os << " side " << std::dec << x.side;
            if (x.r == x.p)
                os << " proj";
        }
    }

    return os.str();
}

/* return valid sagemath code that describes the ideal */
std::string renumber_t::debug_data_sagemath(index_t i) const
{
    p_r_side x = p_r_from_index (i);
    if (is_additional_column (i)) {
        if (get_nb_polys() == 2 && get_sides_of_additional_columns().size() == 2) {
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
            return fmt::format("OK{0}.ideal({1})", x.side, x.p);
        } else {
            /* XXX if p divides the discriminant, make sure that we treat
             * ramified primes correctly. Otherwise a simple and stupid
             * approach can work.
             */
            cxx_mpz_poly f(cpoly->pols[x.side]);
            cxx_gmp_randstate state;
            for(auto const & y : small_primes[x.side]) {
                if (x.p == y.first)
                    return generic_sagemath_string(f, x.side, x.p, x.r, state);
            }

            /* back to the easy case */
            if (x.r == x.p) {
                return fmt::format("(OK{0}.ideal({1})+J{0})",
                        x.side, x.p);
            } else {
                return fmt::format("(OK{0}.fractional_ideal({1},alpha{0}-{2})*J{0})",
                        x.side, x.p, x.r);
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

    std::string format_string = "flat";

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
        os << ", on sides";
        for(auto s : get_sides_of_additional_columns())
            os << " " << s;
        if (get_nb_polys() == 2 && get_sides_of_additional_columns().size() == 2) {
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
    os << P << "size = " << get_size() << "\n";
}

void renumber_t::builder_declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "renumber", "output file for renumbering table");
    param_list_decl_usage(pl, "renumber_format", "format of the renumbering table (\"flat\")");
}

void renumber_t::builder_lookup_parameters(cxx_param_list & pl)
{
    param_list_lookup_string(pl, "renumber");
    param_list_lookup_string(pl, "renumber_format");
}

/* This is the core of the renumber table building routine. Part of this
 * code used to exist in freerel.cpp file.
 */
struct renumber_t::builder{/*{{{*/
    struct prime_chunk {/*{{{*/
        bool preprocess_done = false;
        std::vector<unsigned long> primes;
        std::vector<renumber_t::cooked> C;
        prime_chunk(std::vector<unsigned long> && primes) : primes(primes) {}
        private:
        prime_chunk() = default;
    };/*}}}*/

    renumber_t & R;
    std::ostream * os_p;
    renumber_t::hook * hook;
    stats_data_t stats;
    uint64_t nprimes = 0; // sigh... *must* be ulong for stats().
    index_t R_max_index; // we *MUST* follow it externally, since we're not storing the table in memory.
    builder(renumber_t & R, std::ostream * os_p, renumber_t::hook * hook)
        : R(R)
        , os_p(os_p)
        , hook(hook)
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
    index_t operator()();
    void preprocess(prime_chunk & P, gmp_randstate_ptr rstate);
    void postprocess(prime_chunk & P);
};/*}}}*/

void renumber_t::builder::preprocess(prime_chunk & P, gmp_randstate_ptr rstate)/*{{{*/
{
    ASSERT_ALWAYS(!P.preprocess_done);
    /* change x (list of input primes) into the list of integers that go
     * to the renumber table, and then set "done" to true.
     * This is done asynchronously.
     */
    for(auto p : P.primes) {
        std::vector<std::vector<unsigned long>> all_roots;
        for (int side = 0; side < R.get_nb_polys(); side++) {
            std::vector<unsigned long> roots;
            mpz_poly_srcptr f = R.get_poly(side);

            if (UNLIKELY(p >> R.get_lpb(side))) {
                all_roots.emplace_back(roots);
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
                    && mpz_divisible_ui_p(f->coeff[f->deg], p))
                roots.push_back(p);

            /* take off bad ideals from the list, if any. */
            if (p <= R.bad_ideals_max_p) { /* can it be a bad ideal ? */
                for (size_t i = 0; i < roots.size() ; i++) {
                    unsigned long r = roots[i];
                    if (!R.is_bad(p, r, side))
                        continue;
                    /* bad ideal -> remove this root from the list */
                    roots.erase(roots.begin() + i);
                    i--;
                }
            }
            all_roots.emplace_back(roots);
        }

        /* Data is written in the temp buffer in a way that is not quite
         * similar to the renumber table, but still close enough.
         */
        P.C.emplace_back(R.cook(p, all_roots));
    }
#pragma omp atomic write
    P.preprocess_done = true;
}/*}}}*/

void renumber_t::builder::postprocess(prime_chunk & P)/*{{{*/
{
    ASSERT_ALWAYS(P.preprocess_done);
    /* put all entries from x into the renumber table, and also print
     * to freerel_file any free relation encountered. This is done
     * synchronously.
     *
     * (if freerel_file is NULL, store only into the renumber table)
     */
    for(size_t i = 0; i < P.primes.size() ; i++) {
        p_r_values_t p = P.primes[i];
        renumber_t::cooked & C = P.C[i];

        if (hook) (*hook)(R, p, R_max_index, C);

        if (os_p) {
            R_max_index = R.use_cooked_nostore(R_max_index, p, C);
            (*os_p) << C.text;
        } else {
            ASSERT_ALWAYS(R_max_index == R.get_max_index());
            R.use_cooked(p, C);
            R_max_index = R.get_max_index();
        }

        nprimes++;
    }
    /* free memory ! */
    P.primes.clear();
    P.C.clear();
    progress();
}/*}}}*/

index_t renumber_t::builder::operator()()/*{{{*/
{
    /* Generate the renumbering table. */

    constexpr const unsigned int granularity = 1024;

    std::vector<cxx_gmp_randstate> rstate_per_thread(omp_get_max_threads());
#pragma omp parallel default(none) shared(rstate_per_thread)
    {
#pragma omp single
        {
            prime_info pi;
            prime_info_init(pi);
            std::list<prime_chunk> inflight;
            unsigned long lpbmax = 1UL << R.get_max_lpb();
            unsigned long p = 2;
            for (; p <= lpbmax || !inflight.empty() ;) {
                if (p <= lpbmax) {
                    std::vector<unsigned long> pp;
                    pp.reserve(granularity);
                    for (; p <= lpbmax && pp.size() < granularity;) {
                        pp.push_back(p);
                        p = getprime_mt(pi); /* get next prime */
                    }
                    inflight.emplace_back(std::move(pp));
                    /* do not use a c++ reference for the omp
                     * firstprivate construct. It does not do what we
                     * want. (I saw a _copy_ !)
                     */
                    prime_chunk * latest(&inflight.back());
#pragma omp task firstprivate(latest) default(none) shared(rstate_per_thread)
                    {
                        preprocess(*latest, rstate_per_thread[omp_get_thread_num()]);
                    }
                } else {
#pragma omp taskwait
                }

                for ( ; !inflight.empty() ; ) {
                    bool ready;
                    prime_chunk & next(inflight.front());
#pragma omp atomic read
                    ready = next.preprocess_done;
                    if (!ready)
                        break;
                    postprocess(next);
                    inflight.pop_front();
                }
            }
            prime_info_clear(pi);
        }
    }

    return R_max_index;
}/*}}}*/

index_t renumber_t::build(int for_dl, hook * f)
{
    cxx_param_list pl;
    return build(pl, for_dl, f);
}

index_t renumber_t::build(cxx_param_list & pl, int for_dl, hook * f)
{
    const char * renumberfilename = param_list_lookup_string(pl, "renumber");
    const char * format_string = param_list_lookup_string(pl, "renumber_format");

    if (format_string == NULL) {
        set_format(format_flat);
    } else if (strcmp(format_string, "flat") == 0) {
        format = format_flat;
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

        write_header(*out);
        write_bad_ideals(*out);
    }

    index_t ret = builder(*this, out.get(), f)();

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
    if (i < table.above_add) {
        /* See comment in p_r from index about the special case with 2
         * non-monic sides */
        index_t j = 0;
        for(auto side : table.get_sides_of_additional_columns()) {
            if (j++ == i)
                return { 0, 0, side };
        }
    }
    if (i < table.above_bad) {
        /* annoying. we don't exactly have the pointer to the bad
         * ideal, we have to recover it. */
        index_t ii0 = i;
        for(auto const & I : table.bad_ideals) {
            if (ii0 < (index_t) I.second.nbad)
                return I.first;
            ii0 -= I.second.nbad;
        }
    }
    if (i == table.get_max_index()) {
        return p_r_side {
            std::numeric_limits<p_r_values_t>::max(),
            std::numeric_limits<p_r_values_t>::max(),
            0 };
    }
    {
        p_r_values_t p = table.flat_data[i-table.above_bad][0];
        p_r_values_t vr = table.flat_data[i-table.above_bad][1];
        return table.compute_p_r_side_from_p_vr(p, vr);
    }
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
