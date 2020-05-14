#include "cado.h"
#include <sstream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <exception>
#include <iomanip>
#include <algorithm>
#include <map>
#include "gzip.h"       // ifstream_maybe_compressed
#include "renumber.hpp"
#include "badideals.hpp"
#include "utils.h"

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
filter/dup2.cpp:        pr[i].h = renumber_get_index_from_p_r(renumber_tab, pr[i].p, r,
sieve/las-dlog-base.cpp:        index_t h = renumber_get_index_from_p_r(tab, p, r, side);
sieve/fake_rels.cpp:                index = renumber_get_index_from_p_r(ren_tab, p, r, side);
sieve/fake_rels.cpp:      index_t indq = renumber_get_index_from_p_r(ren_table, q, roots[0], sqside);
sieve/fake_rels.cpp:      indq = renumber_get_index_from_p_r(ren_table, q, roots[0], sqside);

filter/filter_galois.cpp:      renumber_get_p_r_from_index(tab, &p, &rr, &side, i, cpoly);// XXX forward iter
filter/reconstructlog.cpp:        renumber_get_p_r_from_index (tab, &p, &r, &side, i, poly);// XXX forward iter
sieve/fake_rels.cpp:        int side = renumber_get_side_from_index(ren_tab, i, cpoly);// XXX forward iter
sieve/fake_rels.cpp:                renumber_get_p_r_from_index(ren_tab, &p, &r, &side, i, cpoly);// XXX forward iter
filter/reconstructlog.cpp:          int side = renumber_get_side_from_index (data.renum_tab, h, data.poly);//XXX MNFS only
misc/convert_rels.cpp:        renumber_get_p_r_from_index(renumber_table, &p, &r, &side, i, poly);

filter/dup2.cpp:          && renumber_is_bad (&nb, &first_index, renumber_tab, pr[i].p, r,
(followed by handle_bad_ideals which returns a vector of valuations -- the (a,b) pair is required. The norm valuation is needed as input as well, I believe.)
sieve/las-dlog-base.cpp:        if (renumber_is_bad (NULL, NULL, tab, p, r, side))
(which acknowledges that the very bit of code that is in dup2 is missing here, and we would really like to have it)
sieve/fake_rels.cpp:            if (renumber_is_bad (&nb, &index, ren_tab, p, r, side)) {
(there, we do something very simple, and only pick a random index for a bad ideal, in [index, index+nb) )

sieve/fake_rels.cpp:                renumber_badideal_get_p_r_below(ren_tab, &p, &r, &side, i);

filter/reconstructlog.cpp:        if (renumber_is_additional_column (tab, i))
sieve/fake_rels.cpp:        if (renumber_is_additional_column(ren_tab, i)) {
 * }}} */

/* In the above uses, the lookup (side,p,r) -> index is critical in the
 * dup2 and descent steps, while the lookup index -> (side,p,r) is most
 * often used only for forward iterators, where the context can be used
 * to carry some extra information.
 */

/* There are various strategies that can be used to have an in-memory
 * structure that makes these lookups possible.
 */

constexpr const int renumber_format_traditional = 20130603; /*{{{*/
/* The traditional way that cado-nfs has been using for quite a while now
 * (2013-2020) uses the following memory layout.
 * In memory, the table is stored as a sequence of integers.  The info
 * for each prime is stored as a sequence (vp,vr0,vr1,...,vrk) with
 * vp>=vr0>...>vrk. vp is a representative for p.
 * Whenever two adjacent integers in the in-memory table are in strictly
 * increasing order, the latter is a marker that indicates a new p.
 *
 * In the traditional scheme, vr for an ideal (side,p,r) is
 * side'*(p+1)+r, where side' is the side index counted among the
 * non-rational sides. vp is (N-c)*(p+1)-1+c, with N the number of sides
 * and c is 1 if there is a rational side. (In fact, the last "+c" is
 * useless.)
 *
 * The goal of this traditional scheme is that exactly one integer
 * appears in the database for each ideal -- namely, p is always
 * implicit. This means, in particular, that when we have no rational
 * side, there is no room to store _all_ roots: one of them must be
 * sacrificed. The big advantage of this layout is that the
 * index_from_p_r lookup, after the initial dichotomy, is
 * straightforward. The big disadvantage is that p_r_from_index has to
 * recompute the roots mod p in order to return the root mod p that was
 * overwritten in order to store p. (This overhead occurs even in the
 * case of a forward iterator.)
 */
/*}}}*/

constexpr const int renumber_format_variant = 20199999; /*{{{*/
/* Variations around the previous scheme may be tried: we may try to
 * store _all_ roots on non-rational sides.
 * (The rational root is still implicit, if there is a rational side).
 * Basically, we define vp and vr as above, but we include all vr's. This
 * avoids the overhead of recomputing the roots when doing the
 * p_r_from_index lookup.
 * There are two disadvantages.
 *  - First, storage is slightly larger. Instead of one integer per
 *  ideal, we store approximately 1.4 integers. The density of ideals
 *  in either side is given by D, and the number of integers that we
 *  store is given by F, with the following magma code -- assuming
 *  generic Galois groups on both sides:
      D:=func<d1,d2|&+[m+n eq 0 select 0 else m+n where m is #Fix(s) where n is #Fix(t): s in S, t in T] * 1.0/#S/#T where S is SymmetricGroup(d1) where T is SymmetricGroup(d2)>;
      F:=func<d1,d2|&+[m+n eq 0 select 0 else 1+m+n where m is #Fix(s) where n is #Fix(t): s in S, t in T] * 1.0/#S/#T where S is SymmetricGroup(d1) where T is SymmetricGroup(d2)>;
      F(3,4)/D(3,4) is 1.4375 (on the first thought)
 *  - But the real problem is not storage. The index_from_p_r lookup is
 *  destroyed by this layout, since by making p explicit, we have data in
 *  the table that does not correspond to an ideal. This means that we
 *  need to store _more_ stuff to find the consecutive index, once we've
 *  found the position of vr in the table. We may, right after vp, store
 *  vp+I, with I the consecutive index of the first ideal above p. This
 *  means that we have in the table increasing sequences of length 2 at
 *  each p. This keeps the synchronization property. Note that the
 *  p_r_from_index lookup is now _also_ impacted, since we moved from
 *  O(1) to a mildly logarithmic cost (amortized O(1) for forward
 *  iterators). This extra info costs somewhat more. By modifying F
 *  above, we reach:
      F(3,4)/D(3,4) is 1.875 
 */
/*}}}*/

constexpr const int renumber_format_flat = 20200515; /*{{{*/
/* My personal preference would be to bite the bullet and allow for two
 * integers per ideal, always storing both p and r explicitly:
 *  - The _from_index lookups would be trivial (we could keep the encoding
 *    of r as vr = side'*(p+1)+r).
 *  - The _from_p_r lookups would also be doable with very simple-minded
 *    dichotomy (even the standard library can be put to some use here)
 *  - The database in itself would be human readable.
 */
/*}}}*/

constexpr const int renumber_format = renumber_format_traditional;

/* {{{ exceptions */
renumber_t::corrupted_table::corrupted_table(std::string const & x)
    : std::runtime_error(std::string("Renumber table is corrupt: ") + x)
{}
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
static renumber_t::corrupted_table cannot_find_i(index_t i) {/*{{{*/
    std::ostringstream os;
    os << std::hex;
    os << "cannot find data with index 0x" << i;
    return renumber_t::corrupted_table(os.str());
}/*}}}*/
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
static renumber_t::corrupted_table prime_maps_to_garbage(p_r_values_t p, index_t i, p_r_values_t q)/*{{{*/
{
    std::ostringstream os;
    os << std::hex;
    os << "cached index for prime p=0x" << p;
    os << " is " << i << ", which points to q=0x" << q;
    if (renumber_format == renumber_format_flat)
        os << " (should be p)";
    else
        os << " (should be vp)";
    os << " ; note: isprime(p)==" << ulong_isprime(p);
    return renumber_t::corrupted_table(os.str());
}/*}}}*/
static renumber_t::corrupted_table cannot_find_pr(renumber_t::p_r_side x, p_r_values_t vp, p_r_values_t vr)/*{{{*/
{
    std::ostringstream os;
    os << std::hex;
    os << "cannot find p=0x" << x.p << ", r=0x" << x.r << " on side " << x.side;
    os << "; note: vp=0x" << vp << ", vr=0x" << vr;
    return renumber_t::corrupted_table(os.str());
}/*}}}*/
static renumber_t::corrupted_table cannot_find_pr(renumber_t::p_r_side x)/*{{{*/
{
    std::ostringstream os;
    os << std::hex;
    os << "cannot find p=0x" << x.p << ", r=0x" << x.r << " on side " << x.side;
    return renumber_t::corrupted_table(os.str());
}/*}}}*/
static renumber_t::corrupted_table cannot_lookup_p_a_b_in_bad_ideals(renumber_t::p_r_side x, int64_t a, uint64_t b)
{
    std::ostringstream os;
    os << "failed bad ideal lookup for"
        << " (a,b)=(" << a << "," << b << ")"
        << " at p=0x" << std::hex << x.p
        << " on side " << x.side;
    return renumber_t::corrupted_table(os.str());
}
/* }}} */

/* {{{ helper functions relative to the "traditional" and "variant"
 * formats
 */

/* This one is a helper only, because we use it for both T ==
 * p_r_values_t and T == uint64_t */
template<typename T> inline T vp_from_p(T p, int n, int c)
{
    /* The final "+c" is not necessary, but we keep it for compatibility */
    return (n - c) * (p + 1) - 1 + c;
}

/* only used for renumber_format == renumber_format_{traditional,variant} */
p_r_values_t renumber_t::compute_vp_from_p (p_r_values_t p) const/*{{{*/
{
    int n = get_nb_polys();
    int c = (get_rational_side() >= 0);
    return vp_from_p(p, n, c);
}/*}}}*/

/* Inverse function of compute_vp_from_p */
/* only used for renumber_format == renumber_format_{traditional,variant} */
p_r_values_t renumber_t::compute_p_from_vp (p_r_values_t vp) const/*{{{*/
{
    int n = get_nb_polys();
    int c = (get_rational_side() >= 0);
    return (vp + 1 - c) / (n - c);
}/*}}}*/

/* only used for renumber_format == renumber_format_{traditional,variant} */
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

/* only used for renumber_format == renumber_format_{traditional,variant} */
renumber_t::p_r_side renumber_t::compute_p_r_side_from_p_vr (p_r_values_t p, p_r_values_t vr) const/*{{{*/
{
    /* Note that vr is only used for the encoding of non-rational ideals.
     */
    p_r_side res { p, 0, 0 };

    res.r = vr;
    for(res.side = 0 ; res.side < (int) get_nb_polys() ; res.side++) {
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
        } else if (renumber_format == renumber_format_traditional && traditional_get_largest_nonbad_root_mod_p(res)) {
            /* that is the _really, really_ annoying thing with the
             * traditional format */
            return res;
        }
    }
    throw wrong_entry(p, vr);
}/*}}}*/
/* }}} */

/* sort in decreasing order. Faster than qsort for ~ < 15 values in r[] */
/* only used for renumber_format == renumber_format_{traditional,variant} */
/* XXX This is total legacy, and should go away soon (we only temporarily
 * keep it for measurement). See test-sort in
 * tests/utils. This code is actually slow in most cases, and clearly
 * outperformed by the code in iqsort.h
 */
inline void renumber_sort_ul (unsigned long *r, size_t n)
{
    unsigned long rmin;

    if (UNLIKELY (n < 2))
        return;

    if (UNLIKELY (n == 2)) {
        if (r[0] < r[1]) {
            rmin = r[0];
            r[0] = r[1];
            r[1] = rmin;
        }
        return;
    }

    for (size_t i = n; --i;) {
        size_t min = i;
        rmin = r[min];
        for (size_t j = i; j--;) {
            unsigned long rj = r[j];
            if (UNLIKELY (rj < rmin)) {
                min = j;
                rmin = rj;
            }
        }
        if (LIKELY (min != i)) {
            r[min] = r[i];
            r[i] = rmin;
        }
    }
}

renumber_t::cooked renumber_t::cook(unsigned long p, std::vector<std::vector<unsigned long>> & roots) const
{
    cooked C;
    std::ostringstream os;

    for (unsigned int i = 0; i < get_nb_polys() ; i++)
        C.nroots.push_back(roots[i].size());

    if (renumber_format != renumber_format_flat) {
        for (unsigned int i = 0; i < get_nb_polys() ; i++)
            renumber_sort_ul (&roots[i][0], roots[i].size());

        /* With the traditional format, the root on ratside side becomes vp.
         * If there is no ratside side or not root on ratside side for this
         * prime (i.e., lpb < p), then the largest root becomes vp.
         */
        p_r_values_t vp = compute_vp_from_p (p);

        if (renumber_format == renumber_format_variant) {
            /* The "variant" has the advantage of being fairly simple */
            C.traditional.push_back(vp);
            for (int side = get_nb_polys() - 1; side--; ) {
                if (side == get_rational_side())
                    continue;
                for (auto r : roots[side]) {
                    p_r_side x { (p_r_values_t) p, (p_r_values_t) r, side };
                    C.traditional.push_back(compute_vr_from_p_r_side (x));
                }
            }
        } else if (get_nb_polys() == 1) {
            roots[0][0] = vp;
            for(auto x : roots[0])
                C.traditional.push_back(x);
        } else if (get_nb_polys() == 2 && get_rational_side() >= 0) {
            if (LIKELY(roots[get_rational_side()].size())) /* There is at most 1 rational root. */
                C.traditional.push_back(vp);
            else {
                int algside = 1-get_rational_side();
                roots[algside][0] = vp; /* No rational root */
                for(auto x : roots[algside])
                    C.traditional.push_back(x);
            }
        } else if (get_nb_polys() == 2) { /* Two alg polys */
            if (LIKELY (roots[1].size())) {
                roots[1][0] = vp;
                for(auto x : roots[1])
                    C.traditional.push_back(x + p + 1);
            } else {
                roots[0][0] = vp;
            }
            for(auto x : roots[0])
                C.traditional.push_back(x);
        } else { /* More than two polys (with or without ratside side). */
            bool replace_first = false;

            if (get_rational_side() == -1 || roots[get_rational_side()].empty())
                /* The largest root becomes vp */
                replace_first = true;
            else
                C.traditional.push_back(vp);

            for (int side = get_nb_polys() - 1; side--; ) {
                if (side == get_rational_side())
                    continue;

                for (auto r : roots[side]) {
                    if (UNLIKELY(replace_first)) {
                        C.traditional.push_back(vp);
                        replace_first = false;
                    } else {
                        p_r_side x { (p_r_values_t) p, (p_r_values_t) r, side };
                        C.traditional.push_back(compute_vr_from_p_r_side (x));
                    }
                }
            }
        }
        for(auto x : C.traditional)
            os << x << "\n";
        C.text = os.str();
    } else {
        for (int side = 0; side < (int) get_nb_polys(); ++side) {
            for (auto r : roots[side]) {
                p_r_side x { (p_r_values_t) p, (p_r_values_t) r, side };
                C.flat.emplace_back(
                        std::array<p_r_values_t, 2> {
                            (p_r_values_t) p,
                            compute_vr_from_p_r_side (x)
                        });
            }
        }
        for(auto x : C.flat)
            os << x[0] << " " << x[1] << "\n";
        C.text = os.str();
    }
    return C;
}

/* Only for the traditional format, when there is no rational side.
 *
 * Set x.r to the largest root of f modulo x.p such that (x.p,x.r) corresponds to an
 * ideal on side x.side which is not a bad ideal.
 * Note: If there is a projective root, it is the largest (r = p by convention)
 * Return true if such root mod p exists, else non-zero
 */
bool
renumber_t::traditional_get_largest_nonbad_root_mod_p (p_r_side & x) const
{
    mpz_poly_srcptr f = cpoly->pols[x.side];
    mpz_srcptr lc = f->coeff[f->deg];
    p_r_values_t p = x.p;
    int side = x.side;
    if (mpz_divisible_ui_p (lc, p) && !is_bad ({p, p, side})) {
        x.r = p;
        return true;
    }

    auto roots = mpz_poly_roots(cpoly->pols[side], (unsigned long) p);
    renumber_sort_ul (&roots[0], roots.size()); /* sort in decreasing order */
    for (auto r : roots) {
        if (!is_bad ({ (p_r_values_t) p, (p_r_values_t) r, side})) {
            x.r = r;
            return true;
        }
    }
    return false;
}

/* return j such that min <= j <= i, and j maximal, with
 * traditional_data[j] pointing to a vp value (i.e. one of the rare
 * values j such that traditional_data[j-1] <= traditional_data[j]
 * (in the traditional_variant mode, such increases always come in pairs.
 */
index_t renumber_t::traditional_backtrack_until_vp(index_t i, index_t min) const
{
    for( ; i > min && traditional_data[i-1] > traditional_data[i] ; --i);
    if (renumber_format == renumber_format_variant && i > min && traditional_data[i-1] > traditional_data[i])
        i--;
    return i;
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

/* i-above_bad is in [0..traditional_data.size()]
 */
bool renumber_t::traditional_is_vp_marker(index_t i) const
{
    if (i == traditional_data.size()) return true;
    if (i == 0) return true;
    if (traditional_data[i] > traditional_data[i-1]) {
        if (renumber_format == renumber_format_variant) {
            ASSERT_ALWAYS(i + 1 < traditional_data.size());
            return traditional_data[i] <= traditional_data[i + 1];
        }
        return true;
    }
    return false;
}


/* XXX This is an important part of the index_from lookup. This one needs
 * to do a dichotomy, one way or another.
 *
 * Note that we return the index relative to the internal table, which is
 * shifted by [[above_bad]] compared to the indices that we return in the
 * public interface.
 */
index_t renumber_t::get_first_index_from_p(p_r_side x) const
{
    p_r_values_t p = x.p;
    p_r_values_t side = x.side;
    if (p < index_from_p_cache.size()) {
        index_t i = index_from_p_cache[p];
        if (renumber_format == renumber_format_flat) {
            if (UNLIKELY(flat_data[i][0] != p))
                throw prime_maps_to_garbage(p, i, flat_data[i][0]);
        } else {
            p_r_values_t vp = compute_vp_from_p (p);
            if (UNLIKELY(traditional_data[i] != vp))
                throw prime_maps_to_garbage(p, i, traditional_data[i]);
        }
        return i;
    }

    if (UNLIKELY(p >> lpb[side]))
        throw prime_is_too_large(p);

    if (renumber_format == renumber_format_flat) {
        std::array<p_r_values_t, 2> p0 {{ p, 0 }};
        auto it = std::lower_bound(flat_data.begin(), flat_data.end(), p0);
        if (it == flat_data.end())
            throw prime_is_too_large(p);
        index_t i = it - flat_data.begin();
        if ((*it)[0] != p)
            throw prime_maps_to_garbage(p, i, (*it)[0]);
        return i;
    } else {
        /* A priori, traditional_data has exactly above_all-above_bad
         * entries. Except that it holds _only_ in the
         * renumber_format_traditional case, since we insert extra data
         * in the renumber_format_variant case
         */
        if (renumber_format == renumber_format_traditional)
            ASSERT_ALWAYS(above_all == above_bad + traditional_data.size());
        index_t max = traditional_data.size();
        index_t min = above_cache - above_bad;
        p_r_values_t vp = compute_vp_from_p (p);
        /* We have to look for i such that tab[i] == vp between min and max. */
        for( ; max > min ; ) {
            index_t i = min + (max - min) / 2; /* avoids overflow when
                                                  min + max >= UMAX(index_t) */
            i = traditional_backtrack_until_vp(i, min);
            if (traditional_data[i] == vp)
                return i;
            if (traditional_data[i] < vp) {
                if (UNLIKELY(i == min)) {
                    /* This is a corner case. We're below what we're
                     * looking for, sure, but we're actually looping.
                     * Need to break the corner case. We'll finish soon
                     * anyway, because our middle index landed inside the
                     * range for the smallest p, which is obviously truly
                     * exceptional.
                     */
                    i++;
                    for( ; i < max && !traditional_is_vp_marker(i) ; i++)
                        ;
                    if (traditional_data[i] == vp) return i;
                }
                min = i;
            } else {
                max = i;
            }
        }
        throw cannot_find_p(p);
    }
}

index_t renumber_t::index_from_p_r (p_r_side x) const
{
    index_t i = get_first_index_from_p(x);
    p_r_values_t vr = compute_vr_from_p_r_side(x);

    if (renumber_format == renumber_format_flat) {
        /* The "flat" format has really simple lookups */
        for( ; flat_data[i][0] == x.p ; i++) {
            if (flat_data[i][1] == vr)
                return above_bad + i;
        }
        throw cannot_find_pr(x);
    }

    p_r_values_t vp = traditional_data[i];
    /* Now i points to the beginning of data for p.
     *  tab[i] is vp.
     *  in variant format: then comes j such that the first ideal above p
     *  actually has index j
     *  then come the descriptors for all roots mod p
     */

    int outer_idx;
    if (renumber_format == renumber_format_variant) {
        /* This is the special thing about the "variant" format */
        outer_idx = above_bad + traditional_data[++i] - vp;
    } else {
        outer_idx = above_bad + i;
    }

    /* get first vr in the sequence, once we've skipped vp. Note that
     * if there's no vr at all, i might actually be the next vp already.
     */
    i++;

    /* In the "traditional" format, among all roots, the one with the
     * largest vr (which is the rational one if there is a rational side)
     * is actually missing in the table, and replaced by vp. Which makes
     * for several cases in which we want to return outer_idx
     * immediately:
     *  - an ideal on rational side always corresponds to the first
     *    element of a sequence
     *  - if i is the last index of the table
     *  - if the sequence contains only one value
     *  - if the first sequence element is less than vr
     */
    if (x.side == get_rational_side())
        return outer_idx;
    if (i == traditional_data.size())
        return outer_idx;
    if (vp < traditional_data[i])
        return outer_idx;
    if (vr > traditional_data[i])
        return outer_idx;
    /* otherwise we'll find it eventually. */
    for(unsigned int j = 0 ; traditional_data[i + j] < vp ; j++) {
        if (vr == traditional_data[i + j])
            return outer_idx + j;
    }
    throw cannot_find_pr(x, vp, vr);
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
            std::vector<int> exps(I.second.nbad, 0);

            /* work here, and return */
            for(auto J : I.second.branches) {
                int k = J.k;
                p_r_values_t pk = x.p;
                for( ; --k ; ) {
                    p_r_values_t pk1 = pk * x.p;
                    /* we want to be sure that we don't wrap around ! */
                    ASSERT_ALWAYS(pk1 > pk);
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
                modul_set_uint64 (ma, b, m);
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
                            exps.push_back(-v);
                        }
                    }
                    return { first, exps };
                }
            }
        }
        first += I.second.nbad;
    }
    throw cannot_lookup_p_a_b_in_bad_ideals(x, a, b);
}



/* additional columns _must_ be handled differently at this point */
renumber_t::p_r_side renumber_t::p_r_from_index (index_t i) const
{
    if (i < above_add) {
        for(int side = 0 ; side < (int) get_nb_polys() ; side++) {
            if (mpz_poly_is_monic(cpoly->pols[side]))
                continue;
            if (i == 0)
                return { 0, 0, side };
            i--;
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
    if (renumber_format == renumber_format_flat) {
        p_r_values_t p = flat_data[i][0];
        p_r_values_t vr = flat_data[i][1];
        return compute_p_r_side_from_p_vr(p, vr);
    }
    if (renumber_format == renumber_format_traditional) {
        index_t i0 = traditional_backtrack_until_vp(i, 0);
        index_t vr = traditional_data[i];
        index_t vp = traditional_data[i0];
        index_t p  = compute_p_from_vp(vp);
        return compute_p_r_side_from_p_vr(p, vr);
    }
    if (renumber_format == renumber_format_variant) {
        /* That is the annoying part. lookup at i is probably not right. */
        index_t i0;

        index_t max = traditional_data.size();
        index_t min = above_cache - above_bad;
        /* We have to look for i0 and i1 two consecutive vp markers
         * (possibly with i1==max) such that
         * traditional_data[i0 + 1] <= i < traditional_data[i1 + 1]
         */
        for( ; max > min ; ) {
            index_t middle = min + (max - min) / 2; /* avoids overflow */
            middle = traditional_backtrack_until_vp(middle, min);
            if (traditional_data[middle + 1] == i) {
                i0 = middle;
                break;
            } else if (traditional_data[middle + 1] < i) {
                if (UNLIKELY(middle == min)) {
                    /* This is a corner case. We're below what we're
                     * looking for, sure, but we're actually looping.
                     * Need to break the corner case. We'll finish soon
                     * anyway, because our middle index landed inside the
                     * range for the smallest p, which is obviously truly
                     * exceptional.
                     */
                    index_t j = middle;
                    for(j+=2 ; j < max && !traditional_is_vp_marker(j) ; j++)
                        ;
                    if (j == max || i < traditional_data[j + 1]) {
                        i0 = middle;
                        break;
                    }
                }
                min = middle;
            } else {
                max = middle;
            }
        }
        if (min >= max)
            throw cannot_find_i(i);
        unsigned int di = i - traditional_data[i0 + 1];
        index_t vp = traditional_data[i0];
        index_t vr = di ? traditional_data[i0 + 1 + di] : vp;
        index_t p  = compute_p_from_vp(vp);
        return compute_p_r_side_from_p_vr(p, vr);
    }
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
    uint64_t vp = vp_from_p(p, get_nb_polys(), get_rational_side() >= 0);
    if (nbits (vp) <= 32)
        return 32;
    else
        return 64;
}

void renumber_t::read_header(std::istream& is)
{
    std::istringstream iss;
    std::string s;

    if (renumber_format == renumber_format_traditional) {
        for( ; getline(is, s) && !s.empty() && s[0] == '#' ; ) ;
        iss.str(s);
        unsigned int nbits, nbad, nadd, nonmonic_bitmap, nbpol;
        int ratside;
        iss >> nbits >> ratside >> nbad >> nadd >> nonmonic_bitmap >> nbpol;
        for(auto & x : lpb) iss >> x;
        if (!iss) throw corrupted_table("cannot parse");
        if (nbits != needed_bits()
                || ratside != get_rational_side()
                || nbad != (above_bad - above_add)
                || nadd != above_add
                || nbpol != get_nb_polys())
            throw std::runtime_error("incompatible renumber table");
    } else {
        for( ; getline(is, s) && !s.empty() && s[0] == '#' ; ) ;
        iss.str(s);
        for(auto & x : lpb) iss >> x;
        if (!iss) throw corrupted_table("cannot parse");

        // we only have to parse the large prime bounds
        for( ; getline(is, s) && !s.empty() && s[0] == '#' ; ) ;
        iss.str(s);
        for(auto & x : lpb) iss >> x;
        if (!iss) throw corrupted_table("cannot parse");
    }
}

void renumber_t::read_bad_ideals(std::istream& is)
{
    ASSERT_ALWAYS (renumber_format != renumber_format_traditional);
    for(int side = 0; is && side < (int) get_nb_polys(); side++) {
        int x, n;
        is >> x >> n;
        ASSERT_ALWAYS(x == side);
        for( ; n-- ; ) {
            badideal b(is);
            p_r_side x { (p_r_values_t) mpz_get_ui(b.p), (p_r_values_t) mpz_get_ui(b.r), side };
            bad_ideals.emplace_back(x, b);
        }
    }
}

void renumber_t::write_header(std::ostream& os) const
{
    /* The traditional format doesn't even accept comments in the very
     * first line !
     */
    if (renumber_format == renumber_format_traditional) {
        // the first line
        unsigned long nonmonic_bitmap = 0;
        for (unsigned int i = get_nb_polys(); i-- ; ) {
            nonmonic_bitmap <<= 1;
            nonmonic_bitmap += !mpz_poly_is_monic(cpoly->pols[i]);
        }
        os << needed_bits()
            << " " << get_rational_side()
            << " " << above_bad - above_add
            << " " << above_add
            << " " << nonmonic_bitmap
            << " " << get_nb_polys();
        for(auto x : lpb)
            os << " " << x;
        os << "\n";
    }

    // the comment is never significant. However only the newer formats
    // give us the opportunity to stick in a marker for the file format.
    // So if we want (for the moment) to provide old-format renumber
    // files that can be parsed by old code, we can only convey the
    // format information as an (unparsed) side note.
    os << "# Renumber file using format " << renumber_format;

    /* Write the polynomials as comments */
    for (unsigned int i = 0; i < get_nb_polys() ; i++)
        os << "# pol" << i << ": " << cpoly->pols[i] << "\n";

    if (renumber_format != renumber_format_traditional) {
        os << renumber_format << "\n";
        // number of additional columns is implicit anyway.
        os << "# large prime bounds:\n";
        for (unsigned int i = 0; i < get_nb_polys() ; i++) {
            if (i) os << " ";
            os << lpb[i];
        }
        os << "\n";
    }
}

void renumber_t::write_bad_ideals(std::ostream& os) const
{
    /* Write first the bad ideal information at the beginning of file */
    for(int side = 0; os && side < (int) get_nb_polys(); side++) {
        os << "# bad ideals on side" << side << std::endl;
        unsigned int n = 0;
        for(auto const & b : bad_ideals)
            if (b.first.side == side) n++;
        os << side << " " << n << std::endl;
        for(auto const & b : bad_ideals) {
            if (b.first.side == side)
                os << b.second << std::endl; 
        }
    }

    os << "# renumber table for all indices above " << above_bad << ":\n";
}

void renumber_t::use_additional_columns_for_dl()
{
    ASSERT_ALWAYS(above_all == 0);
    above_add = 0;
    for(unsigned int side = 0 ; side < get_nb_polys() ; side++)
        above_add += !mpz_poly_is_monic(cpoly->pols[side]);
    above_bad = above_add;
    above_cache = above_add;
    above_all = above_add;
}

void renumber_t::compute_bad_ideals()
{
    for(int side = 0 ; side < (int) get_nb_polys() ; side++) {
        cxx_mpz_poly f(cpoly->pols[side]);
        if (f->deg == 1) continue;
        std::vector<badideal> badideals = badideals_for_polynomial(f, side);

        for(auto const & b : badideals_for_polynomial(f, side)) {
            p_r_values_t p = mpz_get_ui(b.p);
            p_r_values_t r = mpz_get_ui(b.r);
            p_r_side x { p, r, side };
            bad_ideals.emplace_back(x, b);
        }
    }
}

index_t renumber_t::use_cooked_nostore(index_t n0, cooked const & C)
{
    for(auto n : C.nroots)
        n0 += n;
    return n0;
}

void renumber_t::use_cooked(cooked const & C)
{
    traditional_data.insert(traditional_data.end(), C.traditional.begin(), C.traditional.end());
    flat_data.insert(flat_data.end(), C.flat.begin(), C.flat.end());
    above_all = use_cooked_nostore(above_all, C);
}

void renumber_t::read_table(std::istream& is)
{
    if (renumber_format == renumber_format_flat) {
        for(p_r_values_t p, r ; is >> p >> r ; ) {
            flat_data.emplace_back(std::array<p_r_values_t, 2> {{ p, r }});
        }
    } else {
        for(p_r_values_t v ; is >> v ; )
            traditional_data.push_back(v);
    }
}

renumber_t::renumber_t(const char * filename)
{
    ifstream_maybe_compressed is(filename);
    read_header(is);
    read_bad_ideals(is);
    read_table(is);
}

void renumber_t::read_bad_ideals_info(std::istream & is)
{
    /* the badidealinfo file of old is slightly annoying to deal with.
     */
    bad_ideals.clear();
    std::map<p_r_side, badideal> met;
    for( ; is ; ) {
        std::string s;
        std::istringstream iss;
        for( ; getline(is, s) && !s.empty() && s[0] == '#' ; ) ;
        p_r_values_t p, rk;
        int k, side;
        is >> p >> k >> rk >> side;
        p_r_values_t r = mpz_get_ui(badideal::r_from_rk(p, k, rk)); 
        p_r_side x { p, r, side };
        badideal b(p, r);
        badideal::branch br;
        br.k = k;
        br.r = rk;
        for(int e ; is >> e ; br.v.push_back(e));
        b.nbad = br.v.size();
        auto it = met.find(x);
        if (it != met.end() && it->second.nbad != b.nbad) {
            std::ostringstream os;
            os << "badidealinfo file is bad ;"
                << " valuation vector found in branch description"
                << " is not consistent above"
                << "(" << p << ", " << r << ", side " << side << ")";
            throw std::runtime_error(os.str());
        } else if (it == met.end()) {
            met[x] = b;
        }
        met[x].branches.emplace_back(std::move(br));
    }
    for(auto const & prb : met) {
        bad_ideals.emplace_back(prb);
    }
}

renumber_t::renumber_t(const char * filename, const char * badidealinfofile)
{
    ifstream_maybe_compressed is(filename);
    read_header(is);
    std::ifstream isi(badidealinfofile);
    read_bad_ideals_info(isi);
    read_table(is);
}

std::string renumber_t::debug_data(index_t i) const
{
    p_r_side x = p_r_from_index (i);
    std::ostringstream os;

    os << "i=" << i;

    if (is_additional_column (i)) {
        os << " tab[i]=#"
            << " added column for side " << x.side;
    } else if (is_bad(i)) {
        os << " tab[i]=#"
            << " bad ideal"
            << " (" << bad_ideals[i - above_add].second.nbad << " branches)"
            << " above"
            << " (" << x.p << "," << x.r << ")"
            << " on side " << x.side;
    } else {
        os << " tab[i]=";
        i -= above_bad;
        if (renumber_format == renumber_format_flat) {
            os << " (" << flat_data[i][0] << "," << flat_data[i][1] << ")";
        } else {
            os << traditional_data[i];
        }
        os << " p=%" << x.p;
        if (x.side == get_rational_side()) {
            os << " rat";
        } else {
            os << " r=" << x.r;
            if (x.r == x.p)
                os << " proj";
        }
        os << " side " << x.side;
    }

    return os.str();
}
