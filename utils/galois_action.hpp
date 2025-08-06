#ifndef CADO_UTILS_GALOIS_ACTION_HPP
#define CADO_UTILS_GALOIS_ACTION_HPP

#include <cstdint>
#include <cstddef>

#include <iostream>
#include <string>
#include <vector>
#include <memory>

#include "renumber.hpp"
#include "special-q.hpp"
#include "typedefs.h"
#include "fmt/base.h"
#include "fmt/ostream.h"

/*
 * Implementation of Galois action.
 * A galois action is represented by a 2x2 integral matrix denoted
 * [[alpha, beta], [gamma, delta]].
 * It corresponds to the action:
 *      (u : v) -> (alpha*u+beta*v : gamma*u+delta*v) in the projective plane
 *
 * When applied to a single value r modulo p (as in the method `apply` below),
 * the integer r represents the vector (r : 1) if r < p and (0 : 1) if r == p.
 * The same convention applies to the output.
 *
 * When applied to an (a,b) pair representing a relation, it corresponds to the
 * action of the 2x2 matrix [[delta, -beta], [-gamma, alpha]] (i.e., the inverse
 * of the above matrix). It means that
 *      (a, b) ->  (a',b') := (delta*a-beta*b, -gamma*a+alpha*b)
 * After the matrix is applied it may be necessary to negate a' and b' in
 * order to have b' >= 0.
 */

/* Interface for galois action implementation */
class galois_action_impl_base
{
public:
    galois_action_impl_base() = default;
    galois_action_impl_base(const galois_action_impl_base&) = default;
    galois_action_impl_base& operator=(const galois_action_impl_base&) = default;
    galois_action_impl_base(galois_action_impl_base &&) = default;
    galois_action_impl_base& operator=(galois_action_impl_base &&) = default;

    virtual unsigned int get_order() const = 0;
    virtual unsigned long apply(unsigned long r, unsigned long p) const = 0;
    virtual special_q apply(special_q const & q) const = 0;
    virtual uint64_t hash_ab(int64_t a, uint64_t b,
                             uint64_t CA, uint64_t CB) const = 0;

    virtual int apply_ab(int64_t &, uint64_t &) const = 0;
    virtual unsigned long apply_ab_cofactor(unsigned long r, unsigned long p) const = 0;

    virtual void print_action(std::ostream& os) const = 0;
    virtual void print_name(std::ostream& os) const = 0;

    virtual ~galois_action_impl_base() = default;
};

/* class representing Galois action */
class galois_action
{
public:
    /* Init the class using a predefined name ("none", "id", "_y", "1/y",
     * "autom3.1g", ...).
     */
    explicit galois_action(const std::string & action);

    /* Same as above but with const char *. Passing nullptr is okay and
     * corresponds to the galois action "none".
     */
    explicit galois_action(char const * action)
        : galois_action(std::string(action == nullptr ? "none" : action)) {
    }

    /* Default constructor corresponds the galois action "none" */
    galois_action() : galois_action("none") {
    }

    /* Return the order of the automorphism */
    unsigned int get_order() const;

    /* Return true for actual galois action (i.e., with order > 1) */
    explicit operator bool() const { return get_order() > 1; }

    /* Assume 0 <= r <= p
     * Return (alpha*r+beta)/(gamma*r+delta) modulo p       if r != p
     *        alpha/gamma modulo p                          if r == p
     * Note: oo is represented by having r = p for the input and by returning p
     * for the output.
     *
     * TODO: We rely on the fact that p is not a (non-trivial) prime
     * power, here. Dealing with the wider range of possible values would
     * require some extra work.
     */
    unsigned long apply(unsigned long r, unsigned long p) const;

    /* This won't work with composite special-q, currently */
    special_q apply(special_q const & q) const;

    /* This applies to a,b pairs, not to roots. If a-b*alpha is some
     * element with coprime a,b and b>0, replace (a,b) by (a',b') 
     * such that sigma(a-b*alpha) = a-b*sigma(alpha) is also essentially
     * (a'-b'*alpha), but up to a multiplier. More specifically, the
     * multiplier is the product of something that is constant for the
     * Galois action, and some integer m that we return. The equation
     * is:
     * a-b*sigma(alpha) = m * (a'-b'*alpha) / c_G
     * m may be 1, -1, 2, -2, 3, -3. c_G can be typically 1, alpha,
     * alpha-1
     *
     * for any prime p and any r in Z/p, encoded by an integer in
     * [0,p-1], we must have:
     *
     * a - b * apply(r, p) == m * (a' - b' * r) / apply_ab_cofactor(r, p)
     *
     * (note though that apply_ab_cofactor(r, p) might be zero)
     */
    int apply_ab(int64_t &, uint64_t &) const;

    unsigned long apply_ab_cofactor(unsigned long r, unsigned long p) const;

    /* Compute the hash of the (a,b) pair (using CA and CB as hash constants).
     * The returned value is the same for all (a, b) pairs in the orbit of the
     * galois action.
     * Assumes that a and b are coprime (in fact assumes that a != 0, b != 0,
     * that |a| == b implies (a,b) = (-1, 1) or (1, 1) and that a,b are not both
     * even).
     */
    uint64_t hash_ab(int64_t a, uint64_t b, uint64_t CA, uint64_t CB) const;

    /* Store into sigma the full action of the galois action on the index of the
     * ideals. For each orbit of ideals (I0, ..., Ik) with index (i0, ..., ik),
     * a single representative is chosen (say ij) and sigma[i0], ..., sigma[ik]
     * is set to ij.
     * Return the number of orbits found.
     */
    size_t compute_action_on_index(std::vector<index_t> &sigma,
                                   renumber_t const & tab) const;

    friend std::ostream& operator<<(std::ostream &os, const galois_action &g);

private:
    std::unique_ptr<galois_action_impl_base> impl;
};

namespace fmt {
    template <> struct formatter<galois_action>: ostream_formatter {};
}

#endif  /* CADO_UTILS_GALOIS_ACTION_HPP */
