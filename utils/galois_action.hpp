#ifndef CADO_UTILS_GALOIS_ACTION_HPP
#define CADO_UTILS_GALOIS_ACTION_HPP

#include <iostream>
#include <string>

#include "renumber.hpp"     // for renumber_t
#include "typedefs.h"       // for index_t
#include "fmt/format.h"     // for fmt::formatter
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
 * After the matrix is applied it may be necessary to negate a' and b' in order
 * to vector an ab pair b' >= 0.
 */

/* Interface for galois action implementation */
class galois_action_impl_base
{
public:
    galois_action_impl_base() = default;
    galois_action_impl_base(const galois_action_impl_base&) = delete;
    galois_action_impl_base& operator=(const galois_action_impl_base&) = delete;

    virtual unsigned int get_order() const = 0;
    virtual unsigned long apply(unsigned long r, unsigned long p) const = 0;
    virtual uint64_t hash_ab(int64_t a, uint64_t b,
                             uint64_t CA, uint64_t CB) const = 0;

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
    galois_action(const std::string &action);

    /* Return the order of the automorphism */
    unsigned int get_order() const;

    /* Assume 0 <= r <= p
     * Return (alpha*r+beta)/(gamma*r+delta) modulo p       if r != p
     *        alpha/gamma modulo p                          if r == p
     * Note: oo is represented by having r = p for the input and by returning p
     * for the output.
     */
    unsigned long apply(unsigned long r, unsigned long p) const;

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

    ~galois_action();

private:
    using impl_ptr = galois_action_impl_base *;
    impl_ptr impl;
};

namespace fmt {
    template <> struct formatter<galois_action>: ostream_formatter {};

}

#endif  /* CADO_UTILS_GALOIS_ACTION_HPP */
