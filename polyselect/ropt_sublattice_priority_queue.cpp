#include "cado.h" // IWYU pragma: keep

#include <cstddef>

#include <utility>

#include <gmp.h>

#include "ropt_sublattice_priority_queue.h"
#include "cxx_mpz.hpp"
#include "min_max_heap.hpp"
/**
 * Priority queue for sublattices over a product of p^e.
 */

struct sublattice_info {
    cxx_mpz u, v;
    cxx_mpz modulus;
    sublattice_info() = default;
    sublattice_info(mpz_srcptr u, mpz_srcptr v, mpz_srcptr modulus) {
        mpz_set(this->u, u);
        mpz_set(this->v, v);
        mpz_set(this->modulus, modulus);
    }
    sublattice_info(cxx_mpz const & u, cxx_mpz const & v, cxx_mpz const & modulus):
        u(u), v(v), modulus(modulus) {}
    sublattice_info(cxx_mpz && u, cxx_mpz && v, cxx_mpz && modulus):
        u(u), v(v), modulus(modulus) {}
    bool operator<(sublattice_info const & o) const {
        if (u < o.u) return true;
        if (u > o.u) return false;
        if (v < o.v) return true;
        if (v > o.v) return false;
        if (modulus < o.modulus) return true;
        if (modulus > o.modulus) return false;
        return false;
    }
};

struct sublattice_pq_impl : min_max_heap<std::pair<float, sublattice_info>> {
    size_t max_count = 0;
};

void sublattice_priority_queue_init(sublattice_priority_queue_ptr q,
                         size_t max_len)
{
    auto * qi = new sublattice_pq_impl;
    if (max_len == 0)
        max_len = 1;
    qi->max_count = max_len;
    q->impl = (void *) qi;

}

void sublattice_priority_queue_push(sublattice_priority_queue_ptr q,
                            mpz_srcptr u,
                            mpz_srcptr v,
                            mpz_srcptr mod, 
                            float val )
{
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast)
    auto * qi = (sublattice_pq_impl *) q->impl;
    qi->push(std::pair<float, sublattice_info>(val, { u, v, mod }));
    if (qi->size() > qi->max_count)
        qi->popMin();
}

void sublattice_priority_queue_clear (sublattice_priority_queue_ptr q)
{
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast)
    auto * qi = (sublattice_pq_impl *) q->impl;
    delete qi;
}

size_t sublattice_priority_queue_size (sublattice_priority_queue_srcptr q)
{
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast)
    auto const * qi = (const sublattice_pq_impl *) q->impl;
    return qi->size();
}

int sublattice_priority_queue_empty (sublattice_priority_queue_srcptr q)
{
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast)
    auto const * qi = (const sublattice_pq_impl *) q->impl;
    return qi->empty();
}

/* This calls the fucntion f on all the queue elements, in order. We
 * could as well pass each element's rank and score, but presently we do
 * not.
 */
void sublattice_priority_queue_do(sublattice_priority_queue_srcptr q, void (*f)(mpz_srcptr u, mpz_srcptr v, mpz_srcptr modulus, void * arg), void *arg)
{
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast)
    auto const & qi = * (const sublattice_pq_impl *) q->impl;
    /* copy, and consume the copy */
    sublattice_pq_impl cpi = qi;
    for( ; !cpi.empty() ; ) {
        auto x = cpi.popMax();
        (*f)(x.second.u, x.second.v, x.second.modulus, arg);
    }
}

