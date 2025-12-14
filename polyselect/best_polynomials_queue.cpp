#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cmath>
#include <cstddef>

#include <vector>
#include <sstream>
#include <utility>

#include "best_polynomials_queue.h"
#include "min_max_heap.hpp"
#include "cado_poly.h"

template<typename T>
struct compare_by_score {
    using value_type = std::pair<double, T>;
    bool operator()(value_type const & a, value_type const & b) const {
        return a.first < b.first;
    }
};

struct best_polynomials_queue_impl {
    size_t max_count = 0;
    using queue_value_type = compare_by_score<cxx_cado_poly>::value_type;
    min_max_heap<queue_value_type, std::vector<queue_value_type>, compare_by_score<cxx_cado_poly>> q;
};

void best_polynomials_queue_init(best_polynomials_queue_ptr b, int queue_lenth)
{
    auto * pi = new best_polynomials_queue_impl;
    pi->max_count = queue_lenth;
    b->pimpl = (void *) pi;
}

void best_polynomials_queue_clear(best_polynomials_queue_ptr b)
{
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast)
    auto * pi = (best_polynomials_queue_impl *) b->pimpl;
    delete pi;
}

double best_polynomials_queue_get_best_score(best_polynomials_queue_srcptr b)
{
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast)
    auto const & pi = * (best_polynomials_queue_impl const *) b->pimpl;
    if (pi.q.empty())
        return NAN;
    return pi.q.findMax().first;
}

double best_polynomials_queue_get_worst_score(best_polynomials_queue_srcptr b)
{
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast)
    auto const & pi = * (best_polynomials_queue_impl const *) b->pimpl;
    if (pi.q.empty())
        return NAN;
    return pi.q.findMin().first;
}

size_t best_polynomials_queue_get_count(best_polynomials_queue_srcptr b)
{
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast)
    auto const & pi = * (best_polynomials_queue_impl const *) b->pimpl;
    return pi.q.size();
}

void best_polynomials_queue_try_push(best_polynomials_queue_ptr b, cado_poly_srcptr cpoly, double score)
{
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast)
    auto & pi = * (best_polynomials_queue_impl *) b->pimpl;
    best_polynomials_queue_impl::queue_value_type v;
    v.first = score;
    cado_poly_set(v.second, cpoly);
    pi.q.push(v);
    if (pi.q.size() > pi.max_count) {
        pi.q.popMin();
    }
}

void best_polynomials_queue_print(best_polynomials_queue_srcptr b, FILE * f, const char * prefix)
{
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast)
    auto const & pi = * (best_polynomials_queue_impl const *) b->pimpl;
    /* copy, and consume the copy */
    best_polynomials_queue_impl cpi = pi;
    for(int i = 0 ; !cpi.q.empty() ; i++) {
        auto x = cpi.q.popMax();
        std::ostringstream os;
        os << x.second;
        fprintf(f, "%s%d-th best polynomial (score = %.2f):%s\n", prefix, i, x.first, os.str().c_str());
    }
}

void best_polynomials_queue_do(best_polynomials_queue_srcptr b, void (*f)(int, double, cado_poly_ptr, void *), void *arg)
{
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-cstyle-cast)
    auto const & pi = * (best_polynomials_queue_impl const *) b->pimpl;
    /* copy, and consume the copy */
    best_polynomials_queue_impl cpi = pi;
    for(int i = 0 ; !cpi.q.empty() ; i++) {
        auto x = cpi.q.popMax();
        (*f)(i, x.first, x.second, arg);
    }
}


