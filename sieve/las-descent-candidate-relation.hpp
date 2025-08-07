#ifndef SIEVE_LAS_DESCENT_CANDIDATE_RELATION_HPP_
#define SIEVE_LAS_DESCENT_CANDIDATE_RELATION_HPP_

#include <cmath>

#include <algorithm>
#include <vector>
#include <utility>

#include "special-q.hpp"
#include "relation.hpp"
#include "timing.h"

#ifdef isfinite
/* isfinite is c99 and std::isfinite is c++11 ; it's not totally clear
 * that #include <cmath> + accessing std::isfinite works.
 *
 * Under some conditions, we can get #define'd C functions, which
 * obviously invalidate the C++ prototype (icc version 16.0.3 based on
 * gcc-6.1.0 on CentOS 7.2.1511)
 */
#undef isfinite
#endif

/* For descent mode: we compute the expected time to finish given the
 * factor sizes, and deduce a deadline.  Assuming that not all
 * encountered factors are below the factor base bound, if we expect
 * an additional time T to finish the decomposition, we keep looking
 * for a better decomposition for a grace time which is computed as
 * x*T, for some configurable ratio x (one might think of x=0.2 for
 * instance. x is the grace_time_ratio member), which defines a
 * ``deadline'' for next step.  [If all factors happen to be smooth,
 * the deadline is immediate, of course.] If within the grace period,
 * a new relation is found, with an earlier implied deadline, the
 * deadline is updated. We say that the "decision is taken" when the
 * deadline passes, and the las machinery is told to decide that it
 * should proceed with the descent, and stop processing the current
 * special-q.
 */
struct descent_candidate_relation {
    relation rel;
    std::vector<special_q> outstanding;
    double time_left = INFINITY;
    double deadline = INFINITY;
    // bool marked_taken;      /* false until we take the decision */
    descent_candidate_relation() = default;
    descent_candidate_relation& operator=(descent_candidate_relation const& o) {
        if (this == &o) return *this;
        /* nothing very fancy, except that we keep the old deadline.
         * */
        rel = o.rel;
        outstanding = o.outstanding;
        time_left = o.time_left;
        deadline = std::min(deadline, o.deadline);
        return *this;
    }
    descent_candidate_relation(descent_candidate_relation const & o) = default;
    descent_candidate_relation& operator=(descent_candidate_relation && o) noexcept {
        rel = o.rel;
        outstanding = o.outstanding;
        time_left = o.time_left;
        deadline = std::min(deadline, o.deadline);
        return *this;
    }
    descent_candidate_relation(descent_candidate_relation && o) noexcept = default;
    ~descent_candidate_relation() = default;

    bool operator<(descent_candidate_relation const& b) const
    {
        if (!rel) return false;
        if (!b.rel) return true;
        if (std::isfinite(time_left)) { return time_left < b.time_left; }
        return outstanding.size() < b.outstanding.size();
    }
    explicit operator bool() const { return (bool) rel; }
    bool wins_the_game() const {
        return (bool) rel && (outstanding.empty() || seconds() >= deadline);
    }
    void set_time_left(double t, double grace_time_ratio) {
        time_left = t;
        if (outstanding.empty()) {
            deadline = seconds();
        } else {
            deadline = seconds() + grace_time_ratio * t;
        }
    }
};



#endif	/* SIEVE_LAS_DESCENT_CANDIDATE_RELATION_HPP_ */
