#ifndef SIEVE_LAS_SPECIAL_Q_TASK_SIMPLE_HPP_
#define SIEVE_LAS_SPECIAL_Q_TASK_SIMPLE_HPP_

#include <utility>

#include "special-q.hpp"
#include "las-special-q-task.hpp"
#include "macros.h"
#include "relation.hpp"

struct las_info;

/* A "simple" special-q task is what we use most of the time. It's
 * basically just a special-q, with a status field (perhaps we want to
 * add some statistics fields to it as well, starting with the time
 * spent, at least).
 *
 * This exists in order to provide a simple mirror of what we have for
 * the much more complex special_q_task_tree, which is used in the
 * descent. So most of the interfaces here look thin or even pointless,
 * but that is on purpose.
 *
 * As in the special_q_task_tree case, we have a state machine of
 * statuses, although it's much simpler because we have no recursion
 * involved.
 *
 * pending -> in_progress -> {done, abandoned}
 *
 */
struct special_q_task_simple : public special_q_task {

    explicit special_q_task_simple(special_q q)
        : special_q_task(std::move(q))
          {
          }

    private:
    friend struct special_q_task_collection;
    special_q_task_simple() = default;
    public:

    // NOLINTBEGIN(readability-convert-member-functions-to-static)
    bool must_take_decision() const override { return false; }
    // NOLINTEND(readability-convert-member-functions-to-static)

    void update_status(status_code before, status_code after) override {
        ASSERT_ALWAYS(status == before);
        status = after;
    }

};


#endif	/* SIEVE_LAS_SPECIAL_Q_TASK_SIMPLE_HPP_ */
