#ifndef SIEVE_LAS_SPECIAL_Q_TASK_HPP_
#define SIEVE_LAS_SPECIAL_Q_TASK_HPP_

#include <utility>

#include "special-q.hpp"
// #include "relation.hpp"

struct las_info;

/* A special-q task is either a special_q_task_simple or a
 * special_q_task_tree.
 */
struct special_q_task : public special_q {
    special_q const & sq() const { return *this; }

    enum status_code {
        PENDING,
        IN_PROGRESS,
        IN_RECURSION,
        DONE,
        ABANDONED
    };

    status_code status = PENDING;

    protected:
    special_q_task() = default;
    explicit special_q_task(special_q q)
        : special_q(std::move(q))
    { }


    public:
    special_q_task(special_q_task const &) = default;
    special_q_task(special_q_task &&) = default;
    special_q_task& operator=(special_q_task const &) = default;
    special_q_task& operator=(special_q_task &&) = default;

    // NOLINTBEGIN(readability-convert-member-functions-to-static)
    int depth() const { return 0; }
    // NOLINTEND(readability-convert-member-functions-to-static)

    virtual void update_status(status_code, status_code) = 0;
    virtual bool must_take_decision() const = 0;
    virtual ~special_q_task() = default;
};
#endif	/* SIEVE_LAS_SPECIAL_Q_TASK_HPP_ */
