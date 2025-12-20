#ifndef UTILS_WORK_QUEUE_HPP_
#define UTILS_WORK_QUEUE_HPP_

#include <cstddef>

#include <mutex>
#include <condition_variable>
#include <vector>
#include <list>
#include <tuple>
#include <thread>
#include <memory>
#include <utility>
#include <type_traits>
#include <concepts>
#include <stdexcept>

#include "antilock.hpp"

/* Type-safe implemntation of a work queue.
 *
 * This is a proof-of-concept that can serve as a drop-in replacement for
 * the simple work queues that we had in sqrt/crtalgsqrt.cpp. However,
 * the more interesting goal is to replace utils/threadpool.hpp
 *
 * Any function can be used as a task. No return types are handled though
 * (use references to variables if you wish). Of course, all parameters
 * mut remain live until the task is complete.
 *
 * Tasks are processed in order: if task a is pushed before task b, then
 * task a will be picked up by a worker thread before task b is. This
 * does not say anything of the completion times, of course.
 *
 * Tasks derive from the task_base type, which provides the bare
 * essentials.
 *
 * A task can be waited for by calling the join_task() method (basically,
 * it is the only method that makes sense using on task_base).
 *
 * A task is modeled as a shared_ptr<task_base> (aka
 * work_queue::task_handle). The reasons why we use a shared_ptr and not
 * a unique_ptr are:
 *  - The task must remain in the work queue's data structures for
 *    bookkeeping.
 *  - The control flow that pushes the task quite often wants to wait for
 *    its completion. For this, returning the shared_ptr seems
 *    appropriate.
 *
 */

namespace cado {

struct work_queue {
    struct task_base {
        private:
            bool done = false;
            std::mutex m_;
            std::condition_variable c_;
        public:
            /* default construction is possible */
            task_base() = default;
            /* because there is a mutex, we don't want the task to move.  */
            task_base(task_base const&) = delete;
            task_base(task_base &&) = delete;
            task_base& operator=(task_base const&) = delete;
            task_base& operator=(task_base &&) = delete;
            /* a virtual dtor is good to have! */
            virtual ~task_base() = default;

            virtual void call() = 0;

            /* wait until this task is complete */
            void join_task() {
                std::unique_lock<std::mutex> u(m_);
                for( ; !done ; )
                    c_.wait(u);
            }
        protected:
            /* used by worker threads to mark this task as complete */
            void complete() {
                const std::unique_lock<std::mutex> u(m_);
                done = true;
                c_.notify_all();
            }
    };

    /* this is the type-safe derived type that calls a given function.
     * The function and its argument are stored in this object until the
     * object is consumed by call()
     */
    template<typename Callable, typename... Args>
        struct task : public task_base {
            Callable f;
            std::tuple<Args...> args;
            /* reasonably standard mechanics to expand a tuple into a
             * sequence of arguments */
            template<size_t... S>
                void call_internal(const std::index_sequence<S...>) {
                    f(std::get<S>(args)...);
                }
            void call() override {
                call_internal(std::make_index_sequence<sizeof...(Args)>());
                task_base::complete();
            }
            explicit task(Callable&& f, Args&&... args)
                : f(std::forward<Callable>(f))
                , args(std::forward<Args>(args)...)
            {
            }
        };

    /* the mutex guards accesses to the task list. The worker list needs
     * no mutex */
    std::mutex m;

    /* this is used as a "new work" signal */
    std::condition_variable c;

    std::vector<std::thread> workers;

    using task_handle = std::shared_ptr<task_base>;

    /* obviously, each task must remain alive somewhere in memory until it
     * is picked up.
     */
    std::list<task_handle> tasks;

    template<typename Callable, typename... Args>
        task_handle push_task(Callable && f, Args... args)
        requires std::invocable<Callable, Args...>
        {
            const std::lock_guard<std::mutex> u(m);
            auto t = std::make_shared<task<Callable, std::decay_t<Args>...>>(std::forward<Callable>(f), std::decay_t<Args>(args)...);
            tasks.emplace_back(t);
            c.notify_one();
            return t;
        }

    /* Do a task from the work queue if there is one available. This can
     * be used to "steal" work from the workers if there are workers
     * available. It also makes it possible to use these work queues
     * without workers, if some threads are regularly consuming tasks
     * with this call.
     *
     * return true if a task was called.
     */
    bool do_one_task() {
        task_handle t;
        {
            const std::lock_guard<std::mutex> u(m);
            if (tasks.empty())
                return false;
            if (!tasks.front()) {
                /* If we have a nullptr here, it's an end marker that is
                 * meant to cause termination of registered workers.
                 * There's no reason to consume it here.
                 */
                return false;
            }
            t = std::move(tasks.front());
            tasks.pop_front();
        }
        t->call();
        return true;
    }

    /* This version correctly deals with the use case where we hold a
     * synchronization lock at call time. The unique_lock is on a mutex
     * that is different from this->m.
     *
     * It a task is available on the work queue, it is atomically removed
     * from the queue (with this->m held). If this is the case:
     *  - the invocable f is called while u is still locked. [buggy,
     *  dropped]
     *  - u is temporarily unlocked while the task is called.
     *
     * return true if a task was called.
     */
    bool do_one_task(std::unique_lock<std::mutex> & u
            /* , std::invocable auto const & f = [](){} */) {
        task_handle t;
        {
            const std::lock_guard<std::mutex> u(m);
            if (tasks.empty())
                return false;
            if (!tasks.front()) {
                /* If we have a nullptr here, it's an end marker that is
                 * meant to cause termination of registered workers.
                 * There's no reason to consume it here.
                 */
                return false;
            }
            t = std::move(tasks.front());
            tasks.pop_front();
        }
        //  f();
        /* temporarily release u in an exception-safe way */
        const cado::antilock a(u);
        t->call();
        return true;
    }

    /* the worker thread picks up work, and terminates when it receives
     * a null pointer */
    void worker() {
        for(;;) {
            task_handle t;
            {
                std::unique_lock<std::mutex> u(m);
                for( ; tasks.empty() ; )
                    c.wait(u);
                t = std::move(tasks.front());
                tasks.pop_front();
            }
            if (!t) break;
            t->call();
        }
    }

    explicit work_queue(size_t n = 0)
    {
        add_workers(n);
    }
    void add_workers(size_t n)
    {
        for(size_t i = 0 ; i < n ; i++)
            workers.emplace_back([this](){worker();});
    }

    /* We can't have any copy (of course), but move must also be
     * disabled: the worker threads have a pointer to this, so we'll go
     * right to catastrophe if we move
     */
    work_queue(work_queue const&) = delete;
    work_queue(work_queue &&) = delete;
    work_queue& operator=(work_queue const&) = delete;
    work_queue& operator=(work_queue &&) = delete;

    /* destruction is harder than it seems. We must make sure that all
     * workers complete. For this, we post as many end-of-tasks markers
     * as we have threads, from which it follows that all threads will be
     * ready to be joined.
     */
    ~work_queue() {
        {
            const std::lock_guard<std::mutex> u(m);
            for(auto const & w [[maybe_unused]] : workers) {
                /* add an empty task. It will be recognized by the worker */
                tasks.emplace_back();
                c.notify_one();
            }
        }
        for(auto & w : workers)
            w.join();
    }
};
} /* namespace cado */

#endif	/* UTILS_WORK_QUEUE_HPP_ */

