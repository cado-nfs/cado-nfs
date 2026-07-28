#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cstdint>

#include <condition_variable>
#include <exception>
#include <functional>
#include <mutex>
#include <queue>
#include <thread>
#include <utility>
#include <vector>

#include "macros.h"
#include "threadpool.hpp"
#include "timing.h"

worker_thread::worker_thread(thread_pool & _pool, size_t const _preferred_queue,
                             bool several_threads)
    : pool(_pool)
    , preferred_queue(several_threads ? _preferred_queue : SIZE_MAX)
{
    if (!several_threads) {
        return;
    }
    thread = std::jthread([this]() {
        pool.thread_work_on_tasks(*this);
    });
}

int worker_thread::rank() const
{
    return static_cast<int>(this - pool.threads.data());
}
int worker_thread::nthreads() const
{
    return static_cast<int>(pool.threads.size());
}
bool worker_thread::is_synchronous() const
{
    return pool.is_synchronous();
}

thread_pool::thread_pool(size_t const nr_threads, double & store_wait_time,
                         size_t const nr_queues, bool sync_thread_pool)
    : tasks(nr_queues)
    , exceptions(nr_queues)
    , created(nr_queues, 0)
    , joined(nr_queues, 0)
    , kill_threads(false)
    , store_wait_time(store_wait_time)
{
    ASSERT_ALWAYS(nr_threads == 1 || !sync_thread_pool);
    threads.reserve(nr_threads);
    for (size_t i = 0; i < nr_threads; i++)
        threads.emplace_back(*this, 0, !is_synchronous());
};

thread_pool::~thread_pool()
{
    drain_all_queues();
    {
        auto lock = get_lock();
        kill_threads = true;
        for (auto & T: tasks)
            T.not_empty.notify_all();
    }
    drain_all_queues();
    threads.clear();
    for (auto const & T: tasks)
        ASSERT_ALWAYS_NOTHROW(T.empty());
    for (auto const & E: exceptions)
        ASSERT_ALWAYS_NOTHROW(E.empty());
    store_wait_time += cumulated_wait_time;
}

void thread_pool::thread_work_on_tasks(worker_thread & I)
{
    ASSERT_ALWAYS(!is_synchronous());
    double tt = -seconds_thread();
    while (true) {
        size_t queue = I.preferred_queue;
        thread_task task = get_task(queue);
        if (task.is_terminal())
            break;

        tt += seconds_thread();
        try {
            task(&I);
        } catch (...) {
            std::scoped_lock guard(pool_mutex);
            exceptions[queue].push(std::current_exception());
        }

        {
            auto lock = get_lock();
            joined[queue]++;
            tasks[queue].task_done.notify_all();
        }
        tt -= seconds_thread();
    }
    tt += seconds_thread();
    std::scoped_lock const dummy(mm_cumulated_wait_time);
    cumulated_wait_time += tt;
}

bool thread_pool::all_task_queues_empty() const
{
    for (auto const & T: tasks)
        if (!T.empty())
            return false;
    return true;
}

void thread_pool::enqueue_task(std::function<void(worker_thread*)> task_fn, int id, size_t queue, double cost) 
{
    if (is_synchronous()) {
        created[queue]++;
        try {
            task_fn(threads.data());
        } catch (...) {
            const std::scoped_lock guard(pool_mutex);
            exceptions[queue].push(std::current_exception());
        }
        joined[queue]++;
        return;
    }

    ASSERT_ALWAYS(queue < tasks.size());

    auto lock = get_lock();
    ASSERT_ALWAYS(!kill_threads);

    tasks[queue].push(thread_task(std::move(task_fn), id, cost));
    created[queue]++;

    size_t i = queue;
    if (tasks[i].nr_threads_waiting == 0) {
        for (i = 0; i < tasks.size() && tasks[i].nr_threads_waiting == 0; i++) {}
    }
    if (i < tasks.size())
        tasks[i].not_empty.notify_one();
}

void thread_pool::drain_queue(size_t const queue, bool blocking)
{
    auto lock = get_lock();
    for (size_t const cr = created[queue]; joined[queue] < cr;) {
        if (!blocking) break;
        tasks[queue].task_done.wait(lock);
    }
}

thread_task thread_pool::get_task(size_t & preferred_queue)
{
    ASSERT(!is_synchronous());

    auto lock = get_lock();

    while (!kill_threads && all_task_queues_empty()) {
        /* No work -> tell this thread to wait until work becomes available.
           We also leave the loop when the thread needs to die.
           The while() protects against spurious wake-ups that can fire even if
           the queue is still empty. */
        tasks[preferred_queue].nr_threads_waiting++;
        tasks[preferred_queue].not_empty.wait(lock);
        tasks[preferred_queue].nr_threads_waiting--;
    }
    thread_task task(true);
    if (!(kill_threads && all_task_queues_empty())) {
        size_t & i(preferred_queue);
        if (tasks[i].empty()) {
            for (i = 0; i < tasks.size() && tasks[i].empty(); i++) {}
        }
        /* There must have been a non-empty queue or we'd still be in the
           while() loop above */
        ASSERT_ALWAYS(i < tasks.size());
        task = tasks[i].top();
        tasks[i].pop();
    }

    return task;
}

void thread_pool::drain_all_queues()
{
    for (size_t queue = 0; queue < created.size(); ++queue) {
        drain_queue(queue);
    }
}
