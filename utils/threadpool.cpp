#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cstdint>

#include <chrono>
#include <condition_variable>
#include <exception>
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
    , results(nr_queues)
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
    for (auto const & R: results)
        ASSERT_ALWAYS_NOTHROW(R.empty());
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
        task(&I);
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

void thread_pool::add_task(task_function_t func, task_parameters * params,
                           int const id, size_t const queue, double cost)
{
    enqueue_task([func, params, id](worker_thread* w) {
        task_result* res = func(w, params, id);
        delete res;
    }, id, queue, cost);
}

void thread_pool::drain_queue(size_t const queue, bool blocking)
{
    auto lock = get_lock();
    for (size_t const cr = created[queue]; joined[queue] < cr;) {
        if (results[queue].empty()) {
            if (!blocking) break;
            tasks[queue].not_empty.wait(lock);
            continue;
        }

        if (!blocking && results[queue].front().wait_for(std::chrono::seconds(0)) != std::future_status::ready)
            break;


        auto fut = std::move(results[queue].front());
        results[queue].pop();
        joined[queue]++;

        lock.unlock();
        try {
            fut.get();
        } catch (...) {
            std::scoped_lock guard(pool_mutex);
            exceptions[queue].push(std::current_exception());
        }
        lock.lock();
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
    for (size_t queue = 0; queue < results.size(); ++queue) {
        drain_queue(queue);
    }
}
