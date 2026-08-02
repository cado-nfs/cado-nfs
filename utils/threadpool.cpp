#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cstdint>

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <exception>
#include <functional>
#include <map>
#include <mutex>
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
    , exceptions_mutex(nr_queues)
    , exceptions(nr_queues)
    , created(nr_queues)
    , joined(nr_queues)
    , store_wait_time(store_wait_time)
{
    ASSERT_ALWAYS(nr_threads == 1 || !sync_thread_pool);
    sync_mode = sync_thread_pool;

    worker_queues.resize(nr_threads);
    for (size_t i = 0; i < nr_threads; ++i) {
        worker_queues[i] = std::vector<worker_queue>(nr_queues);
    }

    threads.reserve(nr_threads);
    for (size_t i = 0; i < nr_threads; i++)
        threads.emplace_back(*this, 0, !is_synchronous());
};

thread_pool::~thread_pool()
{
    drain_all_queues();
    {
        const std::scoped_lock lock(pool_mutex);
        kill_threads = true;
    }
    work_cv.notify_all();

    threads.clear();

    for (auto const & E: exceptions)
        ASSERT_ALWAYS_NOTHROW(E.empty());
    store_wait_time += cumulated_wait_time;
}

bool thread_pool::all_task_queues_empty() const
{
    for (size_t q = 0; q < tasks.size(); ++q) {
        if (joined[q].load(std::memory_order_relaxed) < created[q].load(std::memory_order_relaxed))
            return false;
    }
    return true;
}

bool thread_pool::pop_local_or_steal(size_t worker_id, size_t preferred_queue, thread_task & task)
{
    size_t const nq = tasks.size();
    size_t const nw = threads.size();

    // 1. Try owner's local queue for preferred_queue (LIFO pop for local
    // cache warmth)
    {
        auto & wq = worker_queues[worker_id][preferred_queue];
        // Use lock-free atomic load for the fast-path check
        if (wq.count.load(std::memory_order_acquire) > 0) {
            const std::scoped_lock lock(wq.mx);
            if (!wq.tasks.empty()) {
                task = std::move(wq.tasks.back());
                wq.tasks.pop_back();
                wq.count.store(wq.tasks.size(), std::memory_order_release);
                return true;
            }
        }
}

    // 2. Try owner's local queue for other queues
    for (size_t q = 0; q < nq; ++q) {
        if (q == preferred_queue) continue;
        auto & wq = worker_queues[worker_id][q];
        if (wq.count.load(std::memory_order_acquire) > 0) {
            const std::scoped_lock lock(wq.mx);
            if (!wq.tasks.empty()) {
                task = std::move(wq.tasks.back());
                wq.tasks.pop_back();
                wq.count.store(wq.tasks.size(), std::memory_order_release);
                return true;
            }
        }
    }

    // 3. Steal from other workers (FIFO steal from front)
    for (size_t offset = 1; offset < nw; ++offset) {
        size_t const victim = (worker_id + offset) % nw;
        for (size_t q = 0; q < nq; ++q) {
            size_t const actual_q = (preferred_queue + q) % nq;
            auto & wq = worker_queues[victim][actual_q];
            if (wq.count.load(std::memory_order_acquire) > 0) {
                const std::scoped_lock lock(wq.mx);
                if (!wq.tasks.empty()) {
                    task = std::move(wq.tasks.front());
                    wq.tasks.pop_front();
                    wq.count.store(wq.tasks.size(), std::memory_order_release);
                    return true;
                }
            }
        }
    }

    return false;
}

void thread_pool::thread_work_on_tasks(worker_thread & I)
{
    ASSERT_ALWAYS(!is_synchronous());
    auto const worker_id = static_cast<size_t>(I.rank());
    double tt = -seconds_thread();

    while (true) {
        size_t const queue = I.preferred_queue;
        thread_task task;

        if (!pop_local_or_steal(worker_id, queue, task)) {
            // No work found -> wait on condition variable
           std::unique_lock<std::mutex> lock(pool_mutex);
           nr_threads_waiting.fetch_add(1, std::memory_order_seq_cst);

            work_cv.wait(lock, [this, worker_id, queue, &task]() {
                return kill_threads || pop_local_or_steal(worker_id, queue, task);
            });

           nr_threads_waiting.fetch_sub(1, std::memory_order_seq_cst);

            if (kill_threads && task.is_terminal() && all_task_queues_empty()) {
                break;
            }
        }

        if (task.is_terminal()) {
            if (kill_threads) break;
            continue;
        }

        tt += seconds_thread();
        try {
            task(&I);
        } catch (...) {
            const std::scoped_lock guard(exceptions_mutex[task.id % tasks.size()]);
            exceptions[task.id % tasks.size()].push(std::current_exception());
        }

        if (task.group) {
            task.group->notify_joined();
        }

        size_t const queue_id = (task.id < 0) ? 0 : static_cast<size_t>(task.id) % tasks.size();
        size_t const current_joined = joined[queue_id].fetch_add(1, std::memory_order_acq_rel) + 1;
        size_t const target_created = created[queue_id].load(std::memory_order_acquire);

        if (current_joined >= target_created) {
            const std::scoped_lock lock(tasks[queue_id].mx);
            tasks[queue_id].task_done.notify_all();
        }

        tt -= seconds_thread();
    }

    tt += seconds_thread();
    const std::scoped_lock dummy(mm_cumulated_wait_time);
    cumulated_wait_time += tt;
}

void thread_pool::enqueue_task(std::function<void(worker_thread*)> task_fn, size_t queue, double cost, task_group * tg)
{
    if (tg) {
        tg->created_count.fetch_add(1, std::memory_order_relaxed);
    }

    created[queue].fetch_add(1, std::memory_order_relaxed);

    if (is_synchronous()) {
        try {
            worker_thread synchronous_worker(*this, queue, false);
            task_fn(&synchronous_worker);
        } catch (...) {
            const std::scoped_lock guard(exceptions_mutex[queue]);
            exceptions[queue].push(std::current_exception());
        }
        joined[queue].fetch_add(1, std::memory_order_release);
        if (tg) {
            tg->notify_joined();
        }
        return;
    }

    thread_task t(std::move(task_fn), static_cast<int>(queue), cost, tg);

    size_t const victim = round_robin_counter.fetch_add(1, std::memory_order_relaxed) % threads.size();
    auto & wq = worker_queues[victim][queue];
    {
        const std::scoped_lock lock(wq.mx);
        if (!wq.tasks.empty()) {
            auto it = std::ranges::lower_bound(wq.tasks, t,
                [](thread_task const & a, thread_task const & b) {
                    return a.cost < b.cost;
                });
            wq.tasks.insert(it, std::move(t));
        } else {
            wq.tasks.push_back(std::move(t));
        }
        // Update atomic counter under lock
        wq.count.store(wq.tasks.size(), std::memory_order_release);
    }

    std::atomic_thread_fence(std::memory_order_seq_cst);

    if (nr_threads_waiting.load(std::memory_order_seq_cst) > 0) {
        const std::scoped_lock lock(pool_mutex);
        work_cv.notify_one();
    }
}

void thread_pool::drain_queue(size_t const queue, bool blocking)
{
    if (!blocking) return;

    size_t const target = created[queue].load(std::memory_order_acquire);
    if (joined[queue].load(std::memory_order_acquire) >= target)
        return;

    std::unique_lock<std::mutex> lock(tasks[queue].mx);
    tasks[queue].task_done.wait(lock, [this, queue]() {
        return joined[queue].load(std::memory_order_acquire) >=
               created[queue].load(std::memory_order_acquire);
    });
}

void thread_pool::drain_all_queues()
{
    for (size_t queue = 0; queue < tasks.size(); ++queue) {
        drain_queue(queue);
    }
}
void thread_pool::collect_traces(
        std::map<size_t, std::vector<chronograms::bubble>> & destination,
        size_t thread_index_offset)
{
    for(size_t i = 0 ; i < threads.size() ; i++) {
        const size_t j = thread_index_offset + i;
        ASSERT_ALWAYS(!destination.contains(j));
        destination[j] = std::move(threads[i].chronogram);
    }
}
