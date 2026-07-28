#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cstdint> // for SIZE_MAX

#include <condition_variable>
#include <mutex>
#include <queue>   // for queue, priority_queue
#include <vector>
#include <utility>
#include <thread>
#include <exception>
#include <functional>
#include <tuple>

#include "macros.h"
#include "threadpool.hpp"
#include "timing.h" // for seconds_thread
#include "utils_cxx.hpp"

/*
  With multiple queues, when new work is added to a queue, we need to be able
  to wake up one of the threads that prefer work from that queue. Thus we need
  multiple condition variables. If no threads that prefer work from that queue
  are currently waiting, we need to wake up some other thread.

  With k queues, we need k condition variables c[] and k semaphores s[].
  When a thread that prefers queue i waits for work, in increases s[i] and
  starts waiting on c[i]. When a thread that was waiting is woken up, it
  decreases s[i]. When work is added to queue j, it checks whether s[j] is
  non-zero:
    - if so, it signals c[j]
    - if not, it tests whether any other c[l] is non-zero
      - if so, it signals c[l]
      - if not, then no threads are currently sleeping, and nothing needs to be
  done

  We use a simple size_t variable as the semaphore; accesses are
  mutex-protected.
*/

worker_thread::worker_thread(thread_pool & _pool, size_t const _preferred_queue,
                             bool several_threads)
    : pool(_pool)
    , preferred_queue(several_threads ? _preferred_queue : SIZE_MAX)
{
    if (!several_threads) {
        // the "thread" member is uninitialized, but that is fine since
        // it's only used in the dtor anyway, under the condition that
        // preferred_queue != SIZE_MAX.
        // coverity[uninit_member]
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

class thread_task
{
  public:
    std::function<task_result*(worker_thread*)> func;
    int id = 0;
    double cost = 0.0; // costly tasks are scheduled first.

    bool is_terminal() const { return !func; }

    thread_task() = default;
    thread_task(std::function<task_result*(worker_thread*)> f, int id, double cost)
        : func(std::move(f))
        , id(id)
        , cost(cost)
    {}
    explicit thread_task(bool) {} // Terminal signal task

    task_result* operator()(worker_thread* w) const {
        if (func) return func(w);
        return nullptr;
    }

    bool operator<(thread_task const & y) const
    {
        return std::tie(cost, id) < std::tie(y.cost, y.id);
    }
};

class tasks_queue
    : public std::priority_queue<thread_task, std::vector<thread_task>>
    , private NonCopyable
{
  public:
    std::mutex mx;
    std::condition_variable not_empty;
    size_t nr_threads_waiting = 0;
    tasks_queue() = default;
};

class results_queue
    : public std::queue<task_result *>
    , private NonCopyable
{
  public:
    std::condition_variable not_empty;
};

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
    /* Threads start accessing the queues as soon as they run */
    threads.reserve(nr_threads);
    for (size_t i = 0; i < nr_threads; i++)
        threads.emplace_back(*this, 0, !is_synchronous());
};

// ok, if an exception is raised we'll die abruptly.
// coverity[exn_spec_violation]
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
    /* we removed the per-thread timer, because that goes in the way
     * of our intent to make threads more special-q agnostic: timers are
     * attached to the nfs_aux structure, now. This implies that all
     * routines that are called as worker threads must activate their timer
     * on entry.
     *
     */
    ASSERT_ALWAYS(!is_synchronous());
    double tt = -seconds_thread();
    while (true) {
        size_t queue = I.preferred_queue;
        thread_task task = get_task(queue);
        if (task.is_terminal())
            break;
        try {
            tt += seconds_thread();
            task_result * result = task(&I);
            tt -= seconds_thread();
            add_result(queue, result ? result : new task_result());
        } catch (...) {
            tt -= seconds_thread();
            add_exception(queue, std::current_exception());
            add_result(queue, new task_result()); // Keep task count synchronized
        }
    }
    tt += seconds_thread();
    std::lock_guard<std::mutex> const dummy(mm_cumulated_wait_time);
    cumulated_wait_time += tt;
}

bool thread_pool::all_task_queues_empty() const
{
    for (auto const & T: tasks)
        if (!T.empty())
            return false;
    return true;
}

void thread_pool::enqueue_task(std::function<task_result*(worker_thread*)> task_fn,
                               int const id, size_t const queue, double cost)
{
    if (is_synchronous()) {
        created[queue]++;
        try {
            task_result * result = task_fn(threads.data());
            results[queue].push(result ? result : new task_result());
        } catch (...) {
            add_exception(queue, std::current_exception());
            results[queue].push(new task_result());
        }
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

void thread_pool::add_task(task_function_t func, task_parameters * params,
                           int const id, size_t const queue, double cost)
{
    auto task_fn = [func, params, id](worker_thread* w) -> task_result* {
        return func(w, params, id);
    };
    enqueue_task(std::move(task_fn), id, queue, cost);
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
    if (kill_threads && all_task_queues_empty()) {
        /* then the default object above is appropriate for signaling all
         * workers so that they terminate.
         */
    } else {
        /* Find a non-empty task queue, starting with the preferred one */
        size_t & i(preferred_queue);
        if (tasks[i].empty()) {
            for (i = 0; i < tasks.size() && tasks[i].empty(); i++) {
            }
        }
        /* There must have been a non-empty queue or we'd still be in the
           while() loop above */
        ASSERT_ALWAYS(i < tasks.size());
        task = tasks[i].top();
        tasks[i].pop();
    }

    return task;
}

void thread_pool::add_result(size_t const queue, task_result * const result)
{
    ASSERT(!is_synchronous()); // synchronous case: see add_task
    ASSERT_ALWAYS(queue < results.size());
    auto lock = get_lock();
    results[queue].push(result);
    results[queue].not_empty.notify_one();
}

void thread_pool::add_exception(size_t const queue, std::exception_ptr e)
{
    ASSERT(!is_synchronous());
    ASSERT_ALWAYS(queue < results.size());
    auto lock = get_lock();
    exceptions[queue].push(std::move(e));
    results[queue].not_empty.notify_one();
}

std::exception_ptr thread_pool::get_exception(size_t const queue)
{
    ASSERT_ALWAYS(queue < exceptions.size());
    auto lock = get_lock();
    if (!exceptions[queue].empty()) {
        auto e = exceptions[queue].front();
        exceptions[queue].pop();
        return e;
    }
    return nullptr;
}

void thread_pool::rethrow_exceptions(size_t const queue)
{
    if (auto e = get_exception(queue))
        std::rethrow_exception(e);
}

/* Get a result from the specified results queue. If no result is available,
   waits with blocking=true, and returns nullptr with blocking=false. */
task_result * thread_pool::get_result(size_t const queue, bool const blocking)
{
    task_result * result;
    ASSERT_ALWAYS(queue < results.size());

    /* works both in synchronous and non-synchronous case */
    auto lock = get_lock();
    if (!blocking && results[queue].empty()) {
        result = nullptr;
    } else {
        while (results[queue].empty())
            results[queue].not_empty.wait(lock);
        result = results[queue].front();
        results[queue].pop();
        joined[queue]++;
    }
    return result;
}

void thread_pool::drain_queue(size_t const queue, void (*f)(task_result *))
{
    /* works both in synchronous and non-synchronous case */
    auto lock = get_lock();
    for (size_t const cr = created[queue]; joined[queue] < cr;) {
        while (results[queue].empty())
            results[queue].not_empty.wait(lock);
        task_result * result = results[queue].front();
        results[queue].pop();
        joined[queue]++;
        if (f)
            f(result);
        delete result;
    }
}

void thread_pool::drain_all_queues()
{
    for (size_t queue = 0; queue < results.size(); ++queue) {
        drain_queue(queue);
    }
}
