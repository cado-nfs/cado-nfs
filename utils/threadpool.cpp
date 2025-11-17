#include "cado.h" // IWYU pragma: keep

#include <cstddef>
#include <cstdint> // for SIZE_MAX

#include <condition_variable>
#include <mutex>
#include <queue>   // for queue, priority_queue
#include <vector>
#include <pthread.h>

#include "clonable-exception.hpp" // for clonable_exception
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
        // the "pthread" member is uninitialized, but that is fine since
        // it's only used in the dtor anyway, under the condition that
        // preferred_queue != SIZE_MAX.
        // coverity[uninit_member]
        return;
    }
    int const rc = pthread_create(
        &thread, nullptr, thread_pool::thread_work_on_tasks_static, this);
    ASSERT_ALWAYS(rc == 0);
}

worker_thread::~worker_thread()
{
    if (preferred_queue == SIZE_MAX)
        return;
    int const rc = pthread_join(thread, nullptr);
    // timer.self will be essentially lost here. the best thing to do is to
    // create a phony task which collects the busy wait times for all
    // threads, at regular intervals, so that timer.self will be
    // insignificant.
    // pool.timer += timer;
    ASSERT_ALWAYS_NOTHROW(rc == 0);
}

size_t worker_thread::rank() const
{
    return this - pool.threads.data();
}
size_t worker_thread::nthreads() const
{
    return pool.threads.size();
}
bool worker_thread::is_synchronous() const
{
    return pool.is_synchronous();
}

class thread_task
{
  public:
    task_function_t func = nullptr;
    task_parameters * parameters = nullptr;
    int id = 0;
    double cost = 0.0; // costly tasks are scheduled first.

    bool is_terminal() const { return func == nullptr; }

    thread_task(task_function_t _func, task_parameters * _parameters, int _id,
                double _cost)
        : func(_func)
        , parameters(_parameters)
        , id(_id)
        , cost(_cost)
    {
    }
    thread_task(bool) {}
    task_result * operator()(worker_thread * w) const
    {
        return (*func)(w, parameters, id);
    }
};

struct thread_task_cmp {
    bool operator()(thread_task const & x, thread_task const & y) const
    {
        if (x.cost < y.cost)
            return true;
        if (x.cost > y.cost)
            return false;
        // if costs are equal, compare ids (they should be distinct)
        return x.id < y.id;
    }
};

class tasks_queue
    : public std::priority_queue<thread_task, std::vector<thread_task>,
                                 thread_task_cmp>
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

class exceptions_queue
    : public std::queue<clonable_exception *>
    , private NonCopyable
{
  public:
    std::condition_variable not_empty;
};

thread_pool::thread_pool(size_t const nr_threads, double & store_wait_time,
                         size_t const nr_queues, bool sync_thread_pool)
    : monitor_or_synchronous(sync_thread_pool)
    , tasks(nr_queues)
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
        my_unique_lock const u(*this);
        kill_threads = true;
        for (auto & T: tasks)
            broadcast(T.not_empty); /* Wakey wakey, time to die */
    }
    drain_all_queues();
    threads.clear(); /* does pthread_join */
    for (auto const & T: tasks)
        ASSERT_ALWAYS_NOTHROW(T.empty());
    for (auto const & R: results)
        ASSERT_ALWAYS_NOTHROW(R.empty());
    for (auto const & E: exceptions)
        ASSERT_ALWAYS_NOTHROW(E.empty());
    store_wait_time += cumulated_wait_time;
}

void * thread_pool::thread_work_on_tasks_static(void * arg)
{
    auto * I = static_cast<worker_thread *>(arg);
    I->pool.thread_work_on_tasks(*I);
    return nullptr;
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
            if (result != nullptr)
                add_result(queue, result);
        } catch (clonable_exception const & e) {
            tt -= seconds_thread();
            add_exception(queue, e.clone());
            /* We need to wake the listener... */
            add_result(queue, nullptr);
        }
    }
    tt += seconds_thread();
    /* tt is now the wall-clock time spent really within this function,
     * waiting for mutexes and condition variables...  */
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

void thread_pool::add_task(task_function_t func, task_parameters * params,
                           int const id, size_t const queue, double cost)
{
    if (is_synchronous()) {
        /* Execute the function right away, simulate the action of a
         * secondary thread fetching it from the task queue */
        created[queue]++;
        try {
            task_result * result = func(threads.data(), params, id);
            if (result != nullptr)
                results[queue].push(result);
        } catch (clonable_exception const & e) {
            exceptions[queue].push(e.clone());
            /* We do this in the asynchronous case. It isn't clear that
             * we need to do the same in the syncronous case. */
            results[queue].push(nullptr);
        }
        return;
    }
    ASSERT_ALWAYS(queue < tasks.size());

    my_unique_lock const u(*this);

    ASSERT_ALWAYS(!kill_threads);
    tasks[queue].push(thread_task(func, params, id, cost));
    created[queue]++;

    /* Find a queue with waiting threads, starting with "queue" */
    size_t i = queue;
    if (tasks[i].nr_threads_waiting == 0) {
        for (i = 0; i < tasks.size() && tasks[i].nr_threads_waiting == 0; i++) {
        }
    }
    /* If any queue with waiting threads was found, wake up one of them */
    if (i < tasks.size())
        signal(tasks[i].not_empty);
}

thread_task thread_pool::get_task(size_t & preferred_queue)
{
    ASSERT(!is_synchronous());

    my_unique_lock u(*this);

    while (!kill_threads && all_task_queues_empty()) {
        /* No work -> tell this thread to wait until work becomes available.
           We also leave the loop when the thread needs to die.
           The while() protects against spurious wake-ups that can fire even if
           the queue is still empty. */
        tasks[preferred_queue].nr_threads_waiting++;
        wait(tasks[preferred_queue].not_empty, u);
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
    my_unique_lock const u(*this);
    results[queue].push(result);
    signal(results[queue].not_empty);
}

void thread_pool::add_exception(size_t const queue, clonable_exception * e)
{
    ASSERT(!is_synchronous()); // synchronous case: see add_task
    ASSERT_ALWAYS(queue < results.size());
    my_unique_lock const u(*this);
    exceptions[queue].push(e);
    // do we use it ?
    signal(results[queue].not_empty);
}

/* Get a result from the specified results queue. If no result is available,
   waits with blocking=true, and returns nullptr with blocking=false. */
task_result * thread_pool::get_result(size_t const queue, bool const blocking)
{
    task_result * result;
    ASSERT_ALWAYS(queue < results.size());

    /* works both in synchronous and non-synchronous case */
    my_unique_lock u(*this);
    if (!blocking and results[queue].empty()) {
        result = nullptr;
    } else {
        while (results[queue].empty())
            wait(results[queue].not_empty, u);
        result = results[queue].front();
        results[queue].pop();
        joined[queue]++;
    }
    return result;
}

void thread_pool::drain_queue(size_t const queue, void (*f)(task_result *))
{
    /* works both in synchronous and non-synchronous case */
    my_unique_lock u(*this);
    for (size_t const cr = created[queue]; joined[queue] < cr;) {
        while (results[queue].empty())
            wait(results[queue].not_empty, u);
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

/* get an exception from the specified exceptions queue. This is
 * obviously non-blocking, because exceptions are exceptional. So when no
 * exception is there, we return nullptr. When there is one, we return a
 * pointer to a newly allocated copy of it.
 */
clonable_exception * thread_pool::get_exception(size_t const queue)
{
    clonable_exception * e = nullptr;
    ASSERT_ALWAYS(queue < exceptions.size());
    /* works both in synchronous and non-synchronous case */
    my_unique_lock const u(*this);
    if (!exceptions[queue].empty()) {
        e = exceptions[queue].front();
        exceptions[queue].pop();
    }
    return e;
}
