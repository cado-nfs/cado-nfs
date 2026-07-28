#ifndef CADO_THREADPOOL_HPP
#define CADO_THREADPOOL_HPP

#include <cstddef>

#include <condition_variable>
#include <exception>
#include <future>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "utils_cxx.hpp"
#include "macros.h"

class task_parameters {
public:
    virtual ~task_parameters() = default;
};

class task_result {
public:
    virtual ~task_result() = default;
};

class worker_thread;
class thread_pool;

class thread_task
{
  public:
    std::function<void(worker_thread*)> func;
    int id = 0;
    double cost = 0.0; // costly tasks are scheduled first.

    bool is_terminal() const { return !func; }

    thread_task() = default;
    thread_task(std::function<void(worker_thread*)> f, int id, double cost)
        : func(std::move(f))
        , id(id)
        , cost(cost)
    {}
    explicit thread_task(bool) {}

    void operator()(worker_thread* w) const {
        if (func) func(w);
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
    std::condition_variable task_done;
    size_t nr_threads_waiting = 0;
    tasks_queue() = default;
};

class worker_thread {
    friend class thread_pool;
    thread_pool &pool;
    std::jthread thread;
    const size_t preferred_queue;

public:
    int rank() const;
    int nthreads() const;
    thread_pool & get_pool() { return pool; }
    worker_thread(thread_pool &, size_t, bool = true);
    bool is_synchronous() const;
};

using task_function_t = task_result *(*)(worker_thread * worker, task_parameters *, int id);

class thread_pool : private NonCopyable {
public:
    static constexpr int QUEUE_GENERIC = 0;
    static constexpr int QUEUE_ECM = 1;
    static constexpr int QUEUE_MISC = 2;
    static constexpr int NQUEUES = 3;

private:
    friend class worker_thread;

    bool sync_mode = false;
    mutable std::mutex pool_mutex;

    std::unique_lock<std::mutex> get_lock() const {
        std::unique_lock lock(pool_mutex, std::defer_lock);
        if (!sync_mode)
            lock.lock();
        return lock;
    }

    std::vector<worker_thread> threads;
    std::vector<tasks_queue> tasks;

    // Type-erased completion handles (std::future<void>) for drain_queue tracking
    std::vector<std::queue<std::exception_ptr>> exceptions;

    std::vector<size_t> created;
    std::vector<size_t> joined;

    bool kill_threads;
    double & store_wait_time;

    void thread_work_on_tasks(worker_thread &);
    thread_task get_task(size_t& queue);
    bool all_task_queues_empty() const;

    // Universal enqueuer supporting any return type R
    void enqueue_task(std::function<void(worker_thread*)> task_fn, int id, size_t queue, double cost);

public:
    bool is_synchronous() const { return sync_mode; }
    double cumulated_wait_time = 0;
    std::mutex mm_cumulated_wait_time;

    size_t size() {
        auto lock = get_lock();
        return threads.size();
    }

    thread_pool(size_t _nr_threads, double & store_wait_time,
                size_t nr_queues = 1, bool sync_thread_pool = false);
    ~thread_pool();

    void drain_queue(size_t queue, bool blocking = true);
    void drain_all_queues();

    template<typename T>
    std::vector<T> get_exceptions(size_t const queue = 0) {
        auto lock = get_lock();
        std::vector<T> res;
        std::queue<std::exception_ptr> remaining;

        while (!exceptions[queue].empty()) {
            auto ex_ptr = std::move(exceptions[queue].front());
            exceptions[queue].pop();

            try {
                if (ex_ptr)
                    std::rethrow_exception(ex_ptr);
            } catch (const T& e) {
                res.push_back(e);
            } catch (...) {
                remaining.push(std::move(ex_ptr));
            }
        }

        exceptions[queue] = std::move(remaining);
        return res;
    }

    // Unified add_task_lambda: supports any return type R, and signatures (worker*, int) or (worker*)
    template<typename T>
    auto add_task_lambda(T && f, const int id = 0, const size_t queue = 0, double cost = 0.0)
    {
        std::function<void(worker_thread*)> task_fn;

        if constexpr (std::is_invocable_v<T, worker_thread*, int>) {
            task_fn = [f = std::forward<T>(f), id](worker_thread* w) {
                f(w, id);
            };
        } else {
            task_fn = [f = std::forward<T>(f)](worker_thread* w) {
                f(w);
            };
        }
        enqueue_task(task_fn, id, queue, cost);
    }

    // add_shared_task
    template<typename T>
    using shared_task = std::shared_ptr<T>;

    template<typename T, typename... Args>
    static shared_task<T> make_shared_task(Args&&... args) {
        return std::make_shared<T>(std::forward<Args>(args)...);
    }

    template<typename T>
    auto add_shared_task(shared_task<T> const & task_ptr, const int id = 0, const size_t queue = 0, double cost = 0.0)
    {
        return add_task_lambda([task_ptr](worker_thread* w, int task_id) {
            return (*task_ptr)(w, task_id);
        }, id, queue, cost);
    }

    // Legacy function pointer overload
    void add_task(task_function_t func, task_parameters * params, int id, size_t queue = 0, double cost = 0.0);

    template <typename F>
    auto add_future_task(size_t queue, double cost, F&& task_fn)
    {
        using R = std::invoke_result_t<F, worker_thread*>;

        auto pkg = std::make_shared<std::packaged_task<R(worker_thread*)>>(
                std::forward<F>(task_fn)
                );
        std::future<R> fut = pkg->get_future();

        // Enqueue as a fire-and-forget task that invokes the packaged task
        add_task_lambda([pkg](worker_thread* w) {
                (*pkg)(w);
                }, 0, queue, cost);

        return fut;
    }

};

#endif  /* CADO_THREADPOOL_HPP */
