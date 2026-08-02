#ifndef CADO_THREADPOOL_HPP
#define CADO_THREADPOOL_HPP

#include <cstddef>

#include <atomic>
#include <condition_variable>
#include <deque>
#include <exception>
#include <functional>
#include <future>
#include <map>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "chronograms.hpp"
#include "utils_cxx.hpp"

class worker_thread;
class thread_pool;
class task_group;

class thread_task
{
  public:
    std::function<void(worker_thread*)> func;
    int id = 0;
    double cost = 0.0; // costly tasks are scheduled first.
    task_group * group = nullptr;

    bool is_terminal() const { return !func; }

    thread_task() = default;
    thread_task(std::function<void(worker_thread*)> f, int id, double cost, task_group * tg = nullptr)
        : func(std::move(f))
        , id(id)
        , cost(cost)
        , group(tg)
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

class task_group : private NonCopyable {
    friend class thread_pool;
    mutable std::mutex mx;
    std::condition_variable cv;

    alignas(64) std::atomic<size_t> created_count{0};
    alignas(64) std::atomic<size_t> finished_count{0};
    alignas(64) std::atomic<size_t> joined_count{0};

    std::function<void()> completion_cb;

    void notify_joined() {
        // 1. Mark task as finished
        size_t const f = finished_count.fetch_add(1, std::memory_order_acq_rel) + 1;
        size_t const c = created_count.load(std::memory_order_acquire);

        // 2. If this is the last task, extract the callback
        std::function<void()> cb;
        if (f >= c) {
            const std::scoped_lock lock(mx);
            // Re-verify under lock in case new tasks were added
            if (finished_count.load(std::memory_order_relaxed) >=
                created_count.load(std::memory_order_relaxed))
            {
                cb = std::move(completion_cb);
                completion_cb = nullptr;
            }
        }

        // 3. Execute the callback completely lock-free
        if (cb) {
            cb();
        }

        // 4. Mark task as joined (and callback as completed)
        size_t const j = joined_count.fetch_add(1, std::memory_order_release) + 1;
        if (j >= created_count.load(std::memory_order_acquire)) {
            // Wake up wait(). The lock ensures we don't miss a wakeup if wait()
            // is just about to sleep.
            const std::scoped_lock lock(mx);
            cv.notify_all();
        }
    }

public:
    task_group() = default;

    void wait() {
        // std::condition_variable fundamentally requires a lock to wait safely.
        // This is only held while blocking, never during task execution.
        std::unique_lock lock(mx);
        cv.wait(lock, [this]() {
            return joined_count.load(std::memory_order_acquire) >=
                   created_count.load(std::memory_order_acquire);
        });
    }

    void on_complete(std::function<void()> cb) {
        {
            const std::scoped_lock lock(mx);
            if (finished_count.load(std::memory_order_relaxed) >=
                created_count.load(std::memory_order_relaxed))
            {
                // Group is already done. Fake a task to prevent wait()
                // from returning, so that we can run the callback
                // lock-free.
                created_count.fetch_add(1, std::memory_order_relaxed);
            } else {
                completion_cb = std::move(cb);
                return;
            }
        }

        // Execute late callback exactly as if it were a worker task
        cb();
        finished_count.fetch_add(1, std::memory_order_acq_rel);
        size_t const j = joined_count.fetch_add(1, std::memory_order_release) + 1;
        if (j >= created_count.load(std::memory_order_acquire)) {
            const std::scoped_lock lock(mx);
            cv.notify_all();
        }
    }

    size_t created() const {
        return created_count.load(std::memory_order_relaxed);
    }

    size_t joined() const {
        return joined_count.load(std::memory_order_relaxed);
    }
};

class tasks_queue : private NonCopyable {
public:
    std::mutex mx;
    std::condition_variable task_done;
    tasks_queue() = default;
};

// Per-worker, per-queue local task container
struct worker_queue {
    alignas(64) std::mutex mx;
    std::deque<thread_task> tasks;
    std::atomic<size_t> count;  /* use this to provide lock-free info about emptiness */
};

class worker_thread {
    friend class thread_pool;
    thread_pool &pool;
    std::jthread thread;
    const size_t preferred_queue;
    std::vector<chronograms::bubble> chronogram;

public:
    auto trace(chronograms::bubble_info b) {
        return chronograms::bubble_guard(chronogram, b);
    }
    static auto trace(worker_thread * wrk, chronograms::bubble_info b)
    {
        return wrk ? wrk->trace(b) : chronograms::bubble_guard::dummy();
    }
    int rank() const;
    int nthreads() const;
    thread_pool & get_pool() { return pool; }
    worker_thread(thread_pool &, size_t, bool = true);
    bool is_synchronous() const;
};

template <typename F, typename... Args>
concept is_task_callable = std::is_invocable_v<F, worker_thread*, Args...> || std::is_invocable_v<F, Args...>;

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

    std::vector<worker_thread> threads;
    std::vector<tasks_queue> tasks;

    // Distributed per-worker, per-queue work-stealing deques
    // worker_queues[worker_id][queue_id]
    std::vector<std::vector<worker_queue>> worker_queues;
    
    std::vector<std::mutex> exceptions_mutex;
    std::vector<std::queue<std::exception_ptr>> exceptions;

    // Lock-free atomic counters aligned to 64 bytes
    // Cacheline-aligned atomic counter to eliminate false sharing across CPU sockets
    // note that in c++20, std::atomic<> is value-initialized.
    struct alignas(64) padded_atomic : std::atomic<size_t> { };

    std::vector<padded_atomic> created;
    std::vector<padded_atomic> joined;

    std::condition_variable work_cv;
    std::atomic<size_t> nr_threads_waiting{0};
    std::atomic<size_t> round_robin_counter{0};

    bool kill_threads{false};
    double & store_wait_time;

    void thread_work_on_tasks(worker_thread &);
    bool pop_local_or_steal(size_t worker_id, size_t preferred_queue, thread_task & task);
    bool all_task_queues_empty() const;

    void enqueue_task(std::function<void(worker_thread*)> task_fn, size_t queue, double cost, task_group * tg = nullptr);

public:
    auto trace_on_leader(chronograms::bubble_info b)
    {
        return threads.front().trace(b);
    }
    /* use this when you want to record on a _pointer_, without knowing
     * if the pointer can be dereferenced.
     */
    static auto trace_on_leader(thread_pool * pool, chronograms::bubble_info b)
    {
        return pool ? pool->trace_on_leader(b) : chronograms::bubble_guard::dummy();
    }
    void collect_traces(std::map<size_t, std::vector<chronograms::bubble>> & destination, size_t thread_index_offset);

    bool is_synchronous() const { return sync_mode; }
    double cumulated_wait_time = 0;
    std::mutex mm_cumulated_wait_time;

    size_t size() const { return threads.size(); }

    thread_pool(size_t _nr_threads, double & store_wait_time,
                size_t nr_queues = 1, bool sync_thread_pool = false);
    ~thread_pool();

    void drain_queue(size_t queue, bool blocking = true);
    void drain_all_queues();

    template<typename T>
    std::vector<T> get_exceptions(size_t const queue = 0) {
        const std::scoped_lock lock(exceptions_mutex[queue]);
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

    /* {{{ Variadic add_task overloads */
    template <typename F, typename... Args>
    requires is_task_callable<F, Args...>
    void add_task(task_group & tg, size_t queue, double cost, F&& f, Args&&... args)
    {
        auto task_fn = [f = std::forward<F>(f), ...args = std::forward<Args>(args)](worker_thread* w) mutable {
            if constexpr (std::is_invocable_v<F, worker_thread*, Args...>) {
                std::invoke(f, w, std::forward<Args>(args)...);
            } else {
                std::invoke(f, std::forward<Args>(args)...);
            }
        };

        enqueue_task(std::move(task_fn), queue, cost, &tg);
    }

    template <typename F, typename... Args>
    requires is_task_callable<F, Args...>
    void add_task(task_group & tg, size_t queue, F&& f, Args&&... args)
    {
        add_task(tg, queue, 0.0, std::forward<F>(f), std::forward<Args>(args)...);
    }

    template <typename F, typename... Args>
    requires is_task_callable<F, Args...>
    void add_task(task_group & tg, F&& f, Args&&... args)
    {
        add_task(tg, QUEUE_GENERIC, 0.0, std::forward<F>(f), std::forward<Args>(args)...);
    }

    template <typename F, typename... Args>
    requires is_task_callable<F, Args...>
    void add_task(size_t queue, double cost, F&& f, Args&&... args)
    {
        auto task_fn = [f = std::forward<F>(f), ...args = std::forward<Args>(args)](worker_thread* w) mutable {
            if constexpr (std::is_invocable_v<F, worker_thread*, Args...>) {
                std::invoke(f, w, std::forward<Args>(args)...);
            } else {
                std::invoke(f, std::forward<Args>(args)...);
            }
        };

        enqueue_task(std::move(task_fn), queue, cost, nullptr);
    }

    template <typename F, typename... Args>
    requires is_task_callable<F, Args...>
    void add_task(size_t queue, F&& f, Args&&... args)
    {
        add_task(queue, 0.0, std::forward<F>(f), std::forward<Args>(args)...);
    }

    template <typename F, typename... Args>
    requires is_task_callable<F, Args...>
    void add_task(F&& f, Args&&... args)
    {
        add_task(QUEUE_GENERIC, 0.0, std::forward<F>(f), std::forward<Args>(args)...);
    }
    /* }}} */

    /* {{{ add_future_task for value-returning producer tasks */
    template <typename F, typename... Args>
    requires is_task_callable<F, Args...>
    auto add_future_task(size_t queue, double cost, F&& f, Args&&... args)
    {
        if constexpr (std::is_invocable_v<F, worker_thread*, Args...>) {
            using R = std::invoke_result_t<F, worker_thread*, Args...>;
            auto pkg = std::make_shared<std::packaged_task<R(worker_thread*, Args...)>>(
                std::forward<F>(f)
            );
            std::future<R> fut = pkg->get_future();
            add_task(queue, cost, [pkg](worker_thread* w, auto&&... a) {
                (*pkg)(w, std::forward<decltype(a)>(a)...);
            }, std::forward<Args>(args)...);
            return fut;
        } else {
            using R = std::invoke_result_t<F, Args...>;
            auto pkg = std::make_shared<std::packaged_task<R(Args...)>>(
                std::forward<F>(f)
            );
            std::future<R> fut = pkg->get_future();
            add_task(queue, cost, [pkg](auto&&... a) {
                (*pkg)(std::forward<decltype(a)>(a)...);
            }, std::forward<Args>(args)...);
            return fut;
        }
    }

    template <typename F, typename... Args>
    requires is_task_callable<F, Args...>
    auto add_future_task(F&& f, Args&&... args)
    {
        return add_future_task(QUEUE_GENERIC, 0.0, std::forward<F>(f), std::forward<Args>(args)...);
    }
    /* }}} */
};

#endif  /* CADO_THREADPOOL_HPP */
