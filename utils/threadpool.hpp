#ifndef CADO_THREADPOOL_HPP
#define CADO_THREADPOOL_HPP

#include <cstddef>

#include <memory>
#include <mutex>
#include <type_traits>
#include <vector>
#include <utility>
#include <functional>

#include <pthread.h>

#include "utils_cxx.hpp"

struct clonable_exception;

/* C++11 already has classes for mutex and condition_variable */
/* All the synchronization stuff could be moved to the implementation if
   thread_pool used monitor as a dynamically allocated object. Tempting. */

/* Base for classes that hold parameters for worker functions */
class task_parameters {
  public:
  virtual ~task_parameters() = default;
};

/* Base for classes that hold results produced by worker functions */
class task_result {
  public:
  virtual ~task_result() = default;
};

class thread_task;
class tasks_queue;
class results_queue;
class exceptions_queue;
class thread_pool;


class worker_thread {
  friend class thread_pool;
  thread_pool &pool;
  pthread_t thread;
  const size_t preferred_queue;
public:
  worker_thread(worker_thread const &) = delete;
  worker_thread& operator=(worker_thread const &) = delete;

  // move is ok
  worker_thread(worker_thread&&) = default;
  // worker_thread& operator=(worker_thread&&) = default;
  int rank() const;
  int nthreads() const;
  /* It doesn't seem that unholy to me to have a thread access the pool
   * it originates from. It's possibly a good way to do continuations,
   * for example.
   */
  thread_pool & get_pool() { return pool; }
  worker_thread(thread_pool &, size_t, bool = true);
  ~worker_thread();
  bool is_synchronous() const;
};

typedef task_result *(*task_function_t)(worker_thread * worker, task_parameters *, int id);

class thread_pool : private NonCopyable {
    public:
        /* queue 0: main
         * queue 1: ECM
         * queue 2: things that we join almost immediately, but are
         * multithreaded nevertheless: alloc buckets, ...
         */
        static constexpr int QUEUE_GENERIC = 0;
        static constexpr int QUEUE_ECM = 1;
        static constexpr int QUEUE_MISC = 2;
        static constexpr int NQUEUES = 3;

    private:
        friend class worker_thread;

        bool sync_mode = false;
        mutable std::mutex pool_mutex;

        // Helper returning a locked unique_lock in async mode, or an unlocked one in sync mode
        std::unique_lock<std::mutex> get_lock() const {
            std::unique_lock<std::mutex> lock(pool_mutex, std::defer_lock);
            if (!sync_mode)
                lock.lock();
            return lock;
        }

        std::vector<worker_thread> threads;
        std::vector<tasks_queue> tasks;
        std::vector<results_queue> results;
        std::vector<exceptions_queue> exceptions;
        std::vector<size_t> created;
        std::vector<size_t> joined;

        bool kill_threads; /* If true, hands out kill tasks once work queues are empty */
        double & store_wait_time;

        static void * thread_work_on_tasks_static(void *worker);
        void thread_work_on_tasks(worker_thread &);
        thread_task get_task(size_t& queue);
        void add_result(size_t queue, task_result *result);
        void add_exception(size_t queue, clonable_exception * e);
        bool all_task_queues_empty() const;

        void enqueue_task(std::function<task_result*(worker_thread*)> task_fn, int id, size_t queue, double cost);

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
        task_result *get_result(size_t queue = 0, bool blocking = true);
        void drain_queue(size_t queue, void (*f)(task_result*) = NULL);
        void drain_all_queues();
        clonable_exception * get_exception(size_t queue = 0);
        template<typename T>
            T * get_exception(const size_t queue = 0) {
                return dynamic_cast<T*>(get_exception(queue));
            }
        template<typename T>
            std::vector<T> get_exceptions(const size_t queue = 0) {
                std::vector<T> res;
                for(T * e ; (e = get_exception<T>(queue)) != NULL; ) {
                    res.push_back(*e);
                    delete e;
                }
                return res;
            }

        /* {{{ add_task is the simplest interface. It does not even specify who has
         * ownership of the params object. Two common cases can be envisioned.
         *  - either the caller retains ownership, in which case it obviously
         *    has to join all threads before deletion.
         *  - or ownership is transferred to the callee, in which case it is
         *    obviously not shared: we have one params object per task spawned
         *    (possibly at some cost), even if all params objects are distinct.
         * 
         * In the latter case, the id field is only of limited use.
         */
        void add_task(task_function_t func, task_parameters * params, int id, size_t queue = 0, double cost = 0.0);
        /* }}} */

        /* {{{ add_task_lambda.
         *
         * This adds a task to process exactly one lambda function. The lambda
         * function is expected to take the worker thread as only argument.
         * The lambda
         * object is copied. As usual, any references held by the lambda at the
         * time of capture must still be alive at the time of execution, or
         * chaos ensues. This must be guaranteed by the caller.
         *
         * E.g. this is not safe:
         *    {
         *            int foo;
         *            pool.add_task_lambda([&foo](worker_thread*) { frob(foo); });
         *    }
         */
        template<typename T>
        void add_task_lambda(T f, const int id = 0, const size_t queue = 0, double cost = 0.0)
        {
            std::function<task_result*(worker_thread*)> task_fn;

            if constexpr (std::is_invocable_v<T, worker_thread*, int>) {
                task_fn = [f = std::move(f), id](worker_thread* w) -> task_result* {
                    f(w, id);
                    return new task_result;
                };
            } else {
                task_fn = [f = std::move(f)](worker_thread* w) -> task_result* {
                    f(w);
                    return new task_result;
                };
            }

            enqueue_task(std::move(task_fn), id, queue, cost);
        }
        /* }}} */

        /* {{{ add task_class.
         *
         * This creates a copy ff of the class object f of type T, and eventually
         * calls ff(worker, id), deleting ff afterwards. As f itself is copied,
         * only the caller has to care about its deletion.
         *
         * This interface is somewhat less useful than the next one, because
         * there is only limited potential for using the id argument.
         * Furthermore, it happily duplicates the argument descriptors.
         */
    private:
        template<typename T>
            static task_result * call_class_operator(worker_thread * worker, task_parameters * _param, int id) {
                auto clean_param = call_dtor([_param]() { delete _param; });
                (*static_cast<T*>(_param))(worker, id);
                return new task_result;
            }
    public:
        template<typename T>
            void add_task_class(T const & f, const int id, const size_t queue = 0, double cost = 0.0)
            {
                static_assert(std::is_base_of_v<task_parameters, T>, "type must inherit from task_parameters");
                add_task(thread_pool::call_class_operator<T>, new T(f), id, queue, cost);
            }
        /* }}} */

        /* {{{ add_shared_task - wraps std::shared_ptr inside a lambda */
        template<typename T>
            using shared_task = std::shared_ptr<T>;

        template<typename T, typename... Args>
            static shared_task<T> make_shared_task(Args&&... args) {
                return std::make_shared<T>(std::forward<Args>(args)...);
            }

        template<typename T>
            void add_shared_task(shared_task<T> const & task_ptr, const int id = 0, const size_t queue = 0, double cost = 0.0)
            {
                add_task_lambda([task_ptr](worker_thread* w, int task_id) {
                        (*task_ptr)(w, task_id);
                        }, id, queue, cost);
            }
        /* }}} */
};

#endif  /* CADO_THREADPOOL_HPP */
