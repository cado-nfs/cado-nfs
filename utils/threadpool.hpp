#ifndef CADO_THREADPOOL_HPP
#define CADO_THREADPOOL_HPP

#include <cerrno>        // for EBUSY
#include <cstddef>       // for size_t, NULL

#include <memory>         // for shared_ptr, make_shared
#include <mutex>          // for mutex
#include <type_traits>    // for is_base_of
#include <vector>         // for vector
#include <condition_variable>

#include <pthread.h>      // for pthread_cond_broadcast, pthread_cond_destroy

#include "macros.h"       // for ASSERT_ALWAYS
#include "utils_cxx.hpp"  // for call_dtor, NonCopyable

struct clonable_exception;

/* C++11 already has classes for mutex and condition_variable */
/* All the synchronization stuff could be moved to the implementation if
   thread_pool used monitor as a dynamically allocated object. Tempting. */

class monitor {
protected:
  std::mutex m;
public:
  struct my_unique_lock : public std::unique_lock<std::mutex> {
      explicit my_unique_lock(monitor & m) : std::unique_lock<std::mutex>(m.m) {}
  };
  // void enter() {m.lock();}
  // void leave() {m.unlock();}
  static void signal(std::condition_variable &cond) {cond.notify_one();}
  static void broadcast(std::condition_variable &cond){cond.notify_all();}
  static void wait(std::condition_variable &cond, my_unique_lock & u) {cond.wait(u);}
};

class monitor_or_synchronous : public monitor {
    bool sync = false;
    public:
    struct my_unique_lock : public std::unique_lock<std::mutex> {
        explicit my_unique_lock(monitor_or_synchronous & m)
            : std::unique_lock<std::mutex>(m.m, std::defer_lock_t())
        {
            if (!m.sync) lock();
        }
    };
    explicit monitor_or_synchronous(bool sync = false)
        : sync(sync)
    {}
    ATTRIBUTE_NODISCARD bool is_synchronous() const { return sync; }
    // void enter() { if (!sync) monitor::enter();}
    // void leave() { if (!sync) monitor::leave();}
    void signal(std::condition_variable &cond) const { if (!sync) monitor::signal(cond); }
    void broadcast(std::condition_variable &cond) const { if (!sync) monitor::broadcast(cond);}
    void wait(std::condition_variable &cond, my_unique_lock & u) const {
        if (sync)
            ASSERT_ALWAYS(0);
        else
            cond.wait(u);
    }
};

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
  size_t rank() const;
  size_t nthreads() const;
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

class thread_pool : private monitor_or_synchronous, private NonCopyable {
  friend class worker_thread;

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
public:
  // bool is_synchronous() const { return monitor_or_synchronous::is_synchronous(); }
  double cumulated_wait_time = 0;
  std::mutex mm_cumulated_wait_time;

  thread_pool(size_t _nr_threads, double & store_wait_time,
          size_t nr_queues = 1, bool sync_thread_pool = false);
  ~thread_pool();
  task_result *get_result(size_t queue = 0, bool blocking = true);
  void drain_queue(const size_t queue, void (*f)(task_result*) = NULL);
  void drain_all_queues();
  clonable_exception * get_exception(const size_t queue = 0);
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
  void add_task(task_function_t func, task_parameters * params, const int id, const size_t queue = 0, double cost = 0.0);
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
private:
  template<typename T>
      struct task_parameters_lambda : public task_parameters {
          T f;
          task_parameters_lambda(T const& f) : f(f) {}
      };
  template<typename T>
      static
      task_result * do_task_parameters_lambda(worker_thread * worker, task_parameters * _param, int id) {
          auto clean_param = call_dtor([_param]() { delete _param; });
          static_cast<task_parameters_lambda<T>*>(_param)->f(worker, id);
          return new task_result;
      }
public:
  template<typename T>
      void add_task_lambda(T f, const int id, const size_t queue = 0, double cost = 0.0)
      {
          add_task(thread_pool::do_task_parameters_lambda<T>, new task_parameters_lambda<T>(f), id, queue, cost);
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

#if 1
  /* {{{ add_shared_task -- NOT SATISFACTORY. Do not use.
   *
   * This last interface is an attempt at being more useful. We would
   * like to solve the ownership conflict that lurks in the design of
   * add_task. Here we explicitly say that ownership of the T object is
   * shared between the caller which creates it, and the (one or several)
   * task(s) that are to use it. The caller does not have to join all
   * threads before the object of type shared_ptr<T> goes out of scope,
   * since the pool queue will still have enough referenced items to
   * guarantee that the object stays alive.
   *
   * Here, proper use of the id field can lead to efficient data sharing,
   * albeit at the expense of the shared_ptr management.
   *
   * Alas, since the thread_task objects only have raw pointers to the
   * parameter object, there's not much we can do but create an extra
   * level of indirection, which kinds of defeats the purpose...
   *
   * The next step toward making it more useful would be to convert
   * thread_task to also embed shared_ptr's.
   */

  template<typename T>
      struct shared_task : public std::shared_ptr<T>, public task_parameters {
          using super = std::shared_ptr<T>;
          shared_task(super c) : super(c) {}
          T& operator*() { return *(super&)(*this); }
          T const & operator*() const { return *(super const&)(*this); }
      };
  template<typename T, typename... Args>
  static shared_task<T> make_shared_task(Args&&... args) { return shared_task<T>(std::make_shared<T>(args...)); }

private:
  template<typename T>
      static task_result * call_shared_task(worker_thread * worker, task_parameters * _param, int id) {
          auto clean_param = call_dtor([_param]() { delete _param; });
          auto * param = static_cast<thread_pool::shared_task<T>*>(_param);
          (**param)(worker, id);
          return new task_result;
      }
public:
  template<typename T>
      void add_shared_task(shared_task<T> const & f, const int id, const size_t queue = 0, double cost = 0.0)
      {
          add_task(thread_pool::call_shared_task<T>, new shared_task<T>(f), id, queue, cost);
      }
  /* }}} */
#endif
};

#endif  /* CADO_THREADPOOL_HPP */
