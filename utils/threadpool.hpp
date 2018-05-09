#ifndef THREADPOOL_H_
#define THREADPOOL_H_

#include <pthread.h>
#include <queue>
#include <vector>
#include <exception>
#include <stdarg.h>
#include <errno.h>
#include <memory>

#include "macros.h"
#include "utils_cxx.hpp"
#include "tdict.hpp"

#if 0
/* Verbosely log all mutex and condition variable operations */
#include <typeinfo>
#include <verbose.h>
static inline void
thread_log(const char *c, const char *m, const void *p)
{
  unsigned long id;
  const pthread_t pid = pthread_self();
  memcpy(&id, &pid, MIN(sizeof(id), sizeof(pid)));
  verbose_output_print(0, 0, "Thread %lx: %s.%s(%p)\n", id, c, m, p);
  fflush(stdout);
}
#define LOG do{thread_log(typeid(this).name(), __func__, this);}while(0)
#else
#define LOG
#endif

/* C++11 already has classes for mutex and condition_variable */
/* All the synchronization stuff could be moved to the implementation if
   thread_pool used monitor as a dynamically allocated object. Tempting. */

class mutex {
  friend class condition_variable;
  pthread_mutex_t m;
  public:
  mutex(const pthread_mutexattr_t *mutexattr = NULL) {
    LOG;
    pthread_mutex_init(&m, mutexattr);
  }
  ~mutex() {LOG; ASSERT_ALWAYS(pthread_mutex_destroy(&m) == 0);}
  void lock(){LOG; ASSERT_ALWAYS(pthread_mutex_lock(&m) == 0);}
  bool try_lock() {
    LOG;
    int rc = pthread_mutex_trylock(&m);
    ASSERT_ALWAYS(rc == 0 || rc == EBUSY);
    return rc;
  }
  void unlock(){LOG; ASSERT_ALWAYS(pthread_mutex_unlock(&m) == 0);}
};

class condition_variable {
  pthread_cond_t cv;
  public:
  condition_variable(pthread_condattr_t *cond_attr = NULL) {
    LOG;
    ASSERT_ALWAYS(pthread_cond_init(&cv, cond_attr) == 0);
  }
  ~condition_variable() {LOG; pthread_cond_destroy(&cv);}
  void signal() {LOG; ASSERT_ALWAYS(pthread_cond_signal(&cv) == 0);}
  void broadcast(){LOG; ASSERT_ALWAYS(pthread_cond_broadcast(&cv) == 0);}
  void wait(mutex &m) {LOG; ASSERT_ALWAYS(pthread_cond_wait(&cv, &m.m) == 0);}
};

class monitor {
  mutex m;
public:
  void enter() {m.lock();}
  void leave() {m.unlock();}
  void signal(condition_variable &cond) {cond.signal();}
  void broadcast(condition_variable &cond){cond.broadcast();}
  void wait(condition_variable &cond) {cond.wait(m);}
};

/* Base for classes that hold parameters for worker functions */
class task_parameters {
  public:
  virtual ~task_parameters(){};
};

/* Base for classes that hold results produced by worker functions */
class task_result {
  public:
  virtual ~task_result(){};
};

struct clonable_exception : public std::exception {
    virtual clonable_exception * clone() const = 0;
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
  worker_thread(worker_thread const &) = delete;
  worker_thread& operator=(worker_thread const &) = delete;
public:
  worker_thread(worker_thread&&) = default;
  worker_thread& operator=(worker_thread&&) = default;
  int rank() const;
  int nthreads() const;
  /* It doesn't seem that unholy to me to have a thread access the pool
   * it originates from. It's possibly a good way to do continuations,
   * for example.
   */
  thread_pool & get_pool() { return pool; }
  worker_thread(thread_pool &, size_t);
  ~worker_thread();
};

typedef task_result *(*task_function_t)(worker_thread * worker, task_parameters *, int id);

class thread_pool : private monitor, private NonCopyable {
  friend class worker_thread;

  std::vector<worker_thread> threads;
  std::vector<tasks_queue> tasks;
  std::vector<results_queue> results;
  std::vector<exceptions_queue> exceptions;
  std::vector<size_t> created;
  std::vector<size_t> joined;

  bool kill_threads; /* If true, hands out kill tasks once work queues are empty */

  static void * thread_work_on_tasks(void *pool);
  thread_task get_task(size_t& queue);
  void add_result(size_t queue, task_result *result);
  void add_exception(size_t queue, clonable_exception * e);
  bool all_task_queues_empty() const;
public:
  void display_time_charts() const;
  static double cumulated_wait_time;
  thread_pool(size_t _nr_threads, size_t nr_queues = 1);
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
   *            pool.add_task_lambda([&foo](worker_thread*) { frob(foo); *            });
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
          task_parameters_lambda<T> *param = static_cast<task_parameters_lambda<T>*>(_param);
          param->f(worker, id);
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
          T *param = static_cast<T*>(_param);
          (*param)(worker, id);
          return new task_result;
      }
public:
  template<typename T>
      void add_task_class(T const & f, const int id, const size_t queue = 0, double cost = 0.0)
      {
          static_assert(std::is_base_of<task_parameters, T>::value, "type must inherit from task_parameters");
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
   * parameter object, there's not much we can do but to create an extra
   * level of indirection, which kinds of defeats the purpose...
   *
   * The next step toward making it more useful would be to convert
   * thread_task to also embed shared_ptr's.
   */

  template<typename T>
      struct shared_task : public std::shared_ptr<T>, public task_parameters {
          typedef std::shared_ptr<T> super;
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
          thread_pool::shared_task<T> *param = static_cast<thread_pool::shared_task<T>*>(_param);
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

#endif
