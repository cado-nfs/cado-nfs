#ifndef THREADPOOL_H_
#define THREADPOOL_H_

#include <pthread.h>
#include <queue>
#include <vector>
#include <exception>
#include <stdarg.h>
#include <errno.h>

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
class worker_thread;
class tasks_queue;
class results_queue;
class exceptions_queue;
class thread_pool;

typedef task_result *(*task_function_t)(const worker_thread * worker, const task_parameters *);

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
  mutable timetree_t timer;
  int rank() const;
  worker_thread(thread_pool &, size_t);
  ~worker_thread();
  void push_chart_item(time_bubble_chaser const & tt) const { timer.chart.push_back(tt); }
};


class thread_pool : private monitor, private NonCopyable {
  friend class worker_thread;

  std::vector<worker_thread> threads;
  std::vector<tasks_queue> tasks;
  std::vector<results_queue> results;
  std::vector<exceptions_queue> exceptions;

  bool kill_threads; /* If true, hands out kill tasks once work queues are empty */

  static void * thread_work_on_tasks(void *pool);
  thread_task *get_task(size_t queue);
  void add_result(size_t queue, task_result *result);
  void add_exception(size_t queue, clonable_exception * e);
  bool all_task_queues_empty() const;
public:
  void push_chart_item(int rank, time_bubble_chaser const & tt) { threads[rank].timer.chart.push_back(tt); }
  void clear_time_charts() { for(auto & T : threads) T.timer.chart.clear(); }
  void display_time_charts() const;
  thread_pool(size_t _nr_threads, size_t nr_queues = 1);
  ~thread_pool();
  void add_task(task_function_t func, const task_parameters * params, const int id, const size_t queue = 0, double cost = 0.0);
  task_result *get_result(size_t queue = 0, bool blocking = true);
  clonable_exception * get_exception(const size_t queue = 0);

  /* All threads in a thread pool have their respective timer active at
   * all times. We have two different ways to collect the timings they've
   * recorder over the course of their execution. Those depend on whether
   * we wish to count the times spent on _tasks_ (active time), or the
   * time spent _waiting_ for tasks (wait time). This is the purpose of
   * the two methods below.
   */

  /*
   * For accumulate_and_clear_active_time, we expect that all threads are
   * currently waiting, and all recorded child timings are transferred to
   * [rep]. The worker thread's timer is clear of all of its child
   * timings, but its wait time is unchanged.
   */
  void accumulate_and_clear_active_time(timetree_t::super & rep);

  /* For accumulate_and_reset_wait_time, we also expect that all threads
   * are currently waiting, but we also mandate that none has any
   * recorded child timings. That is, the function
   * accumulate_and_clear_active_time() must have been called previously
   * to stow the worker threads' children timings somewhere.
   *
   * accumulate_and_reset_wait_time counts the outstanding wait time for
   * the worker thread, and accumulates it to the target timer. After the
   * call, all threads are waiting again, but with their wait time reset
   * to zero (as if they had just started afresh).
   */
  void accumulate_and_reset_wait_time(timetree_t::super & rep);
};

#endif
