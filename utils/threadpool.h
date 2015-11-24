#ifndef THREADPOOL_H_
#define THREADPOOL_H_

#include <pthread.h>
#include <queue>
#include <vector>
#include <stdarg.h>
#include <errno.h>

#include "macros.h"

class ThreadNonCopyable {
 protected:
   ThreadNonCopyable() {}
   ~ThreadNonCopyable() {}
 private:
   ThreadNonCopyable(const ThreadNonCopyable&);
   ThreadNonCopyable& operator=(const ThreadNonCopyable&);
};

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

typedef task_result *(*task_function_t)(const task_parameters *);

class thread_task;
class worker_thread;
class tasks_queue;
class results_queue;


class thread_pool : private monitor, private ThreadNonCopyable {
  friend class worker_thread;
  typedef worker_thread *worker_thread_ptr;

  worker_thread_ptr *threads;
  tasks_queue *tasks;
  results_queue *results;

  const size_t nr_threads, nr_queues;

  bool kill_threads; /* If true, hands out kill tasks once work queues are empty */

  static void * thread_work_on_tasks(void *pool);
  thread_task *get_task(size_t queue);
  void add_result(size_t queue, task_result *result);
  bool all_task_queues_empty() const;
public:
  thread_pool(size_t _nr_threads, size_t nr_queues = 1);
  ~thread_pool();
  void add_task(task_function_t func, const task_parameters *const params, const int id, const size_t queue = 0, double cost = 0.0);
  task_result *get_result(size_t queue = 0, bool blocking = true);
};

#endif
