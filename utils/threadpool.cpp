#include "cado.h"
#include <string.h>
#include <stdio.h>
#include "threadpool.h"

/*
  With multiple queues, when new work is added to a queue, we need to be able
  to wake up one of the threads that prefer work from that queue. Thus we need
  multiple condition variables. If no threads that prefer work from that queue
  are currently waiting, we need to wake up some other thread.

  With k queues, we need k condition variables c[] and k semaphores s[].
  When a thread that prefers queue i waits for work, in increases s[i] and starts waiting on c[i].
  When a thread that was waiting is woken up, it decreases s[i].
  When work is added to queue j, it checks whether s[j] is non-zero:
    - if so, it signals c[j]
    - if not, it tests whether any other c[l] is non-zero
      - if so, it signals c[l]
      - if not, then no threads are currently sleeping, and nothing needs to be done

  We use a simple size_t variable as the semaphore; accesses are mutex-protected.
*/

class worker_thread : private ThreadNonCopyable {
  friend class thread_pool;
  thread_pool &pool;
  pthread_t thread;
  const size_t preferred_queue;
public:
  worker_thread(thread_pool &, size_t);
  ~worker_thread();
};

worker_thread::worker_thread(thread_pool &_pool, const size_t _preferred_queue)
  : pool(_pool), preferred_queue(_preferred_queue)
{
  int rc = pthread_create(&thread, NULL, pool.thread_work_on_tasks, this);
  ASSERT_ALWAYS(rc == 0);
}

worker_thread::~worker_thread() {
  int rc = pthread_join(thread, NULL);
  ASSERT_ALWAYS(rc == 0);
}

class thread_task {
public:
  const task_function_t func;
  const int id;
  const task_parameters * const parameters;
  const bool please_die;
  const size_t queue;
  const double cost; // costly tasks are scheduled first.
  task_result *result;

  thread_task(task_function_t _func, int _id, const task_parameters *_parameters, size_t _queue, double _cost) :
    func(_func), id(_id), parameters(_parameters), please_die(false), queue(_queue), cost(_cost), result(NULL) {};
  thread_task(bool _kill)
    : func(NULL), id(0), parameters(NULL), please_die(true), queue(0), cost(0.0), result(NULL) {
    ASSERT_ALWAYS(_kill);
  }
};

class thread_task_cmp
{
public:
  thread_task_cmp() {}
  bool operator() (const thread_task *x, const thread_task *y) const {
    if (x->cost < y->cost)
      return true;
    if (x->cost > y->cost)
      return false;
    // if costs are equal, compare ids (they should be distinct)
    return x->id < y->id;
  }
};

class tasks_queue : public std::priority_queue<thread_task *, std::vector<thread_task *>, thread_task_cmp>, private ThreadNonCopyable {
  public:
  condition_variable not_empty;
  size_t nr_threads_waiting;
  tasks_queue() : nr_threads_waiting(0){};
};

class results_queue : public std::queue<task_result *>, private ThreadNonCopyable {
  public:
  condition_variable not_empty;
};


thread_pool::thread_pool(const size_t _nr_threads, const size_t _nr_queues)
  : nr_threads(_nr_threads), nr_queues(_nr_queues), kill_threads(false)
{
  /* Threads start accessing the queues as soon as they run */
  tasks = new tasks_queue[nr_queues];
  results = new results_queue[nr_queues];

  threads = new worker_thread_ptr[nr_threads];
  for (size_t i = 0; i < nr_threads; i++)
    threads[i] = new worker_thread(*this, 0);
};

thread_pool::~thread_pool() {
  enter();
  kill_threads = true;
  /* Wakey wakey, time to die */
  for (size_t i = 0; i < nr_queues; i++)
    broadcast(tasks[i].not_empty);
  leave();
  for (size_t i = 0; i < nr_threads; i++)
    delete threads[i];
  delete[] threads;
  for (size_t i = 0; i < nr_queues; i++)
    ASSERT_ALWAYS(tasks[i].empty());
  delete[] tasks;
  for (size_t i = 0; i < nr_queues; i++)
    ASSERT_ALWAYS(results[i].empty());
  delete[] results;
}

void *
thread_pool::thread_work_on_tasks(void *arg)
{
  worker_thread *I = (worker_thread *) arg;
  while (1) {
    thread_task *task = I->pool.get_task(I->preferred_queue);
    if (task->please_die) {
      delete task;
      break;
    }
    task_function_t func = task->func;
    const task_parameters *params = task->parameters;
    task_result *result = func(params);
    if (result != NULL)
      I->pool.add_result(task->queue, result);
    delete task;
  }
  return NULL;
}

bool
thread_pool::all_task_queues_empty() const
{
  bool empty = true;
  for (size_t i = 0; i < nr_queues; i++)
    empty = empty && tasks[i].empty();
  return empty;
}


void
thread_pool::add_task(task_function_t func, const task_parameters *const params,
                      const int id, const size_t queue, double cost)
{
    ASSERT_ALWAYS(queue < nr_queues);
    enter();
    ASSERT_ALWAYS(!kill_threads);
    tasks[queue].push(new thread_task(func, id, params, queue, cost));

    /* Find a queue with waiting threads, starting with "queue" */
    size_t i = queue;
    if (tasks[i].nr_threads_waiting == 0) {
      for (i = 0; i < nr_queues && tasks[i].nr_threads_waiting == 0; i++) {}
    }
    /* If any queue with waiting threads was found, wake up one of them */
    if (i < nr_queues)
      signal(tasks[i].not_empty);
    leave();
}
  
thread_task *
thread_pool::get_task(const size_t preferred_queue)
{
  enter();
  while (!kill_threads && all_task_queues_empty()) {
    /* No work -> tell this thread to wait until work becomes available.
       We also leave the loop when the thread needs to die.
       The while() protects against spurious wake-ups that can fire even if
       the queue is still empty. */
    tasks[preferred_queue].nr_threads_waiting++;
    wait(tasks[preferred_queue].not_empty);
    tasks[preferred_queue].nr_threads_waiting--;
  }
  thread_task *task;
  if (kill_threads && all_task_queues_empty()) {
    task = new thread_task(true);
  } else {
    /* Find a non-empty task queue, starting with the preferred one */
    size_t i = preferred_queue;
    if (tasks[i].empty()) {
      for (i = 0; i < nr_queues && tasks[i].empty(); i++) {}
    }
    /* There must have been a non-empty queue or we'd still be in the while()
       loop above */
    ASSERT_ALWAYS(i < nr_queues);
    task = tasks[i].top();
    tasks[i].pop();
  }
  leave();
  return task;
}

void
thread_pool::add_result(const size_t queue, task_result *const result) {
  ASSERT_ALWAYS(queue < nr_queues);
  enter();
  results[queue].push(result);
  signal(results[queue].not_empty);
  leave();
}

/* Get a result from the specified results queue. If no result is available,
   waits with blocking=true, and returns NULL with blocking=false. */
task_result *
thread_pool::get_result(const size_t queue, const bool blocking) {
  task_result *result;
  ASSERT_ALWAYS(queue < nr_queues);
  enter();
  if (!blocking and results[queue].empty()) {
    result = NULL;
  } else {
    while (results[queue].empty())
      wait(results[queue].not_empty);
    result = results[queue].front();
    results[queue].pop();
  }
  leave();
  return result;
}
