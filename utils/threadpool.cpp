#include "cado.h"
#include <string.h>
#include <stdio.h>
#include "threadpool.hpp"
#include "barrier.h"

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

worker_thread::worker_thread(thread_pool &_pool, const size_t _preferred_queue)
  : pool(_pool), preferred_queue(_preferred_queue)
{
  int rc = pthread_create(&thread, NULL, pool.thread_work_on_tasks, this);
  ASSERT_ALWAYS(rc == 0);
}

worker_thread::~worker_thread() {
  int rc = pthread_join(thread, NULL);
  ASSERT_ALWAYS(!timer.running());
  ASSERT_ALWAYS(timer.M.size() == 0);
  // timer.self will be essentially lost here. the best thing to do is to
  // create a phony task which collects the busy wait times for all
  // threads, at regular intervals, so that timer.self will be
  // insignificant.
  // pool.timer += timer;
  ASSERT_ALWAYS(rc == 0);
}

int worker_thread::rank() const { return this - pool.threads; }

class thread_task {
public:
  const task_function_t func;
  const int id;
  const task_parameters * parameters;
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

class exceptions_queue : public std::queue<clonable_exception *>, private ThreadNonCopyable {
  public:
  condition_variable not_empty;
};


thread_pool::thread_pool(const size_t _nr_threads, const size_t _nr_queues)
  : nr_threads(_nr_threads), nr_queues(_nr_queues), kill_threads(false)
{
  /* Threads start accessing the queues as soon as they run */
  tasks = new tasks_queue[nr_queues];
  results = new results_queue[nr_queues];
  exceptions = new exceptions_queue[nr_queues];

  threads = reinterpret_cast<worker_thread *>(malloc(nr_threads * sizeof(worker_thread)));
  for (size_t i = 0; i < nr_threads; i++)
    new (threads + i) worker_thread(*this, 0);
};

thread_pool::~thread_pool() {
  enter();
  kill_threads = true;
  /* Wakey wakey, time to die */
  for (size_t i = 0; i < nr_queues; i++)
    broadcast(tasks[i].not_empty);
  leave();
  for (size_t i = 0; i < nr_threads; i++) {
      threads[i].~worker_thread();
  }
  free(reinterpret_cast<void*>(threads));
  for (size_t i = 0; i < nr_queues; i++)
    ASSERT_ALWAYS(tasks[i].empty());
  delete[] tasks;
  for (size_t i = 0; i < nr_queues; i++)
    ASSERT_ALWAYS(results[i].empty());
  delete[] results;
  for (size_t i = 0; i < nr_queues; i++)
    ASSERT_ALWAYS(exceptions[i].empty());
  delete[] exceptions;
}

void *
thread_pool::thread_work_on_tasks(void *arg)
{
  worker_thread *I = (worker_thread *) arg;
  ACTIVATE_TIMER(I->timer);
  while (1) {
    thread_task *task = I->pool.get_task(I->preferred_queue);
    if (task->please_die) {
      delete task;
      break;
    }
    task_function_t func = task->func;
    const task_parameters *params = task->parameters;
    try {
        task_result *result = func(I, params);
        if (result != NULL)
            I->pool.add_result(task->queue, result);
    } catch (clonable_exception const& e) {
        I->pool.add_exception(task->queue, e.clone());
        /* We need to wake the listener... */
        I->pool.add_result(task->queue, NULL);
    }
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
thread_pool::add_task(task_function_t func, const task_parameters * params,
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

void
thread_pool::add_exception(const size_t queue, clonable_exception * e) {
  ASSERT_ALWAYS(queue < nr_queues);
  enter();
  exceptions[queue].push(e);
  // do we use it ?
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

/* get an exception from the specified exceptions queue. This is
 * obviously non-blocking, because exceptions are exceptional. So when no
 * exception is there, we return NULL. When there is one, we return a
 * pointer to a newly allocated copy of it.
 */
clonable_exception * thread_pool::get_exception(const size_t queue) {
    clonable_exception * e;
    ASSERT_ALWAYS(queue < nr_queues);
    enter();
    if (exceptions[queue].empty()) {
        e = NULL;
    } else {
        e = exceptions[queue].front();
        exceptions[queue].pop();
    }
    leave();
    return e;
}

void thread_pool::accumulate_and_clear_active_time(timetree_t & rep) {
    for (size_t i = 0; i < nr_threads; ++i) {
        /* timers may be running when they're tied to a subthread which
         * is currently doing a pthread_cond_wait or such. So in that
         * case, we want it to continue keeping track on its cpu busy
         * time, while we are interested in the time the thread has spent
         * on the tasks it was assigned.
         *
         * On the other hand, there's a fairly simple use case of
         * merging with a stopped timer. In that case, we'll simply use
         * the += operator.
         *
         * It's quite possible that the tdict::operator+= should have some
         * distinguishing capacity based along these lines, so that the
         * logic here could me moved there.
         */
        if (threads[i].timer.running()) {
            ASSERT_ALWAYS(&threads[i].timer == threads[i].timer.current);
            rep.steal_children_timings(threads[i].timer);
        } else {
            // ASSERT_ALWAYS(0);
            rep += threads[i].timer;
            threads[i].timer = timetree_t ();
        }
    }
}

struct everybody_must_do_that : public task_parameters {
    mutable barrier_t barrier;
    everybody_must_do_that(int nthreads) { barrier_init(&barrier, nthreads); }
};
struct everybody_must_do_that_result : public task_result {
    double v;
    everybody_must_do_that_result(timetree_t & me) : v(me.stop_and_start()){}
};

task_result * everybody_must_do_that_task(const worker_thread * worker, const task_parameters * _param)
{
    const everybody_must_do_that * param = dynamic_cast<const everybody_must_do_that*>(_param);
    barrier_wait(&param->barrier, NULL, NULL, NULL);
    return new everybody_must_do_that_result(worker->timer);
}

void thread_pool::accumulate_and_reset_wait_time(timetree_t & rep) {
    /* we need to create a task so that each thread does what we want it
     * to do. The tricky part is the we really want all thread to block
     * and reach this code. In effect, we want the callee function to
     * embody a barrier_wait
     */
    everybody_must_do_that E(nr_threads);
    for(unsigned int i = 0 ; i < nr_threads ; i++) {
        add_task(everybody_must_do_that_task, &E, 0);
    }
    for(unsigned int i = 0 ; i < nr_threads ; i++) {
        everybody_must_do_that_result * res = dynamic_cast<everybody_must_do_that_result*>(get_result());
        rep.self += res->v;
        delete res;
    }
}
