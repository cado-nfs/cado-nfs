#include "cado.h" // IWYU pragma: keep

#include <cstring>
#include <cstdlib>
#include <cstdio>

#include <thread>
#include <deque>
#include <string>
#include <utility>

#include "fmt/base.h"
#include "fmt/format.h"
#include "fmt/std.h"

#include "tests_common.h"
#include "threadpool.hpp"

// coverity[root_function]
int main(int argc, char const * argv[])
{
  tests_common_cmdline(&argc, &argv, PARSE_ITER);
  unsigned long iter = 10;
  tests_common_get_iter(&iter);
  double wait_time = 0;

  thread_pool pool(5, wait_time, 2);

  std::deque<std::future<size_t>> results;

  const char * message = "Hello world";

  for (unsigned long i = 0; i < iter; i++) {
    const size_t queue = i % 2;
    auto fut = pool.add_future_task(queue, 0.0, [message](worker_thread*) {
                auto s = fmt::format("This is thread {}: {}",
                        std::this_thread::get_id(), message);
                fmt::print("{}\n", s);
                return s.size();
                });

    results.push_back(std::move(fut));
  }

  while(!results.empty()) {
      auto rc = results.front().get();
      results.pop_front();
      fmt::print("I've printed {} characters\n", rc);
  }

  return EXIT_SUCCESS;
}
