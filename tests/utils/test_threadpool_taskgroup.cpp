#include "cado.h" // IWYU pragma: keep

#include <atomic>
#include <chrono>
#include <cstdlib>
#include <thread>

#include "fmt/base.h"

#include "tests_common.h"
#include "threadpool.hpp"
#include "macros.h"

static void test_wait_ordering(thread_pool & pool)
{
    fmt::print("Testing task_group::wait() ordering relation...\n");
    task_group tg;
    constexpr int ntasks = 20;
    std::atomic<int> completed_tasks{0};

    for (int i = 0; i < ntasks; ++i) {
        pool.add_task(tg, [&completed_tasks](worker_thread *) {
            std::this_thread::sleep_for(std::chrono::milliseconds(2));
            completed_tasks.fetch_add(1, std::memory_order_relaxed);
        });
    }

    ASSERT_ALWAYS(tg.created() == static_cast<size_t>(ntasks));

    // Calling wait() MUST block until all ntasks have finished
    tg.wait();

    ASSERT_ALWAYS(completed_tasks.load() == ntasks);
    ASSERT_ALWAYS(tg.joined() == static_cast<size_t>(ntasks));
    fmt::print("  -> Passed task_group::wait() ordering test.\n");
}

static void test_on_complete_callback(thread_pool & pool)
{
    fmt::print("Testing task_group::on_complete() ordering relation...\n");
    task_group tg1;
    task_group tg2;
    constexpr int ntasks = 16;
    std::atomic<int> phase1_count{0};
    std::atomic<int> phase2_count{0};
    std::atomic<bool> callback_ran{false};
    std::atomic<bool> phase1_done_before_callback{false};

    for (int i = 0; i < ntasks; ++i) {
        pool.add_task(tg1, [&phase1_count](worker_thread *) {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
            phase1_count.fetch_add(1, std::memory_order_relaxed);
        });
    }

    // Register on_complete callback on tg1.
    // When tg1 finishes, callback should run and spawn tg2 tasks.
    tg1.on_complete([&pool, &tg2, &phase1_count, &phase2_count, &callback_ran, &phase1_done_before_callback]() {
        phase1_done_before_callback = (phase1_count.load() == ntasks);
        callback_ran = true;

        for (int i = 0; i < ntasks; ++i) {
            pool.add_task(tg2, [&phase2_count](worker_thread *) {
                phase2_count.fetch_add(1, std::memory_order_relaxed);
            });
        }

        /* Test adding to the same task group in the completion callback
         */
        tg2.on_complete([&pool, &tg2, &phase2_count]() {
            for (int i = 0; i < ntasks; ++i) {
                pool.add_task(tg2, [&phase2_count](worker_thread *) {
                    phase2_count.fetch_add(1, std::memory_order_acq_rel);
                });
            }
        });
    });

    // Wait on tg1 (which executes callback and enqueues tg2 tasks), then
    // wait on tg2
    tg1.wait();
    tg2.wait();

    ASSERT_ALWAYS(callback_ran.load());
    ASSERT_ALWAYS(phase1_done_before_callback.load());
    auto n1 = phase1_count.load();
    auto n2 = phase2_count.load();
    fmt::print("{} {}\n", n1, n2);
    ASSERT_ALWAYS(n1 == ntasks);
    ASSERT_ALWAYS(n2 == 2 * ntasks);

    // Now test registering on_complete AFTER tasks have already finished
    std::atomic<bool> late_callback_ran{false};
    tg2.on_complete([&late_callback_ran]() {
        late_callback_ran = true;
    });

    ASSERT_ALWAYS(late_callback_ran.load());
    fmt::print("  -> Passed task_group::on_complete() ordering test.\n");
}

static void test_parallel_chains_ordering(thread_pool & pool)
{
    fmt::print("Testing two parallel chains with strict pre/post completion ordering...\n");
    constexpr int nsides = 2;
    constexpr int ntasks_per_side = 12;

    std::atomic<int> ds_active[nsides] = {0, 0};
    std::atomic<int> ds_completed[nsides] = {0, 0};
    std::atomic<int> fib_completed[nsides] = {0, 0};
    std::atomic<bool> callback_executed[nsides] = {false, false};

    task_group ds_tg[nsides];
    task_group fib_tg[nsides];

    for (int side = 0; side < nsides; ++side) {
        // Enqueue DS tasks for side
        for (int i = 0; i < ntasks_per_side; ++i) {
            pool.add_task(ds_tg[side], [side, &ds_active, &ds_completed](worker_thread *) {
                ds_active[side].fetch_add(1, std::memory_order_relaxed);
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
                ds_active[side].fetch_sub(1, std::memory_order_relaxed);
                ds_completed[side].fetch_add(1, std::memory_order_relaxed);
            });
        }

        // When DS for side completes, trigger FIB for side
        ds_tg[side].on_complete([&pool, &fib_tg, side, &ds_active, &ds_completed, &fib_completed, &callback_executed]() {
            // Check strict ordering: ALL DS tasks for this side MUST be
            // completed before callback
            ASSERT_ALWAYS(ds_completed[side].load() == ntasks_per_side);
            ASSERT_ALWAYS(ds_active[side].load() == 0);
            callback_executed[side].store(true);

            for (int i = 0; i < ntasks_per_side; ++i) {
                pool.add_task(fib_tg[side], [side, &ds_completed, &fib_completed](worker_thread *) {
                    // Check strict ordering: DS MUST be 100% finished
                    // when any FIB task runs
                    ASSERT_ALWAYS(ds_completed[side].load() ==
                            ntasks_per_side);
                    std::this_thread::sleep_for(std::chrono::milliseconds(1));
                    fib_completed[side].fetch_add(1, std::memory_order_relaxed);
                });
            }
        });
    }

    // Wait on DS then FIB task groups for both sides
    for (int side = 0; side < nsides; ++side) {
        ds_tg[side].wait();
        fib_tg[side].wait();
        ASSERT_ALWAYS(callback_executed[side].load());
        ASSERT_ALWAYS(ds_completed[side].load() == ntasks_per_side);
        ASSERT_ALWAYS(fib_completed[side].load() == ntasks_per_side);
    }

    fmt::print("  -> Passed two parallel chains strict ordering test.\n");
}

// coverity[root_function]
int main(int argc, char const * argv[])
{
    tests_common_cmdline(&argc, &argv, PARSE_ITER);
    unsigned long iter = 5;
    tests_common_get_iter(&iter);

    double wait_time = 0;
    thread_pool pool(8, wait_time, 2);

    for (unsigned long i = 0; i < iter; ++i) {
        fmt::print("Iteration {}/{}\n", i + 1, iter);
        test_wait_ordering(pool);
        test_on_complete_callback(pool);
        test_parallel_chains_ordering(pool);
    }

    return EXIT_SUCCESS;
}
