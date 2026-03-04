#ifndef CADO_REPLACEMENT_BARRIER_HPP_
#define CADO_REPLACEMENT_BARRIER_HPP_

/* This is meant to provide a very basic equivalent for std::barrier for
 * libstdc++10, which doesn't have it (only version 11.1 has
 * std::barrier)
 *
 * Our old C version of a barrier did have in and out functions, so
 * presumably it's possible to work our way to use them and recover the
 * completion function mechanism.
 *
 * Likewise, _our_ main usage is arrive_and_wait and it's the only thing
 * we implement. arrive() and arrive_and_drop() are left as exercises.
 */

#include <cstddef>
#include <type_traits>

#include "barrier.h"

namespace cado_std_replacement
{
    struct foo {
        void operator()() const { }
    };
    static_assert(std::is_empty_v<foo>);

struct barrier
{
    barrier_t B {};

    explicit barrier(ptrdiff_t count)
    {
        barrier_init(&B, nullptr, static_cast<int>(count));
    }

    ~barrier()
    {
        barrier_destroy(&B, nullptr);
    }

    barrier(barrier const &) = delete;
    barrier & operator=(barrier const &) = delete;
    barrier(barrier const &&) = delete;
    barrier & operator=(barrier const &&) = delete;

    void arrive_and_wait() {
        barrier_wait(&B, nullptr, nullptr, nullptr);
    }
};
} // namespace cado_std_replacement

#ifndef __cpp_lib_barrier
namespace std {
    /* this will error out, by design, if a non-trivial completion
     * function is passed.
     */
    template<int=0>
    using barrier = cado_std_replacement::barrier;
}
#endif

#endif /* CADO_REPLACEMENT_BARRIER_HPP_ */
