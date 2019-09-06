#ifndef LINGEN_MEMORY_POOL_HPP_
#define LINGEN_MEMORY_POOL_HPP_

#include <mutex>
#include <cstdlib>
#include "macros.h"

namespace memory_pool_details {
    template<bool loose> struct inaccuracy_handler {};
    template<>
        struct inaccuracy_handler<false> {
            void handle_expand(size_t already_allocated, size_t asked, size_t & previously_allowed) {
                ASSERT_ALWAYS(already_allocated + asked <= previously_allowed);
            }
        };
    template<>
        struct inaccuracy_handler<true> {
            size_t cumulated_inaccuracy = 0;
            void handle_expand(size_t already_allocated, size_t asked, size_t & previously_allowed);
        };
}


template<bool loose = false>
struct memory_pool : public memory_pool_details::inaccuracy_handler<loose> {
    std::mutex mm;
    public:
    size_t allowed=0;
    size_t allocated=0;
    size_t peak=0;
    void * alloc(size_t s)
    {
        std::lock_guard<std::mutex> dummy(mm);
        memory_pool_details::inaccuracy_handler<loose>::handle_expand(allocated, s, allowed);
        allocated += s;
        if (allocated > peak) peak = allocated;
        return malloc(s);
    }
    void free(void * p, size_t s)
    {
        std::lock_guard<std::mutex> dummy(mm);
        ASSERT_ALWAYS(allocated >= s);
        allocated -= s;
        ::free(p);
    }
    void * realloc(void * p, size_t s, size_t ns)
    {
        if (s == ns) return p;
        std::lock_guard<std::mutex> dummy(mm);
        /* We allow reallocating stuff that was not allocated at the present
         * recursive level */
        memory_pool_details::inaccuracy_handler<loose>::handle_expand(allocated, ns - s, allowed);
        allocated -= s;
        allocated += ns;
        return ::realloc(p, ns);
    }

    private:
    void report_inaccuracy(size_t diff);

    class guard_base {
        size_t oldsize;
        size_t mysize;
        public:
        guard_base(memory_pool & memory, size_t s) : mysize(s)
        {
            oldsize = memory.allowed;
            if (oldsize == SIZE_MAX || s == SIZE_MAX)
                memory.allowed = SIZE_MAX;
            else
                memory.allowed += s;
            if (oldsize == 0)
                ASSERT_ALWAYS(memory.allocated == 0);
            memory.peak = 0;
        }
        void pre_dtor(memory_pool & memory) {
            if (oldsize == 0)
                ASSERT_ALWAYS(memory.allocated == 0);
            if (memory.allowed != SIZE_MAX)
                memory.allowed -= mysize;
            else
                memory.allowed = oldsize;
            /* We don't do memory.allowed = oldsize unconditionally, because of
             * the inaccuracy tolerance:
             *  - guards that set the memory loose (and subsequent
             *  guards, if any) will never trigger a change of the
             *  .allowed field, of course, and will never report
             *  inaccuracy.  For these, setting to oldsize is naturally
             *  the only thing to do.
             *  - on the other hand, more precise guards may see
             *  memory.allowed grow above the initial value, which we
             *  still record as oldsize.  For these, we must not use
             *  oldsize. We simply subtract the amount of memory we had
             *  provisioned, and the growth that happened since the ctor
             *  will continue to accumulate, which is what we want.
             */
            ASSERT_ALWAYS(memory.allocated <= memory.allowed);
        }
    };

    public:
    template<typename T>
        struct guard : private guard_base {
            guard(size_t s) : guard_base(T::memory, s) {}
            ~guard() { guard_base::pre_dtor(T::memory); }
        };

};

typedef memory_pool<true> memory_pool_loose;
typedef memory_pool<false> memory_pool_strict;

#if 0
/* unused */
/* must be explicitly instantiated */
template<typename T> struct memory_guard;
#endif

#endif	/* LINGEN_MEMORY_POOL_HPP_ */
