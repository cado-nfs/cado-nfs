#ifndef LINGEN_MEMORY_POOL_HPP_
#define LINGEN_MEMORY_POOL_HPP_
// IWYU pragma: no_include <bits/exception.h>
#include <cstdint>          // for SIZE_MAX
#include <new>               // for bad_alloc
#include <mutex>
#include <cstdlib>
#include <cstdio>
#include <exception> // std::exception // IWYU pragma: keep
#include <string>

#define MEMORY_POOL_ALLOC_CHECK(X) memory_pool_details::alloc_check(#X, (X))

namespace memory_pool_details {
    void alloc_check(const char * text, bool condition);
    template<bool loose> struct inaccuracy_handler {};
    template<>
        struct inaccuracy_handler<false> {
            void handle_expand(size_t already_allocated, size_t asked, size_t & previously_allowed) {
                MEMORY_POOL_ALLOC_CHECK(already_allocated + asked <= previously_allowed);
            }
        };
    template<>
        struct inaccuracy_handler<true> {
            size_t cumulated_inaccuracy = 0;
            void handle_expand(size_t already_allocated, size_t asked, size_t & previously_allowed);
        };
}

class memory_pool_exception : public std::exception {
    std::string message;
public:
    memory_pool_exception(std::string const & s);
    virtual const char * what() const noexcept override { return message.c_str(); }
};

template<bool loose = false>
struct memory_pool : public memory_pool_details::inaccuracy_handler<loose> {
    std::mutex mm;
    size_t allowed=0;
    size_t allocated=0;
    size_t peak=0;
    friend class guard_base;
    public:
    void * alloc(size_t s)
    {
        std::lock_guard<std::mutex> dummy(mm);
        memory_pool_details::inaccuracy_handler<loose>::handle_expand(allocated, s, allowed);
        allocated += s;
        if (allocated > peak) peak = allocated;
        void * p = malloc(s);
        if (!p)
            throw std::bad_alloc();
        return p;
    }
    void free(void * p, size_t s)
    {
        std::lock_guard<std::mutex> dummy(mm);
        MEMORY_POOL_ALLOC_CHECK(allocated >= s);
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
        void * x = ::realloc(p, ns);
        if (!x) {
            fprintf(stderr, "Throwing std::bad_alloc after realloc %zu -> %zu\n", s, ns);
            // should we do that or not ? realloc itself won't, on
            // failure. And we're going to die anyway, so...
            // free(p);
            throw std::bad_alloc();
        }
        return x;
    }

    private:
    void report_inaccuracy(size_t diff);

    public:
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
                MEMORY_POOL_ALLOC_CHECK(memory.allocated == 0);
            memory.peak = 0;
        }
        void pre_dtor(memory_pool & memory) {
            if (oldsize == 0)
                MEMORY_POOL_ALLOC_CHECK(memory.allocated == 0);
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
            MEMORY_POOL_ALLOC_CHECK(memory.allocated <= memory.allowed);
        }
    };
};

typedef memory_pool<true> memory_pool_loose;
typedef memory_pool<false> memory_pool_strict;

template<typename ptr, bool is_loose>
class memory_pool_wrapper : public memory_pool<is_loose> {
    typedef memory_pool<is_loose> super;
public:
    inline ptr alloc(size_t s) { return (ptr) super::alloc(s); }
    void free(ptr p, size_t s) { super::free((void*) p, s); }
    ptr realloc(ptr p, size_t s, size_t ns) { return (ptr) super::realloc(p, s, ns); }
};

#endif	/* LINGEN_MEMORY_POOL_HPP_ */
