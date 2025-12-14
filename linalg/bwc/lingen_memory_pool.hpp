#ifndef CADO_LINGEN_MEMORY_POOL_HPP
#define CADO_LINGEN_MEMORY_POOL_HPP

#include <cstdint>
#include <cstdlib>
#include <cstdio>

#include <new>
#include <mutex>
#include <exception>
#include <string>

/* A memory pool is an object from which heap storage may be allocated.
 * No surprises here. The memory pool keeps track of the memory that is
 * allocated and freed.
 * The memory pool only tolerates allocation up to a certain limit. While
 * the memory pool is generally a singleton object, it is possible to
 * scope one or several "memory guards" that refer to this memory pool.
 * Each memory guard adds up to the total accepted allocation size of the
 * memory pool, and this memory is expected to be freed when the guard
 * goes out of scope.
 */
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
    memory_pool_exception(std::string);
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
        guard_base(guard_base const &) = delete;
        guard_base& operator=(guard_base const &) = delete;

#if 0
        guard_base(guard_base && g)
        {
            old_size = g.old_size;
            my_size = g.my_size;
            /* then we should modify the g fields so that its pre-dtor
             * call becomes a no-op. Unfortunately, the interface doesn't
             * let us do that, which is a very clear misfeature!
             */
        }
        guard_base& operator=(guard_base &&) = delete;
#else
        guard_base(guard_base && g) = delete;
        guard_base& operator=(guard_base &&) = delete;
#endif
    };
};

using memory_pool_loose = memory_pool<true>;
using memory_pool_strict = memory_pool<false>;

template<typename ptr, bool is_loose>
class memory_pool_wrapper : public memory_pool<is_loose> {
    using super = memory_pool<is_loose>;
public:
    ptr alloc(size_t s) { return (ptr) super::alloc(s); }
    void free(ptr p, size_t s) { super::free((void*) p, s); }
    ptr realloc(ptr p, size_t s, size_t ns) { return (ptr) super::realloc(p, s, ns); }
};

#endif	/* LINGEN_MEMORY_POOL_HPP_ */
