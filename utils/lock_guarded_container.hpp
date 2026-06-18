#ifndef CADO_LOCK_GUARDED_CONTAINER_HPP
#define CADO_LOCK_GUARDED_CONTAINER_HPP

#include <mutex>

/* This is a bundle of a container (well, in fact, of anything) plus a
 * lock. The copy/move work just the same as they do (or don't) in the
 * original container, while a lock is created afresh for each new
 * object. Moves don't move the mutexes.
 */

template<typename T> struct lock_guarded_container : public T {
    private:
    mutable std::mutex mm;
    T& ctr() { return *this; }
    T const & ctr() const { return *this; }
    public:
    std::mutex & mutex() const { return mm; }
    /* forward the constructors of the embedded container. No lock needed
     * here. */
    template<typename... Args> explicit lock_guarded_container(Args&& ...args)
        : T { std::forward<Args>(args)... }
    {}
    /* add a few of our own to circumvent the lack of copy/move for the
     * mutex.
     * We considered having a dedicated (private) ctor taking a transient
     * lock_guard object, so that we can avoid the cost of
     * constructing+assigning, and construct the T base simply in the
     * initializer list by delegating to this private ctor. Alas, this is
     * incompatible with our perfect forwarding of the constructors of T.
     * I don't know if a workaround is possible.
     */
    lock_guarded_container(lock_guarded_container<T> const & o)
    {
        std::scoped_lock const foo(o.mutex());
        ctr() = o.ctr();
    }
    lock_guarded_container(lock_guarded_container<T> && o) noexcept
    {
        std::scoped_lock const foo(o.mutex());
        std::scoped_lock const bar(mutex());
        std::swap(ctr(), o.ctr());
    }
    lock_guarded_container& operator=(lock_guarded_container<T> const & o) {
        if (this == &o) return *this;
        std::scoped_lock const foo(o.mutex());
        std::scoped_lock const bar(mutex());
        ctr() = o.ctr();
        return *this;
    }
    lock_guarded_container& operator=(lock_guarded_container<T> && o) noexcept {
        std::scoped_lock const foo(o.mutex());
        std::scoped_lock const bar(mutex());
        std::swap(ctr(), o.ctr());
        return *this;
    }
    ~lock_guarded_container() = default;
    /* The .locked() proxy object provides locked access up to the next
     * sequence point where the proxy object is destroyed.
     */
    struct proxy {
        std::scoped_lock<std::mutex> lk;
        T * ctr;
        T * operator->() { return ctr; };
        T & operator*() { return *ctr; };
        proxy(std::mutex & mm, T * ctr) : lk(mm), ctr(ctr) {}
    };
    proxy locked() { return { mm, this }; }
    struct const_proxy {
        std::scoped_lock<std::mutex> lk;
        T const * ctr;
        T const * operator->() { return ctr; };
        T const & operator*() { return *ctr; };
        const_proxy(std::mutex & mm, T const * ctr) : lk(mm), ctr(ctr) {}
    };
    const_proxy locked() const { return { mm, this }; }
};

#endif	/* CADO_LOCK_GUARDED_CONTAINER_HPP */
