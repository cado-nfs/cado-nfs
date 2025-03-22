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
        std::lock_guard<std::mutex> const foo(o.mutex());
        (T&)*this = (T const&) o;
    }
    lock_guarded_container(lock_guarded_container<T> && o) noexcept
    {
        std::lock_guard<std::mutex> const foo(o.mutex());
        std::lock_guard<std::mutex> const bar(mutex());
        std::swap((T&)*this, (T&) o);
    }
    lock_guarded_container& operator=(lock_guarded_container<T> const & o) {
        if (this == &o) return *this;
        std::lock_guard<std::mutex> const foo(o.mutex());
        std::lock_guard<std::mutex> const bar(mutex());
        (T&)*this = (T const&) o;
        return *this;
    }
    lock_guarded_container& operator=(lock_guarded_container<T> && o) noexcept {
        std::lock_guard<std::mutex> const foo(o.mutex());
        std::lock_guard<std::mutex> const bar(mutex());
        std::swap((T&)*this, (T&) o);
        return *this;
    }
    ~lock_guarded_container() = default;
};

#endif	/* CADO_LOCK_GUARDED_CONTAINER_HPP */
