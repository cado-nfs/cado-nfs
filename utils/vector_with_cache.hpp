#ifndef UTILS_VECTOR_WITH_CACHE_HPP_
#define UTILS_VECTOR_WITH_CACHE_HPP_

#include <cstddef>
#include <vector>
#include <array>
#include <type_traits>

#include "macros.h"

/* This is a (very partial) equivalent of std::vector that provides
 * no-indirection storage for its element when it so happens that only a
 * few such elements need to be stored.
 *
 * It's eminently unclear that it actually provides a measurable speedup.
 * The only point is that our old code had this feature baked in, so
 * let's provide an abstraction layer for it before we decide if it is
 * worth keeping or not.
 *
 * Many things aren't implemented.
 */
template<typename T, size_t N = 16>
requires std::is_move_assignable_v<T>
class vector_with_cache
{
    std::vector<T> heap;
    std::array<T, N> immediate = {};
    
    /* an ad hoc implementation that reinvents std::vector would not have
     * to redefine these fields, as it would be possible to subvert the
     * ones from the vector code in order to achieve this goal.
     */
    T * data = nullptr;
    size_t size_ = 0;

    public:

    vector_with_cache()
        : data(immediate.data())
    {}

    size_t size() const { return size_; }
    T & operator[](size_t i) { return data[i]; }
    T const & operator[](size_t i) const { return data[i]; }
    T & last() { return data[size_-1]; }
    T const & last() const { return data[size_-1]; }

    ~vector_with_cache() = default;
    vector_with_cache(vector_with_cache<T, N> const & o)
        : size_(o.size_)
    {
        if (size_ <= N) {
            immediate = o.immediate;
            data = immediate.data();
        } else {
            heap = o.heap;
            data = heap.data();
        }
    }
    vector_with_cache(vector_with_cache<T, N> && o) noexcept
        : size_(o.size_)
    {
        if (size_ <= N) {
            immediate = o.immediate;
            data = immediate.data();
        } else {
            heap = std::move(o.heap);
            data = heap.data();
        }
        o.clear();
    }
    vector_with_cache& operator=(vector_with_cache<T, N> const & o)
    {
        if (this == &o)
            return *this;
        size_ = o.size_;
        if (size_ <= N) {
            immediate = o.immediate;
            data = immediate.data();
        } else {
            heap = o.heap;
            data = heap.data();
        }
        return *this;
    }
    vector_with_cache& operator=(vector_with_cache<T, N> && o) noexcept
    {
        if (this == &o)
            return *this;
        size_ = o.size_;
        if (size_ <= N) {
            immediate = o.immediate;
            data = immediate.data();
        } else {
            heap = std::move(o.heap);
            data = heap.data();
        }
        o.clear();
        return *this;
    }

    using iterator = T *;
    using const_iterator = T const *;

    iterator begin() { return data; }
    iterator end() { return data + size_; }
    const_iterator begin() const { return data; }
    const_iterator end() const { return data + size_; }
    bool empty() const { return size_ == 0; }
    void clear() { size_ = 0; heap.clear(); data = immediate.data(); }
    iterator erase(iterator b, iterator e) {
        auto pos = b - begin();
        if (data == immediate.data()) {
            /* easy! */
            std::move(e, end(), b);
            size_ -= e - b;
        } else {
            heap.erase(heap.begin() + pos, heap.begin() + pos + (e - b));
            size_ = heap.size();
            if (heap.size() <= N) {
                std::move(heap.begin(), heap.end(), immediate.begin());
                data = immediate.data();
            }
        }
        return data + pos;
    }

    private:
    void transfer_before_extend() {
        if (size_ == N) {
            ASSERT(data == immediate.data());
            heap.clear();
            heap.reserve(N);
            for(auto & c : immediate)
                heap.push_back(std::move(c));
            data = heap.data();
        }
    }
    public:
    void push_back(T const & a) {
        transfer_before_extend();
        if (size_ >= N) {
            ASSERT(data == heap.data());
            heap.push_back(a);
            /* make sure that reallocs don't kill our pointer */
            data = heap.data();
            size_++;
        } else {
            ASSERT(data == immediate.data());
            immediate[size_++] = a;
        }
    }
    template<typename... Args>
    void emplace_back(Args&& ...args)
    requires std::is_constructible_v<T, Args...>
    {
        transfer_before_extend();
        if (size_ >= N) {
            ASSERT(data == heap.data());
            heap.emplace_back(std::forward<Args>(args)...);
            /* make sure that reallocs don't kill our pointer */
            data = heap.data();
            size_++;
        } else {
            ASSERT(data == immediate.data());
            immediate[size_++] = T(std::forward<Args>(args)...);
        }
    }
};

#endif	/* UTILS_VECTOR_WITH_CACHE_HPP_ */
