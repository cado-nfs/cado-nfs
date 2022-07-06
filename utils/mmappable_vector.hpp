#ifndef _MMAPPABLE_VECTOR_HPP_
#define _MMAPPABLE_VECTOR_HPP_

#include "mmap_allocator.hpp"
#include <algorithm>
#include <iterator>
#include <limits>
#include <memory>
#include <stdexcept>
#include <stdio.h>
#include <string>
#include <type_traits>
#include <vector>

/* see bdfb2543f for a version that plays ugly tricks with the
 * std::vector internal fields. The version below doesn't do such things.
 */

template<typename T>
struct works_with_mmappable_vector : public std::is_trivially_constructible<T>
{};

template<typename T, typename A = mmap_allocator_details::mmap_allocator<T>>
class mmappable_vector : public A
{
    typedef mmappable_vector<T, A> self;
    static_assert(
      works_with_mmappable_vector<T>::value,
      "mmap_allocator only works with types that are explicitly allowed");

    public:
    typedef T value_type;
    typedef A allocator_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef typename std::allocator_traits<allocator_type>::pointer pointer;
    typedef typename std::allocator_traits<allocator_type>::const_pointer
      const_pointer;
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
    typedef
      typename std::iterator_traits<iterator>::difference_type difference_type;
    typedef typename std::make_unsigned<difference_type>::type size_type;

    private:
    pointer _start;
    pointer _finish;
    pointer _end_of_storage;

    public:
    mmappable_vector()
      : _start(nullptr)
      , _finish(nullptr)
      , _end_of_storage(nullptr)
    {}

    explicit mmappable_vector(A const& alloc)
      : A(alloc)
      , _start(nullptr)
      , _finish(nullptr)
      , _end_of_storage(nullptr)
    {}

    mmappable_vector(size_type n, T const& val, A const& alloc = A())
      : A(alloc)
      , _start(A::allocate(n))
      , _finish(_start + n)
      , _end_of_storage(_finish)
    {
        std::fill_n(_start, n, val);
    }

    explicit mmappable_vector(size_type n)
      : mmappable_vector(n, T())
    {}

    mmappable_vector(const mmappable_vector<T, A>& o, A const& a = A())
      : mmappable_vector(o.begin(), o.end(), a)
    {}

    self& operator=(self const& x)
    {
        if (&x == this)
            return *this;
        self t(x, get_allocator());
        swap(t);
    };

    mmappable_vector(mmappable_vector<T, A>&& o)
      : _start(o._start)
      , _finish(o._finish)
      , _end_of_storage(o._end_of_storage)
    {
        get_allocator_ref() = o;
        o._start = nullptr;
        o._finish = nullptr;
        o._end_of_storage = nullptr;
    }

    mmappable_vector& operator=(mmappable_vector<T, A>&& o) {
        swap(o);
        return *this;
    }

    template<typename Iter,
             typename = typename std::is_convertible<
               typename std::iterator_traits<Iter>::value_type,
               value_type>::type>
    mmappable_vector(Iter first, Iter last, A const& a = A())
      : A(a)
    {
        range_initialize(
          first,
          last,
          typename std::iterator_traits<Iter>::iterator_category());
    }

    private:
    template<typename Iter>
    void range_initialize(Iter first,
                          Iter last,
                          std::random_access_iterator_tag)
    {
        typename std::iterator_traits<Iter>::difference_type n = std::distance(first, last);
        _start = A::allocate(n);
        _finish = _start + n;
        _end_of_storage = _finish;
        std::copy(first, last, _start);
    }

    template<typename Iter>
    void range_initialize(Iter first, Iter last, std::forward_iterator_tag)
    {
        for (; first != last; push_back(*first++))
            ;
    }

    public:
    ~mmappable_vector() { A::deallocate(_start, capacity()); }

    inline size_type capacity() const { return _end_of_storage - _start; }
    inline size_type size() const { return _finish - _start; }

    inline iterator begin() { return _start; }
    inline iterator end() { return _finish; }
    inline reverse_iterator rbegin() { return reverse_iterator(_finish); }
    inline reverse_iterator rend() { return reverse_iterator(_start); }
    inline const_iterator cbegin() const { return _start; }
    inline const_iterator cend() const { return _finish; }
    inline const_reverse_iterator crbegin() const
    {
        return reverse_iterator(_finish);
    }
    inline const_reverse_iterator crend() const
    {
        return reverse_iterator(_start);
    }
    inline const_iterator begin() const { return _start; }
    inline const_iterator end() const { return _finish; }
    inline const_reverse_iterator rbegin() const
    {
        return reverse_iterator(_finish);
    }
    inline const_reverse_iterator rend() const
    {
        return reverse_iterator(_start);
    }

    static inline size_type max_size()
    {
        return std::numeric_limits<size_type>::max();
    }

    inline bool empty() const { return _finish == _start; }

    /* missing:
     *
       void resize(size_type) {}
       void shrink_to_fit() {}
       */

    void reserve(size_type n)
    {
        if (n > capacity()) {
            const size_type s = size();
            pointer t = A::allocate(n);
            std::move(_start, _finish, t);
            A::deallocate(_start, capacity());
            _start = t;
            _finish = t + s;
            _end_of_storage = t + n;
        }
    }
    inline reference operator[](size_type n) { return _start[n]; }
    inline const_reference operator[](size_type n) const { return _start[n]; }
    inline reference at(size_type n)
    {
        if (n >= size())
            throw std::out_of_range("mmappable_vector::at");
        return (*this)[n];
    }
    inline const_reference at(size_type n) const
    {
        if (n >= size())
            throw std::out_of_range("mmappable_vector::at");
        return (*this)[n];
    }
    inline reference front() { return *_start; }
    inline const_reference front() const { return *_start; }
    inline reference back() { return *rbegin(); }
    inline const_reference back() const { return *rbegin(); }
    inline pointer data() { return _start; }
    inline const_pointer data() const { return _start; }

    inline allocator_type get_allocator() { return (allocator_type) * this; }

    private:
    inline allocator_type& get_allocator_ref()
    {
        return (allocator_type&)*this;
    }

    public:

    /* Unimplemented modifiers:
       assign
       pop_back
       insert
       erase
       emplace
       */

    void push_back(value_type const& x)
    {
        if (_finish == _end_of_storage)
            reserve(1 + 2 * capacity());
        *_finish++ = x;
    }

    void push_back(value_type&& x)
    {
        if (_finish == _end_of_storage)
            reserve(1 + 2 * capacity());
        *_finish++ = std::move(x);
    }

    template<typename... Args>
    void emplace_back(Args&&... args)
    {
        push_back(value_type(std::forward<Args>(args)...));
    }

    void clear() { _finish = _start; }

    public:
    /* The four functions below are the main reason why we
     * really need to subclass the container, and subclassing
     * the allocator is not sufficient: We need the vector to
     * become aware that there are many elements already
     * there, and that goes with tinkering with the inner
     * fields of the vector type. Handing over
     * always the same pointer is not enough.
     */
    void mmap(size_t n)
    {
        reserve(n);
        _finish = _start + n;
    }
    void munmap()
    {
        clear();
        A::deallocate(_start, capacity());
        _start = nullptr;
        _finish = nullptr;
        _end_of_storage = nullptr;
    }

    /* Those two are only to align with the original code I
     * found on github, although personally I don't like
     * this kind of shortcut. These two methods are
     * explicitly tinkering with the allocator field */

    /* Adding enable_if because really that only makes sense
     * with our allocator, no other */
    typename std::enable_if<
      std::is_same<mmap_allocator_details::mmap_allocator<T>, A>::value>::type
    mmap_file(const char* filename,
              mmap_allocator_details::access_mode mode,
              mmap_allocator_details::offset_type offset,
              mmap_allocator_details::size_type length)
    {
        if (get_allocator_ref().has_defined_mapping())
            throw mmap_allocator_details::mmap_allocator_exception(
              "already mapped");
        get_allocator_ref() = A(filename, mode, offset, length);
        mmap(length);
    }
    typename std::enable_if<
      std::is_same<mmap_allocator_details::mmap_allocator<T>, A>::value>::type
    munmap_file()
    {
        if (!get_allocator_ref().has_defined_mapping())
            throw mmap_allocator_details::mmap_allocator_exception(
              "not yet mapped");
        munmap();
        get_allocator_ref() = A();
    }

    void swap(mmappable_vector<T, A>& x)
    {
        std::swap((allocator_type&)*this, (allocator_type&)x);
        std::swap(_start, x._start);
        std::swap(_finish, x._finish);
        std::swap(_end_of_storage, x._end_of_storage);
    }
};

template<typename T, typename A>
void
swap(mmappable_vector<T, A>& a, mmappable_vector<T, A>& b)
{
    a.swap(b);
}

#endif /* MMAPPABLE_VECTOR_HPP_ */
