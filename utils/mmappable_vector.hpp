#ifndef _MMAPPABLE_VECTOR_HPP_
#define _MMAPPABLE_VECTOR_HPP_

/* This is inspired from https://github.com/johannesthoma/mmap_allocator
 * License is LGPL.
 * The version here has been trimmed down significantly, look for uptream
 * version for more features */

#include <memory>
#include <string>
#include <stdio.h>
#include <vector>
#include "mmap_allocator.hpp"
#include "cxx_std_vector_ugly_accessor.hpp"

template <typename T, typename A = mmap_allocator_details::mmap_allocator<T> >
class mmappable_vector: public std::vector<T, A> {
    public:
        typedef std::vector<T, A> Base;

        typedef typename Base::const_iterator const_iterator;
        typedef typename Base::iterator iterator;
        typedef T value_type;
        typedef A allocator_type;

        mmappable_vector(): Base() { }

        mmappable_vector(const mmappable_vector<T, A> &) = default;
        mmappable_vector& operator=(const mmappable_vector<T, A> &) = default;

        /* This is not conforming, since the container
         * requirements for this signature stipulate that the
         * items must be value-initialized */
        //explicit mmappable_vector(size_t n): Base() { mmap(n); }

        /* prefer this one, which will boom if the mapping is
         * readonly */
        explicit mmappable_vector(size_t n): Base(n) { }

        explicit mmappable_vector(A alloc): Base(alloc) { }

        mmappable_vector(iterator from, iterator to): Base(from, to) { }

        template <typename Iter>
            mmappable_vector(Iter first, Iter last, A a = A()):
                Base(first, last, a)
    { }

        mmappable_vector(int n, T val, A alloc): Base(n, val, alloc) { }

        mmappable_vector(int n, T val): Base(n, val) { }

        mmappable_vector(std::vector<T,std::allocator<T> > v):
            std::vector<T,std::allocator<T> >(v)
    { }


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
            Base::reserve(n);
            auto & u = cado_nfs::ugly_accessor(*this);
            u._finish() = u._start() + n;
        }
        void munmap()
        {
            Base::clear();
            auto & u = cado_nfs::ugly_accessor(*this);
            // The following two lines are typically what the dtor does,
            // but we should rather call the dtor explicitly, since it
            // has access to private stuff that is a better match to our
            // use of reserve() above. In particular, llvm libcxx with
            // address sanitizer mode would want this.
            //
            // size_t n = Base::size();
            // u._allocator().deallocate(u._start(), n);
            reinterpret_cast<Base *>(this)->~Base();
            u._start() = 0;
            u._finish() = 0;
            u._end_of_storage() = 0;
        }

        /* Those two are only to align with the original code I
         * found on github, although personally I don't like
         * this kind of shortcut. These two methods are
         * explicitly tinkering with the allocator field */

        /* Adding enable_if because really that only makes sense
         * with our allocator, no other */
        typename std::enable_if<std::is_same<mmap_allocator_details::mmap_allocator<T>, A>::value>::type mmap_file(const char * filename, mmap_allocator_details::access_mode mode, mmap_allocator_details::offset_type offset, mmap_allocator_details::size_type length) {
            A & a(cado_nfs::ugly_accessor(*this)._allocator());
            if (a.has_defined_mapping()) throw mmap_allocator_details::mmap_allocator_exception("already mapped");
            a = A(filename, mode, offset, length);
            mmap(length);
        }
        typename std::enable_if<std::is_same<mmap_allocator_details::mmap_allocator<T>, A>::value>::type munmap_file() {
            A & a(cado_nfs::ugly_accessor(*this)._allocator());
            if (!a.has_defined_mapping()) throw mmap_allocator_details::mmap_allocator_exception("not yet mapped");
            munmap();
            a = A();
        }


        void swap(mmappable_vector<T,A> & __x) {
            /* /usr/include/c++/7/bits/stl_vector.h:103:
             * std::_Vector_base<T,A>::_M_swap_data does not seem to do
             * anything about the data fields of the allocator. IDK if
             * it's a bug or a feature, but that sounds definitely
             * worrisome.
             */
            auto & u = cado_nfs::ugly_accessor(*this);
            auto & x = cado_nfs::ugly_accessor(__x);
            std::swap(u._allocator(), x._allocator());
            std::swap(u._start(), x._start());
            std::swap(u._finish(), x._finish());
            std::swap(u._end_of_storage(), x._end_of_storage());
        }

};


template <typename T, typename A>
void swap(mmappable_vector<T,A> & a, mmappable_vector<T,A> & b)
{
    a.swap(b);
}

#endif /* MMAPPABLE_VECTOR_HPP_ */
