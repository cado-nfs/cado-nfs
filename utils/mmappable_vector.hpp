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

template <typename T, typename A = mmap_allocator_details::mmap_allocator<T> >
class mmappable_vector: public std::vector<T, A> {
    public:
        typedef std::vector<T, A> Base;

        typedef typename Base::const_iterator const_iterator;
        typedef typename Base::iterator iterator;
        typedef T value_type;
        typedef A allocator_type;

        mmappable_vector(): Base() { }

        /* I think this does plain crap */
        mmappable_vector(const mmappable_vector<T, A> &other): Base(other) { } 

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
#ifdef HAVE_GLIBC_VECTOR_INTERNALS
            Base::_M_impl._M_finish = Base::_M_impl._M_start + n;
#else
#error "Not GNU libstdc++, please expand code"
#endif
        }
        void munmap()
        {
            size_t n = Base::size();
            Base::clear();
#ifdef HAVE_GLIBC_VECTOR_INTERNALS
            Base::_M_deallocate(Base::_M_impl._M_start, n);
            Base::_M_impl._M_start = 0;
            Base::_M_impl._M_finish = 0;
            Base::_M_impl._M_end_of_storage = 0;
#else
#error "Not GNU libstdc++, please expand code"
#endif
        }

        /* Those two are only to align with the original code I
         * found on github, although personally I don't like
         * this kind of shortcut. These two methods are
         * explicitly tinkering with the allocator field */

        /* Adding enable_if because really that only makes sense
         * with our allocator, no other */
        typename std::enable_if<std::is_same<mmap_allocator_details::mmap_allocator<T>, A>::value>::type mmap_file(const char * filename, mmap_allocator_details::access_mode mode, mmap_allocator_details::offset_type offset, mmap_allocator_details::size_type length) {
#ifdef HAVE_GLIBC_VECTOR_INTERNALS
            A & a(Base::_M_get_Tp_allocator());
#else
#error "Not GNU libstdc++, please expand code"
#endif
            if (a.has_defined_mapping()) throw mmap_allocator_details::mmap_allocator_exception("already mapped");
            a = A(filename, mode, offset, length);
            mmap(length);
        }
        typename std::enable_if<std::is_same<mmap_allocator_details::mmap_allocator<T>, A>::value>::type munmap_file() {
#ifdef HAVE_GLIBC_VECTOR_INTERNALS
            A & a(Base::_M_get_Tp_allocator());
#else
#error "Not GNU libstdc++, please expand code"
#endif
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
#ifdef HAVE_GLIBC_VECTOR_INTERNALS
            A & a(Base::_M_get_Tp_allocator());
            std::swap(a, (A&) __x._M_get_Tp_allocator());
            std::swap(Base::_M_impl._M_start, __x._M_impl._M_start);
            std::swap(Base::_M_impl._M_finish, __x._M_impl._M_finish);
            std::swap(Base::_M_impl._M_end_of_storage, __x._M_impl._M_end_of_storage);
#else
#error "Not GNU libstdc++, please expand code"
#endif
        }

};


template <typename T, typename A>
void swap(mmappable_vector<T,A> & a, mmappable_vector<T,A> & b)
{
    a.swap(b);
}

#endif /* MMAPPABLE_VECTOR_HPP_ */
