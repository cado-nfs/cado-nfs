#ifndef CXX_STD_VECTOR_UGLY_ACCESSOR_HPP_
#define CXX_STD_VECTOR_UGLY_ACCESSOR_HPP_

/* This is ugly. We're creating an accessor that reaches the internal
 * fields of std::vector, so that we can play tricks with the underlying
 * allocator, and use an mmap-backed storage instead. There is no
 * explicit provision for doing that in the standard, so we must resort
 * to unholy stuff.
 */

#include <vector>

namespace cado_nfs {
#ifdef __GLIBCXX__
/* This code works with the gnu libstdc++ */

#ifdef __OpenBSD__
/* openbsd seems to have sort of glibc++, but mmap behaviour is odd. Don't
 * want to bother. (542bc4bce)
 */
#error "forcibly disabling mmappable vectors for openbsd"
#endif

template<typename T, typename Alloc>
struct ugly_accessor_helper : public std::vector<T, Alloc> {
    typedef std::vector<T, Alloc> Base;
    T * & _start() { return Base::_M_impl._M_start; }
    T * & _finish() { return Base::_M_impl._M_finish; }
    T * & _end_of_storage() { return Base::_M_impl._M_end_of_storage; }
    Alloc & _allocator() { return Base::_M_get_Tp_allocator(); };
};

template<typename T, typename Alloc>
ugly_accessor_helper<T, Alloc> & 
ugly_accessor(std::vector<T, Alloc> & c) { return reinterpret_cast<ugly_accessor_helper<T, Alloc> &>(c); }
#endif

#ifdef _LIBCPP_STD_VER
/* This code only works for the llvm libcxx implementation */

template<typename T, typename Alloc>
struct ugly_accessor_helper : public std::__vector_base<T, Alloc> {
    typedef std::__vector_base<T, Alloc> Base;
    T * & _start() { return Base::__begin_; }
    T * & _finish() { return Base::__end_; }
    T * & _end_of_storage() { return Base::__end_cap(); }
    Alloc & _allocator() { return Base::__alloc(); }
};

template<typename T, typename Alloc>
ugly_accessor_helper<T, Alloc> & 
ugly_accessor(std::vector<T, Alloc> & c) {
    typedef std::vector<T, Alloc> Front;
    typedef std::__vector_base<T, Alloc> Base;
    static_assert(std::is_base_of<Base, Front>::value, "This code only works for the llvm libcxx implementation");
    return reinterpret_cast<ugly_accessor_helper<T, Alloc> &>(c);
}
#endif
}



#endif	/* CXX_STD_VECTOR_UGLY_ACCESSOR_HPP_ */
