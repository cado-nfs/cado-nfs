#ifndef MULTITYPED_ARRAY_HPP_
#define MULTITYPED_ARRAY_HPP_

#include <utility>

/* A multityped_array<F, 1, 4> is equivalent to
 *   struct foo {
 *      F<1>::type x1;
 *      F<2>::type x2;
 *      F<3>::type x3;
 *      F<4>::type x4;
 *   };
 * where the xi's are expository only.
 *
 * To get variable x2, from an object X of type multityped_array<F, n> (for
 * any n>=2), do X.get<2>(), and store that in a value of type
 * F<2>::type (or reference, or const reference).
 */

template<template<int> class F, int n0, int n1> struct multityped_array;

namespace multityped_array_details {

    template<typename T, int depth> struct dig : public dig<typename T::super, depth-1> {};
    template<typename T> struct dig<T,0> {
        typedef typename T::type type;
        static inline type & get(T & A) { return A.x; }
        static inline type const & get(T const & A) { return A.x; }
    };
}

template<template<int> class F, int n0, int n1> struct multityped_array : public multityped_array<F, n0, n1-1> {
    typedef multityped_array<F, n0, n1> self;
    typedef multityped_array<F, n0, n1-1> super;
    typedef typename F<n1-1>::type type;
    typename F<n1-1>::type x;
    /* We could use multityped_array_details::dig<self, n1-1-k>::type instead of
     * F<k>::type, but that would ruin the nice diagnostic that
     * static_assert offers here
     */
    template<int k> inline typename F<k>::type& get() {
        static_assert(n0 <= k && k < n1, "attempt to get member of multityped_array out of bounds");
        return multityped_array_details::dig<self, n1-1-k>::get(*this);
    }
    template<int k> inline typename F<k>::type const & get() const {
        static_assert(n0 <= k && k < n1, "attempt to get member of multityped_array out of bounds");
        return multityped_array_details::dig<self, n1-1-k>::get(*this);
    }
    template<typename... Args> multityped_array(Args&&... args)
        : super(std::forward<Args>(args)...)
        , x(std::forward<Args>(args)...)
    {}
};
template<template<int> class F, int n0> struct multityped_array<F, n0, n0> {
    template<typename... Args> multityped_array(Args&&...) {}
};

/* Here's an example of how we can do a for-loop on this sort of things
 *
 * multityped_array_foreach<G> on a multityped_array calls G<i> on the
 * data member of index i (which has type F<i>::type, where F is the one
 * that occurs in the definition of the multityped_array).
 */

#if 0
template<template<int> class G> class multityped_array_foreach {
    template<template<int> class F, typename... Args> struct inner {
        template<int n0>
        void operator()(multityped_array<F, n0, n0> &, Args...) const {}
        template<int n0>
        void operator()(multityped_array<F, n0, n0> const &, Args...) const {}

        template<int n0, int n1>
        void operator()(multityped_array<F, n0, n1> & A, Args... args) const {
            operator()((multityped_array<F, n0, n1-1> &) A, args...);
            G<n1-1>()(A.x, args...);
        }
        template<int n0, int n1>
        void operator()(multityped_array<F, n0, n1> const & A, Args... args) const {
            operator()((multityped_array<F, n0, n1-1> const &) A, args...);
            G<n1-1>()(A.x, args...);
        }
    };
    public:
    template<template<int> class F, int n0, int n1, typename... Args>
        void operator()(multityped_array<F, n0, n1> & A, Args... args) const {
            return inner<F, Args...>()(A, args...);
        }
    template<template<int> class F, int n0, int n1, typename... Args>
        void operator()(multityped_array<F, n0, n1> const & A, Args... args) const {
            return inner<F, Args...>()(A, args...);
        }
};
#endif

template<typename G, template<int> class F> struct multityped_array_foreach_inner {
    template<int n0>
    void operator()(G &, multityped_array<F, n0, n0> &) const {}
    template<int n0>
    void operator()(G &, multityped_array<F, n0, n0> const &) const {}

    template<int n0, int n1>
    void operator()(G & g, multityped_array<F, n0, n1> & A) const {
        operator()(g, (multityped_array<F, n0, n1-1> &) A);
        // g.template operator()<n1-1>(A.x);
        g(A.x);
    }
    template<int n0, int n1>
    void operator()(G & g, multityped_array<F, n0, n1> const & A) const {
        operator()(g, (multityped_array<F, n0, n1-1> const &) A);
        // g.template operator()<n1-1>(A.x);
        g(A.x);
    }
};

template<typename G, template<int> class F, int n0, int n1>
void multityped_array_foreach(G & g, multityped_array<F, n0, n1> & A) {
    return multityped_array_foreach_inner<G, F>()(g, A);
}
template<typename G, template<int> class F, int n0, int n1>
void multityped_array_foreach(G & g, multityped_array<F, n0, n1> const & A) {
    return multityped_array_foreach_inner<G, F>()(g, A);
}
/*
 * old g++ seems to have difficulties with this variant, and is puzzled
 * by the apparent ambiguity -- newer g++ groks it correctly, as does
 * clang
template<typename G, template<int> class F, int n0, int n1>
void multityped_array_foreach(G && g, multityped_array<F, n0, n1> & A) {
    return multityped_array_foreach_inner<G, F>()(g, A);
}
template<typename G, template<int> class F, int n0, int n1>
void multityped_array_foreach(G && g, multityped_array<F, n0, n1> const & A) {
    return multityped_array_foreach_inner<G, F>()(g, A);
}
*/


template<typename G, typename T, template<int> class F> struct multityped_array_fold_inner {
    template<int n0>
    T operator()(G &, T const & t0, multityped_array<F, n0, n0> &) const { return t0; }
    template<int n0>
    T operator()(G &, T const & t0, multityped_array<F, n0, n0> const &) const { return t0; }

    template<int n0, int n1>
    T operator()(G & g, T const & t0, multityped_array<F, n0, n1> & A) const {
        // T t1 = g.template operator()<n1-1>(t0, A.x);
        T t1 = g(t0, A.x);
        return operator()(g, t1, (multityped_array<F, n0, n1-1> &) A);
    }
    template<int n0, int n1>
    T operator()(G & g, T const & t0, multityped_array<F, n0, n1> const & A) const {
        // T t1 = g.template operator()<n1-1>(t0, A.x);
        T t1 = g(t0, A.x);
        return operator()(g, t1, (multityped_array<F, n0, n1-1> &) A);
    }
};

template<typename G, typename T, template<int> class F, int n0, int n1>
T multityped_array_fold(G & g, T const &t0, multityped_array<F, n0, n1> & A) {
    return multityped_array_fold_inner<G, T, F>()(g, t0, A);
}
template<typename G, typename T, template<int> class F, int n0, int n1>
T multityped_array_fold(G & g, T const &t0, multityped_array<F, n0, n1> const & A) {
    return multityped_array_fold_inner<G, T, F>()(g, t0, A);
}
/*
 * old g++ seems to have difficulties with this variant, and is puzzled
 * by the apparent ambiguity -- newer g++ groks it correctly, as does
 * clang
template<typename G, typename T, template<int> class F, int n0, int n1>
T multityped_array_fold(G && g, T const &t0, multityped_array<F, n0, n1> & A) {
    return multityped_array_fold_inner<G, T, F>()(g, t0, A);
}
template<typename G, typename T, template<int> class F, int n0, int n1>
T multityped_array_fold(G && g, T const &t0, multityped_array<F, n0, n1> const & A) {
    return multityped_array_fold_inner<G, T, F>()(g, t0, A);
}
*/


/* This one is overly complicated */
template<typename G> class multityped_array_locate {
    template<template<int> class F> struct inner {
        template<int n0>
        typename G::type operator()(multityped_array<F, n0, n0> &, typename G::key_type &) const { return typename G::type(); }
        template<int n0>
        typename G::type operator()(multityped_array<F, n0, n0> const &, typename G::key_type &) const { return typename G::type(); }

        template<int n0, int n1>
        typename G::type operator()(multityped_array<F, n0, n1> & A, typename G::key_type & k) const {
            typename G::type res = (*this)((multityped_array<F, n0, n1-1> &) A, k);
            if (res != typename G::type())
                return res;
            // return G().template operator()<n1-1>(A.x, k);
            return G()(A.x, k);
        }
        template<int n0, int n1>
        typename G::type operator()(multityped_array<F, n0, n1> const & A, typename G::key_type & k) const {
            typename G::type res = (*this)((multityped_array<F, n0, n1-1> &) A, k);
            if (res != typename G::type())
                return res;
            // return G().template operator()<n1-1>(A.x, k);
            return G()(A.x, k);
        }
    };
    public:
    /* The key is mutable through the process for the inner functions.
     * However we do not wish to expose that to the caller.
     */
    template<template<int> class F, int n0, int n1>
        typename G::type operator()(multityped_array<F, n0, n1> & A, typename G::key_type k) const {
            return inner<F>()(A, k);
        }
    template<template<int> class F, int n0, int n1>
        typename G::type operator()(multityped_array<F, n0, n1> const & A, typename G::key_type k) const {
            return inner<F>()(A, k);
        }
    /*
    template<template<int> class F, int n0, int n1>
        typename G::type operator()(multityped_array<F, n0, n1> & A, typename G::key_type & k) const {
            return inner<F>()(A, k);
        }
    template<template<int> class F, int n0, int n1>
        typename G::type operator()(multityped_array<F, n0, n1> const & A, typename G::key_type & k) const {
            return inner<F>()(A, k);
        }
        */
};

/*
 *
 * Here is one example of how we can use this construction.
 *

#include <vector>
#include <array>
#include <iostream>
#include "utils/multityped_array.hpp"

template<int n> struct type_factory {
    typedef std::array<int, n> type;
};

struct print {
    template<typename T>
    void operator()(T const & x) {
        constexpr size_t n = x.size();
        std::cout << "field x" << n << " has elements";
        for(auto a : x) std::cout << " " << a;
        std::cout << std::endl;
    }
};

struct print2 {
    int k;
    template<typename T>
    void operator()(T const & x) {
        constexpr size_t n = x.size();
        std::cout << "field x" << n << " has elements";
        for(auto a : x) std::cout << " " << a;
        std::cout << " [" << k << "]";
        std::cout << std::endl;
    }
};

struct fill {
    template<typename T>
    void operator()(T & x) {
        for(unsigned int i = 0 ; i < x.size() ; ++i)
            x[i] = i*x.size();
    }
};

struct return_pointer_if_in_subrange {
    typedef int * type;
    typedef int key_type;
    template<typename T>
    type operator()(T & x, int & k) {
        if ((size_t) k < x.size()) {
            return &(x[k]);
        } else {
            k -= x.size();
            return NULL;
        }
    }
    template<typename T>
    type operator()(T const & x, int & k) {
        if ((size_t) k < x.size()) {
            return &(x[k]);
        } else {
            k -= x.size();
            return NULL;
        }
    }
};

struct accumulate_sizes {
    template<typename T>
    int operator()(int t0, T const & a) const {
        return t0 + a.size();
    }
};



int main()
{
    typedef multityped_array<type_factory, 0, 6> A_t;
    A_t A;

    multityped_array_foreach(fill(), A);

    type_factory<3>::type A3 = A.get<3>();
    std::cout << "A3[2] == " << A3[2] << std::endl;

    multityped_array_foreach(print2 {(int) 1}, A);

    int v = 5;
    multityped_array_locate<return_pointer_if_in_subrange>()(A, v);

    int s = multityped_array_fold(accumulate_sizes(), 0, A);

    std::cout << "total " << s << std::endl;
}

 *
 */

#endif	/* MULTITYPED_ARRAY_HPP_ */
