#ifndef CADO_MULTITYPED_ARRAY_HPP
#define CADO_MULTITYPED_ARRAY_HPP

#include <cstddef>

#include <tuple>
#include <type_traits>
#include <utility>

/* A cado::multityped_array<F, 1, 5> is equivalent to
 *   struct foo {
 *      F<1>::type x1;
 *      F<2>::type x2;
 *      F<3>::type x3;
 *      F<4>::type x4;
 *   };
 * where the xi's are expository only.
 *
 * To get variable x2, from an object X of type cado::multityped_array<F,
 * n0, n1> (assuming n0, n1 are such that n0 <= 2 < n1), do X.get<2>(),
 * and store that in a value of type F<2>::type (or reference, or const
 * reference).
 *
 * To iterate a function f on all the xis, do cado::foreach(X, f)
 *
 * To get a pointer to one of the xis such that a predicate f(xi) is
 * true, do cado::find(X, f), and store the result in a value of type
 * std::common_type<F<n0>::type, ...,  F<n1>::type>::type (which is also
 * decltype(X)::common_type). This requires that
 * proper overloads be defined, such as for example
 * template<size_t m, size_t n> struct std::common_type<F<m>, F<n>> {
 * using type = F_base; }; (assuming F<m>::type is F<m> -- otherwise this
 * construct must be adapted)
 */

namespace cado
{

/* as an exercise, here is the multityped_array construction,
 * implemented in a simpler way with this proxy
 */
namespace multityped_array_details
{

/* the multityped_array_impl retains information on the starting
 * index n0, and has info on the differences relative to n0. Its
 * default implementation is nil.
 * TODO: isn't it odd that we don't even define it? Sounds like an
 * incomplete type, right?
 */
template <template <int> class F, int n0,
          typename T = std::integer_sequence<int>>
struct multityped_array_impl;

template<typename T, typename = void>
struct has_type_member : std::false_type {};

template<typename T>
struct has_type_member<T, std::void_t<typename T::type>>
: std::true_type {};
/*
*/

template <template <int> class F, int n0, int... i>
struct multityped_array_impl<F, n0, std::integer_sequence<int, i...>> {
    /* The type expands to a full tuple with all types created by the
     * type functor F
     */
    using type = std::tuple<typename F<n0 + i>::type...>;

    using common_type = std::common_type<typename F<n0 + i>::type...>;

    template<typename Callable>
    using merged_type_on_ref = std::common_type<std::invoke_result_t<Callable, typename F<n0 + i>::type &>...>;
    template<typename Callable>
    using merged_type_on_cref = std::common_type<std::invoke_result_t<Callable, typename F<n0 + i>::type const &>...>;
};

template <template <int> class F, int n0, int n1>
using multityped_array_helper =
    multityped_array_impl<F, n0, std::make_integer_sequence<int, n1 - n0>>;

template <typename Search, typename S = std::integer_sequence<int>>
struct find_helper;



template <template <int> class F, int n0, int n1>
struct multityped_array : multityped_array_helper<F, n0, n1>::type {

#if 0
    using multityped_array_helper<F, n0, n1>::merged_type_on_cref;
    using multityped_array_helper<F, n0, n1>::merged_type_on_ref;
#else
    template<typename Callable>
        using merged_type_on_cref = multityped_array_helper<F, n0, n1>::template merged_type_on_cref<Callable>;
    template<typename Callable>
        using merged_type_on_ref = multityped_array_helper<F, n0, n1>::template merged_type_on_ref<Callable>;
#endif

    using common_type = multityped_array_helper<F, n0, n1>::common_type;
    /* We provide member function helpers for get. However, it makes
     * sense to prefer cado::get<> helpers, which do not need the
     * extra .template syntax!
     *
     * A very important difference with std::get (that works on a
     * tuple) is that here our indices range from n0 to n1, instead
     * of from 0 to n1 - n0.
     *
     * Likewise, we define tuple_element member types that are
     * indexed from n0 to n1
     */
    template <int n> typename F<n>::type & get()
    {
        static_assert(n0 <= n && n < n1,
                      "Cannot call n for a member index out of range");
        return std::get<n - n0>(*this);
    }
    template <int n> typename F<n>::type const & get() const
    {
        static_assert(n0 <= n && n <= n1,
                      "Cannot call n for a member index out of range");
        return std::get<n - n0>(*this);
    }

    template <int n> using tuple_element = F<n>::type;

    template<typename Callable> inline void foreach(Callable&&);
    template<typename Callable> inline void foreach(Callable&&) const;

    template<typename Callable>
    auto
    find(Callable && f)
    requires has_type_member<merged_type_on_ref<Callable>>::value
    {
        return multityped_array_details::find_helper<
            Callable, typename std::make_integer_sequence<int, n1 - n0>>()(
                    std::forward<Callable>(f), *this);
    }

    template<typename Callable>
    auto
    find(Callable && f) const
    requires has_type_member<merged_type_on_cref<Callable>>::value
    {
        return multityped_array_details::find_helper<
            Callable, typename std::make_integer_sequence<int, n1 - n0>>()(
                    std::forward<Callable>(f), *this);
    }

};

/* foreach().
 * f must be callable with all elements of x in turn.
 *
 * We create a function object which recursively calls its lower-order
 * siblings until it reaches the bottom.
 */
template <typename Callable, typename S = std::integer_sequence<int>>
struct foreach_helper {
    template <typename T> void operator()(Callable &&, T &) {}
    template <typename T> void operator()(Callable &&, T const &) {}
};
template <typename Callable, int n0, int... indices>
struct foreach_helper<Callable, std::integer_sequence<int, n0, indices...>> {
    template <typename T> void operator()(Callable && f, T & x) const
    {
        f(std::get<n0>(x));
        foreach_helper<Callable, std::integer_sequence<int, indices...>>()(
            std::forward<Callable>(f), x);
    }
    template <typename T> void operator()(Callable && f, T const & x) const
    {
        f(std::get<n0>(x));
        foreach_helper<Callable, std::integer_sequence<int, indices...>>()(
            std::forward<Callable>(f), x);
    }
};
template<template <int> class F, int n0, int n1>
template <typename Callable>
inline void
multityped_array<F, n0, n1>::foreach(Callable && f)
{
    multityped_array_details::foreach_helper<
        Callable, typename std::make_integer_sequence<int, n1 - n0>>()(
                std::forward<Callable>(f), *this);
}
template<template <int> class F, int n0, int n1>
template <typename Callable>
inline void
multityped_array<F, n0, n1>::foreach(Callable && f) const
{
    multityped_array_details::foreach_helper<
        Callable, typename std::make_integer_sequence<int, n1 - n0>>()(
                std::forward<Callable>(f), *this);
}

/* find0().
 * f must be a predicate that can be evaluated on all elements of
 * x.
 * We return a pointer to the first
 * member that matches the predicate, or nullptr if there is none.
 */
template <typename Predicate, typename S>
struct find0_helper {
    template <typename T> std::nullptr_t operator()(Predicate &&, T const &)
    {
        return nullptr;
    }
};
template <typename Predicate, int n0, int... indices>
struct find0_helper<Predicate, std::integer_sequence<int, n0, indices...>> {
    template <typename T>
    typename T::common_type::type * operator()(Predicate && f, T & x) const
    {
        if (f(std::get<n0>(x)))
            return &std::get<n0>(x);
        else
            return find0_helper<Predicate,
                               std::integer_sequence<int, indices...>>()(
                std::forward<Predicate>(f), x);
    }
    template <typename T>
    typename T::common_type::type * operator()(Predicate && f, T const & x) const
    {
        if (f(std::get<n0>(x)))
            return &std::get<n0>(x);
        else
            return find0_helper<Predicate,
                               std::integer_sequence<int, indices...>>()(
                std::forward<Predicate>(f), x);
    }
};


/* find().
 */
template <typename Caller, int n0>
struct find_helper<Caller, std::integer_sequence<int, n0>> {
    template <typename T>
        auto operator()(Caller && f, T const & x)
    {
        return f(std::get<n0>(x));
    }
};
template <typename Caller, int n0, int... indices>
struct find_helper<Caller, std::integer_sequence<int, n0, indices...>> {
    template <typename T>
    typename T::template merged_type_on_ref<Caller>::type operator()(Caller && f, T & x) const
    {
        auto c = f(std::get<n0>(x));
        if (c)
            return c;
        return find_helper<Caller,
                               std::integer_sequence<int, indices...>>()(
                std::forward<Caller>(f), x);
    }
    template <typename T>
    typename T::template merged_type_on_cref<Caller>::type operator()(Caller && f, T const & x) const
    {
        auto c = f(std::get<n0>(x));
        if (c)
            return c;
        return find_helper<Caller,
                               std::integer_sequence<int, indices...>>()(
                std::forward<Caller>(f), x);
    }
};
} /* namespace multityped_array_details */

using multityped_array_details::multityped_array;

/* get() accessors */
template <int n, template <int> class F, int n0, int n1>
typename F<n>::type & get(multityped_array<F, n0, n1> & x)
{
    return x.template get<n>();
}
template <int n, template <int> class F, int n0, int n1>
typename F<n>::type const & get(multityped_array<F, n0, n1> const & x)
{
    return x.template get<n>();
}


} /* namespace cado */
#endif /* CADO_MULTITYPED_ARRAY_HPP */
