#ifndef UTILS_CXX_HPP_
#define UTILS_CXX_HPP_

#include <cstdio>
#include <cstdlib>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

#include "fmt/format.h"

#include "macros.h"

/* Base class with private copy-constructor and assignment operator.
   Classes which are not copy-constructible can inherit this with:
   private NonCopyable.
   
   See also:

   https://www.boost.org/doc/libs/1_67_0/libs/core/doc/html/core/noncopyable.html
   */
struct NonCopyable {
   NonCopyable() = default;
   ~NonCopyable() = default;
   NonCopyable(NonCopyable&&) = delete;
   NonCopyable& operator=(NonCopyable&&) = delete;
   NonCopyable(NonCopyable const &) = delete;
   NonCopyable& operator=(NonCopyable const &) = delete;
};

/* This is an example, but not a class one can inherit from: inheriting
 * does not have the same effect as having this kind of ctors declared...
 *
class MoveOnly
{
public:
   MoveOnly() {}
   ~MoveOnly() {}
   MoveOnly(const MoveOnly&) = delete;
   MoveOnly& operator=(const MoveOnly&) = delete;
   MoveOnly(MoveOnly&&) = default;
   MoveOnly& operator=(MoveOnly&&) = default;
};
 */

/* Generic tool to cause an arbitrary closure to be called at destructor
 * time. See detached_cofac for a use case.
 */
template<typename T> struct call_dtor_s {
    T x;
    explicit call_dtor_s(T x): x(x) {}
    ~call_dtor_s() { x(); }
    call_dtor_s& operator=(call_dtor_s&&) = delete;
    call_dtor_s(call_dtor_s const &) = delete;
    call_dtor_s& operator=(call_dtor_s const &) = delete;

    // move works
    call_dtor_s(call_dtor_s&&) noexcept = default;
};
template<typename T> call_dtor_s<T> call_dtor(T x) { return call_dtor_s<T>(x); }

/* miscellaneous stuff for C++ */

/* This is handy for range-based for loops that also need a counter.
 */
template<typename T>
struct increment_counter_on_dtor {
    T & a;
    explicit increment_counter_on_dtor(T & a) : a(a) {}
    ~increment_counter_on_dtor() { ++a; }
    increment_counter_on_dtor(increment_counter_on_dtor&&) = delete;
    increment_counter_on_dtor& operator=(increment_counter_on_dtor&&) = delete;
    increment_counter_on_dtor(increment_counter_on_dtor const &) = delete;
    increment_counter_on_dtor& operator=(increment_counter_on_dtor const &) = delete;
};

/* usage example:
 *
    size_t counter = 0;
    for(auto x : foo) {
        increment_counter_on_dtor<size_t> _dummy(counter);
        if (x % 2 == 0) continue;
        std::cout << "[" << counter << "]: " << x << "\n";
    }

 * of course, &x-begin(foo) is also acceptable, but that works only for
 * random-access iterators.
 */


/* A class that can be used to instrument functions to count events.

The final counts are printed at dtor time; if the counter object is static,
it gets printed at program exit.
Example:

void foo(int x) {
  static StaticHistogram calls(1, "Nr of times foo() was called");
  static StaticHistogram hist(100, "Histogram of values foo() was called with", true);
  calls.inc();
  if (0 <= x && x < 100) hist.inc(x);
  [... rest of foo() ...]
}
*/

class StaticHistogram {
    std::string name;
    const size_t len;
    bool print_indices;
    std::unique_ptr<unsigned long[]> c;
public:
    explicit StaticHistogram(const size_t _len, const char *_name = nullptr, const bool _print_indices = false)
        : name(_name)
        , len(_len)
        , print_indices(_print_indices)
        , c(new unsigned long[_len])
    {
        for (size_t i = 0; i < len; i++)
            c[i] = 0;
    }
    ~StaticHistogram() {
        if (!name.empty())
            printf("%s: ", name.c_str());
        if (print_indices) {
            for (size_t i = 0; i < len; i++)
                if (c[i] > 0)
                    printf("%s%zu:%lu", (i > 0) ? " " : "", i, c[i]);
        } else {
            for (size_t i = 0; i < len; i++)
                printf("%s%lu", (i > 0) ? " " : "", c[i]);
        }
        printf("\n");
    }
    void inc(const size_t i = 0) {
        ASSERT_ALWAYS(i < len);
        c[i]++;
    }
    StaticHistogram(StaticHistogram const &) = delete;
    StaticHistogram& operator=(StaticHistogram const &) = delete;
    StaticHistogram& operator=(StaticHistogram&&) = delete;

    // move is ok.
    StaticHistogram(StaticHistogram&&) = default;
};


#if __cplusplus < 201402L
namespace std {
template <bool B, typename T = void>
using enable_if_t = typename std::enable_if<B, T>::type;
}
#endif

/* A type trait that checks whether an integral type T can be cast losslessly
   to an integral type U.

   In particular, both T and U must be integral types, must have the same
   signedness, and the maximal permissible value of type T must be no greater
   than that of U. (We assume value ranges of signed types to be essentially
   symmetric around 0).

   Example use:

   template <typename T, integral_fits_t<T, long> = 0 >
   void print(T v) {printf("%ld\n", (long) v);}
   template <typename T, integral_fits_t<T, unsigned long> = 0 >
   void print(T v) {printf("%lu\n", (unsigned long) v);}
*/

/* Note: with gcc 9.2.1, a debug build can't instantiate the full check
 * "both integral + same sign + compatible maxval" lazily, and therefore
 * we get a warning. So we have to resort to an ugly workaround.
 */

template <typename T, typename U>
struct integral_fits_pre_ {
    static constexpr bool value = std::is_integral<T>::value && std::is_integral<U>::value &&
                                  std::is_signed<T>::value == std::is_signed<U>::value;
};

template<bool pre_flag, typename T, typename U>
struct integral_fits_post;

template<typename T, typename U>
struct integral_fits_post<true, T, U> {
    static constexpr bool value = std::numeric_limits<T>::max() <= std::numeric_limits<U>::max();
};
template<typename T, typename U>
struct integral_fits_post<false, T, U> {
    static constexpr bool value = false;
};

template <typename T, typename U>
struct integral_fits_ {
    static constexpr bool value_pre = integral_fits_pre_<T, U>::value;
    static constexpr bool value = integral_fits_post<value_pre, T, U>::value;
};


template<bool> struct integral_fits_final : std::false_type {};
template<> struct integral_fits_final<true> : std::true_type { typedef bool type; };

template <typename T, typename U>
struct integral_fits : integral_fits_final<integral_fits_<T, U>::value> {};

template <typename T, typename U >
using integral_fits_t = typename integral_fits<T, U>::type;

/* Use this for unique_ptr's of objects allocated with malloc() */
template<typename T>
struct free_delete
{
    void operator()(T* x) {
        free(x);        // NOLINT(cppcoreguidelines-no-malloc,hicpp-no-malloc)
    }
};

/* This makes it possible to replace FILE* by std::unique_ptr<FILE>
 * almost transparently. (note that we can't do that with
 * fclose_maybe_compressed, unfortunately, since it requires to keep
 * track of the file name)
 */
template<>
struct std::default_delete<FILE>
{
    void operator()(FILE* x) { fclose(x); }
};

#define CADO_DEFAULT_CXX_CTOR(T)        \
    T() = default

#define CADO_DEFAULT_COPY_CTOR(T)       \
    T(T const &) = default

#define CADO_DEFAULT_COPY_ASSIGNMENT(T) \
    T& operator=(T const &) = default

#define CADO_DEFAULT_MOVE_CTOR(T)       \
    T(T&&) = default

#define CADO_DEFAULT_MOVE_ASSIGNMENT(T) \
    T& operator=(T&&) = default

#define CADO_DEFAULT_ALL_FIVE(T)        \
    CADO_DEFAULT_CXX_CTOR(T);           \
    CADO_DEFAULT_COPY_CTOR(T);          \
    CADO_DEFAULT_COPY_ASSIGNMENT(T);    \
    CADO_DEFAULT_MOVE_CTOR(T);          \
    CADO_DEFAULT_MOVE_ASSIGNMENT(T)


static inline std::vector<std::string> split(
        const std::string& s,
        const std::string& delimiter)
{
    std::vector<std::string> tokens;
    size_t pos = 0, next;
    for ( ;
            (next = s.find(delimiter, pos)) != std::string::npos ;
            pos = next + delimiter.size()
            )
        tokens.push_back(s.substr(pos, next - pos));
    tokens.push_back(s.substr(pos));

    return tokens;
}

/* do we want these on iterables instead of vectors ? */
template<typename ItemType>
static inline std::string join(std::vector<ItemType> const & items,
        const std::string& delimiter)
{
    std::string ret;
    for(auto const & s : items) {
        if (!ret.empty()) ret += delimiter;
        ret += fmt::format("{}", s);
    }
    return ret;
}

template<>
inline std::string join<std::string>(std::vector<std::string> const & items,
        const std::string& delimiter)
{
    std::string ret;
    for(auto const & s : items) {
        if (!ret.empty()) ret += delimiter;
        ret += s;
    }
    return ret;
}

template<typename Lambda, typename Iterable>
static inline std::string join(Iterable first, Iterable last,
        const std::string& delimiter,
        Lambda lambda)
{
    std::string ret;
    for(Iterable x = first ; x != last ; ++x) {
        if (!ret.empty()) ret += delimiter;
        ret += lambda(*x);
    }
    return ret;
}

/* not sure it's really needed. after all, we might want to have a
 * generic map construction
 */
template<typename Lambda, typename ItemType>
static inline std::string join(std::vector<ItemType> const & items,
        const std::string& delimiter,
        Lambda lambda)
{
    return join(items.begin(), items.end(), delimiter, lambda);
}


#endif	/* UTILS_CXX_HPP_ */
