#ifndef CADO_UTILS_CXX_HPP
#define CADO_UTILS_CXX_HPP

#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <array>
#include <cctype>
#include <ios>
#include <istream>
#include <iterator>
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

struct convert_bool {
    template<typename T>
    bool operator()(T const & x) const { return bool(x); }
};

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
 * track of the file name).
 * However it is probably preferrable to define the deleter explicitly
 * with the std::unique_ptr<FILE, delete_FILE> form, which does not rely
 * on presence/absence of the utils_cxx.hpp header which tampers with the
 * std namespace.
 */
template<>
struct std::default_delete<FILE>
{
    void operator()(FILE* x) { fclose(x); }
};

using delete_FILE = std::default_delete<FILE>;

// these macros are not very useful. The added benefit to having all the
// stuff expanded is minor, or even nonexistent.
// NOLINTBEGIN(bugprone-macro-parentheses)
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

#define CADO_DELETE_COPY_CTOR(T)       \
    T(T const &) = delete

#define CADO_DELETE_COPY_ASSIGNMENT(T) \
    T& operator=(T const &) = delete

#define CADO_DELETE_MOVE_CTOR(T)       \
    T(T&&) = delete

#define CADO_DELETE_MOVE_ASSIGNMENT(T) \
    T& operator=(T&&) = delete

#define CADO_DELETE_ALL_FOUR(T)         \
    CADO_DELETE_COPY_CTOR(T);           \
    CADO_DELETE_COPY_ASSIGNMENT(T);     \
    CADO_DELETE_MOVE_CTOR(T);           \
    CADO_DELETE_MOVE_ASSIGNMENT(T)

// NOLINTEND(bugprone-macro-parentheses)


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

template<typename Iterable>
static inline std::string join(Iterable first, Iterable last,
        const std::string& delimiter)
{
    std::string ret;
    for(Iterable x = first ; x != last ; ++x) {
        if (x != first) ret += delimiter;
        ret += fmt::format("{}", *x);
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
        if (x != first) ret += delimiter;
        ret += lambda(*x);
    }
    return ret;
}

template<typename ItemType>
static inline std::string join(std::vector<ItemType> const & items,
        const std::string& delimiter)
{
    return join(items.begin(), items.end(), delimiter);
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

template<typename ItemType, size_t N>
static inline std::string join(std::array<ItemType, N> const & items,
        const std::string& delimiter)
{
    return join(items.begin(), items.end(), delimiter);
}

/* not sure it's really needed. after all, we might want to have a
 * generic map construction
 */
template<typename Lambda, typename ItemType, size_t N>
static inline std::string join(std::array<ItemType, N> const & items,
        const std::string& delimiter,
        Lambda lambda)
{
    return join(items.begin(), items.end(), delimiter, lambda);
}

namespace strip_details {
    struct isspace {
        bool operator()(char c) const {
            return std::isspace(c);
        }
    };
}

template<typename T>
inline std::string& lstrip(std::string &s, T f)
{
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [f](int ch) {
        return !f(ch);
    }));
    return s;
}

template<typename T>
inline std::string & rstrip(std::string &s, T f)
{
    s.erase(std::find_if(s.rbegin(), s.rend(), [f](int ch) {
        return !f(ch);
    }).base(), s.end());
    return s;
}

template<typename T>
inline std::string & strip(std::string & s, T f)
{
    return lstrip(rstrip(s, f), f);
}

template<typename T>
inline std::string lstrip(std::string const & s, T f)
{
    std::string t = s;
    return lstrip(t, f);
}

template<typename T>
inline std::string rstrip(std::string const & s, T f)
{
    std::string t = s;
    return rstrip(t, f);
}

template<typename T>
inline std::string strip(std::string const & s, T f)
{
    std::string t = s;
    return strip(t, f);
}

inline std::string& lstrip(std::string & s) {
    return lstrip(s, strip_details::isspace());
}
inline std::string& rstrip(std::string & s) {
    return rstrip(s, strip_details::isspace());
}
inline std::string& strip(std::string & s) {
    return strip(s, strip_details::isspace());
}
inline std::string lstrip(std::string const & s) {
    return lstrip(s, strip_details::isspace());
}
inline std::string rstrip(std::string const & s) {
    return rstrip(s, strip_details::isspace());
}
inline std::string strip(std::string const & s) {
    return strip(s, strip_details::isspace());
}


/* use this as: input_stream >> read_container(container, maximum_size)
 * or possibly: input_stream >> read_container(container)
 */
template<typename Container>
struct read_container_impl {
    Container & c;
    typename Container::size_type N;
};
template<typename Container>
read_container_impl<Container>
read_container(Container & c,
    typename Container::size_type N = std::numeric_limits<typename Container::size_type>::max())
{
    return read_container_impl<Container> { c, N };
}

template<typename Container>
std::istream& operator>>(std::istream& is, read_container_impl<Container>&& R)
{
    R.c.clear();
    for(typename Container::size_type i = 0 ; is && i < R.N ; ++i) {
        typename Container::value_type v;
        if (is >> v)
            R.c.emplace_back(std::move(v));
    }
    return is;
}

template<typename T>
void checked_realloc(T * & var, size_t N)
{
    if (!(N)) {
        free(var);
        (var) = nullptr;
    } else {
        auto * p = static_cast<T *>(realloc((var), (N) * sizeof(T)));
        if (!p && (var) != nullptr)
            free((var));
        ASSERT_ALWAYS(p != nullptr);
        (var) = p;
    }
}

static inline std::unique_ptr<FILE, delete_FILE> fopen_helper(std::string const & filename, const char * mode, bool accept_failure = false)
{
    std::unique_ptr<FILE, delete_FILE> f(fopen(filename.c_str(), mode));
    DIE_ERRNO_DIAG(!f && !accept_failure, "fopen(%s)", filename.c_str());
    return f;
}

struct decomposed_path : public std::vector<std::string> {
    decomposed_path() : decomposed_path("/") {}
    explicit decomposed_path(std::string const &);
    explicit decomposed_path(const char * s) : decomposed_path(std::string(s)) {}
    explicit operator std::string() const;
    std::string dirname() const;
    std::string basename() const { return back(); }
    bool is_relative() const;
    bool is_absolute() const { return !is_relative(); }
    std::string extension() const;
};

namespace cado::details {
    template<typename T, typename ... Args>
    struct double_ratio_impl {
        double operator()(T const & t, Args&&...);
    };
    template<typename T>
    struct double_ratio_impl<T> {
        double operator()(T const & t) { return t; }
    };
    template<typename T, typename U, typename ... Args>
    struct double_ratio_impl<T, U, Args...> {
        double operator()(T const & t, U const & u, Args&&... args) {
            return !u ? 0 : double_ratio_impl<T, Args...>()(static_cast<double>(t) / static_cast<double>(u), std::forward<Args>(args)...);
        }
    };
} /* namespace cado::details */

template<typename T, typename ... Args>
double double_ratio(T const & t, Args&&... args)
{
    return cado::details::double_ratio_impl<T, Args...>()(t, std::forward<Args>(args)...);
}

template<int N>
struct expect_s
{
    const char * s;
    explicit expect_s(const char s0[N]) : s(s0) {}
};

template<int N>
static expect_s<N> expect(char const (&s0)[N]) { return expect_s<N> { s0 }; }

template<int N>
static std::istream& operator>>(std::istream& is, expect_s<N> const & e)
{
    char t[N];
    is.get(t, N);  // side-
    if (strcmp(t, e.s) != 0)
        is.setstate(std::ios::failbit);
    return is;
}

struct expect_string {
    std::string s;
};

static inline expect_string expect(std::string const & s) { return { s }; }


static inline std::istream& operator>>(std::istream & is, expect_string const & e)
{
    if(!std::equal(std::begin(e.s), std::end(e.s), std::istreambuf_iterator<char>{is})) {
        is.setstate(is.rdstate() | std::ios::failbit);
    }
    return is;
}

namespace cado {
    /* example:
     *
        foo_type compute_foo() const;
        cado::cached_property<foo_type> cached_foo;
        foo_type const & foo() const {
            return cached_foo([this](){compute_foo();});
        }
     */

    template<typename T> struct cached_property {
        mutable std::unique_ptr<T> c;
        template<typename F>
        T const & operator()(F f) const
        {
            /* TODO: locking.
             */
            if (!c)
                c = std::make_unique<T>(f());
            return *c;
        }
        bool is_cached() const { return c.get(); }
        cached_property() = default;
        ~cached_property() = default;
        /* make these copyable - movable. The cache is not retained when
         * we copy, but it is kept when we move.
         */
        cached_property(cached_property const&) {}
        cached_property& operator=(cached_property const& o) {
            if (this != &o)
                c.reset();
            return *this;
        }
        cached_property(cached_property &&) noexcept = default;
        cached_property& operator=(cached_property && o) noexcept = default;
    };
} /* namespace cado */

namespace cado {
    /* A type trait that checks whether an integral type T can be cast
     * losslessly to an integral type U.
     *
     * In particular, both T and U must be integral types, must have the
     * same signedness, and the maximal permissible value of type T must
     * be no greater than that of U. (We assume value ranges of signed
     * types to be essentially symmetric around 0).
     */

    /* Note: the full check
     * "both integral + same sign + compatible maxval" cannot be
     * evaluated lazily, so we need an ugly workaround.
     */

    template <typename T, typename U, bool =
        std::is_integral_v<T> && // E: Template argument for template t…
        std::is_integral_v<U> && 
        std::is_signed_v<T> == std::is_signed_v<U>>
        struct integral_fits_aux {
            static constexpr bool value = std::numeric_limits<T>::max() <= std::numeric_limits<U>::max();
        };
    template <typename T, typename U>
        struct integral_fits_aux<T, U, false> : public std::false_type {};
    template <typename T, typename U>
        inline constexpr bool integral_fits_v = integral_fits_aux<T, U>::value;

    /* This is an alternative implementation which uses a different
     * mechanism.
     */
    namespace is_narrowing_conversion_detail
    {
        template<typename From, typename To, typename = void>
            struct is_narrowing_conversion_impl : std::true_type {};
        template<typename From, typename To, typename = void>
            struct is_non_narrowing_conversion_impl : std::false_type {};

        template<typename From, typename To>
            struct is_narrowing_conversion_impl<From, To, std::void_t<decltype(To{std::declval<From>()})>> : std::false_type {};
        template<typename From, typename To>
            struct is_non_narrowing_conversion_impl<From, To, std::void_t<decltype(To{std::declval<From>()})>> : std::true_type {};
    }  /* namespace is_narrowing_conversion_detail */

    template<typename From, typename To>
        struct is_narrowing_conversion : is_narrowing_conversion_detail::is_narrowing_conversion_impl<From, To> {};
    template<typename From, typename To>
        struct is_non_narrowing_conversion : is_narrowing_conversion_detail::is_non_narrowing_conversion_impl<From, To> {};

    template<typename From, typename To>
    inline constexpr bool is_narrowing_conversion_v = is_narrowing_conversion<From, To>::value;
    template<typename From, typename To>
    inline constexpr bool is_non_narrowing_conversion_v = is_non_narrowing_conversion<From, To>::value;

    /* this trait can be used in ctors and assignment operators, for
     * instance */
    template<typename T, typename U>
    static constexpr bool converts_via =
        integral_fits_v<T, U> &&
        is_non_narrowing_conversion_v<T, U>;

} /* namespace cado */


#endif	/* CADO_UTILS_CXX_HPP */
