#ifndef CADO_PARAMS_HPP
#define CADO_PARAMS_HPP

#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <charconv>
#include <map>
#include <ostream>
#include <istream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include <array>
#include <mutex>
#include <concepts>
#include <ranges>
#if __cpp_lib_to_chars >= 201611
#else
/* libstdc++-10 does not have std::from_chars for floating-point types
 * (those were introduced in 11.1:
 * https://gcc.gnu.org/onlinedocs/libstdc++/manual/status.html#:~:text=Elementary%20string%20conversions)
 * we fall back to stringstream parsing instead
 */
#include <sstream>
#endif

#include "fmt/base.h"
#include "fmt/format.h"

#include <gmp.h>

#include "macros.h"
#include "utils_cxx.hpp"

namespace cado::params
{

struct parameter_error : public std::runtime_error {
    template<typename... Args>
    explicit parameter_error(fmt::format_string<Args...> s, Args && ...args)
        : std::runtime_error(fmt::format(s, std::forward<Args>(args)...))
    {
    }
    explicit parameter_error(std::string const & s)
        : std::runtime_error(s)
    {
    }
};

enum class origin { FROM_FILE, FROM_CMDLINE };

struct copy_previous_side {
};

template<typename T, typename... Args>
bool parse(std::string const & s, T & value, Args && ...args);

template<typename T> struct parser;

template<std::floating_point T> struct parser<T> {// {{{
    bool operator()(std::string const & s, T & value) const
    {
#if __cpp_lib_to_chars >= 201611
        const char * b = s.data();
        const char * e = s.data() + s.size();
        auto [ptr, ec] = std::from_chars(b, e, value);
        return ptr == e;
#else
        std::istringstream ss(s);
        return ss >> value && ss.eof();
#endif
    }
};
// }}}
template<std::integral T> struct parser<T> {// {{{
    /* integral types get loose parsing: floating point expressions such
     * as 2e6 are parsed as integers as well.
     * More precisely we recognize any integer in one of the following forms:
     *      <some integer>
     *      <somd double that casts properly to T **AND BACK**>
     *
     * note that we want the full string to match!
     */

    bool operator()(std::string const & s, T & value) const
    {
#if __cpp_lib_to_chars >= 201611
        const char * b = s.data();
        const char * e = s.data() + s.size();
        if (auto [ptr, ec] = std::from_chars(b, e, value); ptr == e)
            return true;
#else
        if (std::istringstream ss(s); ss >> value && ss.eof())
            return true;
#endif

        double d;
#if __cpp_lib_to_chars >= 201611
        if (auto [ptr, ec] = std::from_chars(b, e, d); ptr != e)
            return false;
#else
        if (std::istringstream ss(s); !(ss >> d && ss.eof()))
            return false;
#endif
        value = T(d);
        if (double(value) == d)
            return true;
        else
            throw parameter_error("incorrect double-to-integer conversion in parameter");
    }
};
// }}}
template<> struct parser<std::string> {// {{{
    bool operator()(std::string const & s, std::string & value) const
    {
        value = s;
        return true;
    }
};
// }}}
template<typename T> struct parser<std::vector<T>> {// {{{
    bool operator()(std::string const & s, std::vector<T> & value, std::string const & sep = ",") const
    {
        std::vector<T> w;

        for(auto const & t : split(s, sep)) {
            T v {};
            if (!parse(t, v))
                return false;
            w.push_back(v);
        }
        value = std::move(w);
        return true;
    }
};
// }}}
template<typename T, size_t N> struct parser<T[N]> {// {{{
    bool operator()(std::string const & s, T (&value)[N], std::string const & sep = ",") const
    {
        std::vector<T> w;
        parse(s, w, sep);
        if (w.size() != N)
            return false;
        std::copy(w.begin(), w.end(), value);
        return true;
    }
    bool operator()(std::string const & s, T (&value)[N], std::string const & sep, T const & dfl) const
    {
        std::vector<T> w;
        parse(s, w, sep);
        if (w.size() > N)
            return false;
        for( ; w.size() < N ; w.push_back(dfl));
        std::copy(w.begin(), w.end(), value);
        return true;
    }
};
// }}}

/* This one matters because we don't want the generic array parser to
 * catch mpz_t ! The implementation is in gmp_aux2.cpp
 */
template<> struct parser<mpz_t> {// {{{
    bool operator()(std::string const & s, mpz_t value) const;
};
// }}}
template<typename T, size_t N> struct parser<std::array<T, N>> {// {{{
    bool operator()(std::string const & s, std::array<T, N> & value, std::string const & sep = ",") const
    {
        std::vector<T> w;
        parse(s, w, sep);
        if (w.size() != N)
            return false;
        std::ranges::copy(w, value.begin());
        return true;
    }
    bool operator()(std::string const & s, std::array<T, N> & value, std::string const & sep, T const & dfl) const
    {
        std::vector<T> w;
        parse(s, w, sep);
        if (w.size() > N)
            return false;
        for( ; w.size() < N ; w.push_back(dfl));
        std::ranges::copy(w, value.begin());
        return true;
    }
};
// }}}
template<typename T, typename U> struct parser<std::pair<T, U>> {// {{{
    bool operator()(std::string const & s, std::pair<T, U> & value, std::string const & sep = ",") const
    {
        auto toks = split(s, sep);
        if (toks.size() != 2)
            return false;
        std::pair<T, U> w;
        if (!parse(toks[0], w.first) || !parse(toks[1], w.second))
            return false;
        value = std::move(w);
        return true;
    }
};
// }}}
template<> struct parser<bool> {// {{{
    bool operator()(std::string const & s, bool & b) const {
        if (s == "1" || s == "true" || s == "True" || s == "yes" || s == "on") {
            b = true;
            return true;
        } else if (s == "0" || s == "false" || s == "False" || s == "no" || s == "off") {
            b = false;
            return true;
        }
        return false;
    }
};// }}}


/* to declare a parser for a custom type, provide a (full) specialization
 * of cado::params::from_chars. For example:
 *
template<>
struct parser<cxx_mpz> {
    bool operator()(std::string const &, T & value) {
    ...
    }
};
value);
 */

template<typename T, typename... Args>
bool parse(std::string const & s, T & value, Args && ...args)
{
    return parser<T>()(s, value, std::forward<Args>(args)...);
}

struct cxx_param_list {
    // documented parameters
    std::string usage_header;
    struct collate {
        /* compare two strings, intentionally collating - and _ (except at the
         * beginning of the string) */
        bool operator()(std::string const & a, std::string const & b) const
        {
            size_t k;
            for (k = 0; k < a.size() && k < b.size(); k++) {
                int r = (a[k] > b[k]) - (b[k] > a[k]);
                if (k && (a[k] == '-' || a[k] == '_') &&
                    (b[k] == '-' || b[k] == '_'))
                    r = 0;
                if (r)
                    return r < 0;
            }
            ASSERT_ALWAYS(k == a.size() || k == b.size());
            return a.size() < b.size();
            // return (k < b.size()) < (k < a.size());
            // return (a[k] > b[k]) < (b[k] > a[k]);
        }
    };
    std::map<std::string, std::string, collate> documentation;
    std::vector<std::string> documentation_structure;

    struct parameter {
        std::string value;
        origin from;
        bool parsed;
        int seen;
        explicit parameter(std::string value = {},
                           origin from = origin::FROM_FILE, bool parsed = false,
                           int seen = 1)
            : value(std::move(value))
            , from(from)
            , parsed(parsed)
            , seen(seen)
        {
        }
        parameter(parameter const &) = default;
        parameter & operator=(parameter const &) = default;
        parameter(parameter &&) = default;
        parameter & operator=(parameter &&) = default;
        ~parameter() = default;
    };
    std::map<std::string, parameter, collate> p;
    std::map<std::string, std::string, collate> aliases;
    std::map<std::string, int *, collate> switches;

    bool empty() const { return p.empty(); }

    /* We use this to remember the first command line pointer which have
     * been given to us */
    int cmdline_argc0 = 0;
    char const ** cmdline_argv0 = nullptr;
    // did the user use the doc functionality ?
    bool use_doc = false;

    private:
    cxx_param_list(cxx_param_list const& o,
            std::scoped_lock<std::mutex> const &)
        : usage_header(o.usage_header)
        , documentation(o.documentation)
        , documentation_structure(o.documentation_structure)
        , p(o.p)
        , aliases(o.aliases)
        , switches(o.switches)
        , cmdline_argc0(o.cmdline_argc0)
        , cmdline_argv0(o.cmdline_argv0)
        , use_doc(o.use_doc)
    {}
    void swap(cxx_param_list & o,
            std::scoped_lock<std::mutex> const &,
            std::scoped_lock<std::mutex> const &) noexcept
    {
        std::swap(usage_header, o.usage_header);
        std::swap(documentation, o.documentation);
        std::swap(documentation_structure, o.documentation_structure);
        std::swap(p, o.p);
        std::swap(aliases, o.aliases);
        std::swap(switches, o.switches);
        std::swap(cmdline_argc0, o.cmdline_argc0);
        std::swap(cmdline_argv0, o.cmdline_argv0);
        std::swap(use_doc, o.use_doc);
    }
    cxx_param_list(cxx_param_list & o,
            std::scoped_lock<std::mutex> const & ao) noexcept
    {
        swap(o, std::scoped_lock(mutex), ao);
    }

    public:
    /* I believe it's only used once or twice in las_todo_list.cpp */
    cxx_param_list(cxx_param_list const& o)
        : cxx_param_list(o, std::scoped_lock(o.mutex))
    { }

    /* honestly, we only have these to stop clang-tidy from whining */
    cxx_param_list& operator=(cxx_param_list const& o)
    {
        if (this != &o) {
            cxx_param_list c(o, std::scoped_lock(o.mutex));
            swap(c, std::scoped_lock(mutex), std::scoped_lock(c.mutex));
        }
        return *this;
    }
    cxx_param_list(cxx_param_list && o) noexcept
        : cxx_param_list(o, std::scoped_lock(o.mutex))
    { }
    cxx_param_list& operator=(cxx_param_list && o) noexcept
    {
        if (this != &o) {
            cxx_param_list c(o, std::scoped_lock(o.mutex));
            swap(c, std::scoped_lock(mutex), std::scoped_lock(c.mutex));
        }
        return *this;
    }

    ~cxx_param_list() = default;
    cxx_param_list() = default;

    /* Whether the most convenient interface is lookup() or has() or
     * parse(), I'm not really sure
     */
    std::string lookup(std::string const & key)
    {
        if (auto const * c = get_assoc_ptr(key))
            return *c;
        else
            return {};
    }
    std::string const * has(std::string const & key) {
        return get_assoc_ptr(key);
    }
    const char * lookup_old(std::string const & key) {
        auto const * t = get_assoc_ptr(key);
        if (!t) return nullptr;
        return t->c_str();
    }

    void declare_usage_header(std::string const & header) {
        usage_header = header;
    }
    void declare_usage_section(std::string const & header) {
        documentation_structure.push_back(std::string(" ") + header);
    }
    void declare_usage(std::string const & key, std::string const & doc) {
        documentation_structure.push_back(key);
        documentation[key] = doc;
        /* Note that duplicate calls to param_list_decl_usage for the same
         * key will not trigger two distinct prints of the same
         * documentation string. See the collapsing logic in
         * param_list_print_usage.
         *
         * XXX hmmm. does this observation still hold?
         */
        use_doc = true;
    }
    bool is_documented(std::string const & key) const {
        return documentation.find(key) != documentation.end();
    }

    void print_usage(FILE * = stderr) const;
    void print_usage(std::ostream &) const;

    // for debugging.
    void display_debug(FILE * f) const
    {
        for(auto const & [key, value] : p)
            fmt::print(f, "{}={}\n", key, value.value);
    }


    // This function is a shorthand which does employ some hackery put into
    // param lists, which remember their oldest argv, argc pair.
    void print_command_line(FILE * stream) const;

    // A switch is a command-line argument which sets a value by its mere
    // presence. Could be for instance --verbose, or --use-smart-algorithm
    void configure_switch(std::string const & key);
    void configure_switch_old(std::string const & switchname, int * p);

    // This one allows shorthands. Notice that the alias string has to
    // contain the exact form of the wanted alias, which may be either "-x",
    // "--x", or "x=" (a terminating = tells the program that the option is
    // wanted all in one go, like in ./a.out m=42, in contrast to ./a.out -m
    // 42).
    void configure_alias(std::string const & key, std::string const & alias);

    template <typename... Args>
    void fail(fmt::format_string<Args...> s, Args &&... args) const
    {
        std::string ss = fmt::format(s, std::forward<Args>(args)...);
        if (!ss.ends_with("\n"))
            ss += "\n";
        // fmt::print(stderr, "{}", ss);
        print_usage(stderr);
        throw std::runtime_error(ss);
    }

    std::string binary_name() const
    {
        if (auto * t = cmdline_argv0)
            return t[0];
        else
            return {};
    }

    /* This does the following:
     *  - save argv[0]
     *  - argv++, argc--
     *  - update the pl structure with all the arguments that are found.
     *  - process --help to display the help message as configured via
     *  param_list_decl_usage
     *  - return with some trailing arguments if they are allowed, otherwise
     *  error out.
     */
    void process_command_line(int & argc, char const **& argv,
                              bool accept_trailing_args = false);
    void process_command_line_and_extra_parameter_files(int & argc,
                                                        char const **& argv);
    // sees whether the arguments pointed to by argv[0] and (possibly)
    // argv[1] correspond to either -<key> <value>, --<key> <value> or
    // <key>=<value> ; configured switches and aliases for the param list are
    // also checked.
    int update_cmdline(int & argc, char const **& argv);

    // warns against unused command-line parameters. This normally indicates
    // a user error. parameters ignored from config files are considered
    // normal (although note that in some cases, it could be bad as well).
    // the number of unused parameters is returned.
    int warn_unused() const;

    // takes a file in the Cado-NFS params format and stores the dictionary
    // of parameters.
    int read(std::istream & is, bool stop_on_empty_line=false);
    int read(FILE * f, bool stop_on_empty_line=false);
    int read(std::string const & filename)
    {
        std::ifstream is(filename);
        if (!is)
            fail("Cannot open {}", filename);
        return read(is);
    }


    // this one is the ``joker'' call.
    void add_key(std::string const & key, std::string const & value,
                 origin from);

    // removing a key can be handy before saving a config file. Some options
    // are relevant only for one particular invokation, and not for saving.
    void remove_key(std::string const & key) {
        auto it = p.find(key);
        if (it != p.end())
            p.erase(it);
    }

    private:
    template<typename T, typename... Args>
    bool parse_base(bool stealth, std::string const & key, T & r, Args&& ...args)
    {
        std::string value;
        std::string diagnostic;
        if (!get_assoc(key, value, stealth))
            return false;
        try {
            if (cado::params::parse(value, r, std::forward<Args>(args)...))
                return true;
        } catch (parameter_error const & e) {
            diagnostic = e.what();
        }

        auto f = fmt::format("cannot cast parameter {} to type {}",
                key, typeid(T).name());
        if (!diagnostic.empty())
            f += fmt::format(": {}", diagnostic);
        throw parameter_error(f);
    }
    public:
    template<typename T, typename... Args>
    bool parse(std::string const & key, T & r, Args&& ...args)
    {
        return parse_base(false, key, r, std::forward<Args>(args)...);
    }
    template<typename T, typename... Args>
    bool parse_mandatory(std::string const & key, T & r, Args&& ...args)
    {
        if (!parse(key, r, std::forward<Args>(args)...))
            fail("missing mandatory argument {}", key);
        return true;
    }
    template<typename T, typename... Args>
    bool parse_stealth(std::string const & key, T & r, Args&& ...args)
    {
        return parse_base(true, key, r, std::forward<Args>(args)...);
    }

    /* this returns a value-initialized T if the key is absent */
    template <typename T> T parse(std::string const & key)
    {
        T r {};
        parse<T>(key, r);
        return r;
    }

    template <typename T> T parse_optional(std::string const & key, T const & dfl)
    {
        T r = dfl;
        parse<T>(key, r);
        return r;
    }

    template <typename T> T parse_mandatory(std::string const & key)
    {
        T r {};
        parse_mandatory<T>(key, r);
        return r;
    }

    size_t get_list_count(std::string const & key, std::string const & sep = ",");

    private:
    /* This allows a mixed syntax for the specification of a list of values.
     * We accept either:
     *  - a comma-separated list of values, e.g., lpb=27,29
     *    (for historical reasons, lpbs=27,29 also works).
     *  - OR values being specified one by one, e.g., lpb0=27 lpb1=29
     *
     * This function definitely wants all resulting values to be set, so that
     * a policy has to be defined for the missing values (if the items are
     * specified in the "direct" way as with lpb1).
     *
     *  - ARGS_PER_SIDE_DEFAULT_AS_IS: we leave the missing values exactly as
     *    they were on input
     *  - ARGS_PER_SIDE_DEFAULT_COPY_PREVIOUS: if lpb<n> missing and lpb<n-1>
     *    is set, use the latter as a default value.
     *
     * RETURN VALUE: badly specified as of now, please refrain from
     * using.
     */
    template<typename T>
    std::vector<T> parse_per_side_base(std::string const & key0, size_t n)
    {
        auto key = drop_one_or_two_leading_dashes(key0);
        bool has_lpb01 = false;
        std::vector<T> lpb_arg;
        for(size_t side = 0 ; side < n ; side++) {
            auto keyi = fmt::format("{}{}", key, side);
            T r {};
            has_lpb01 = parse_stealth(keyi, r);
            lpb_arg.push_back(r);
        }

        bool has_nlpbs;
        {
            if (!(has_nlpbs = parse(key, lpb_arg, ","))) {
                char const c = std::string(key).back();
                if (c != 'x' && c != 's') {
                    /* try with s. */
                    auto keys = key + "s";
                    has_nlpbs = parse_stealth(keys, lpb_arg, ",");
                }
            }
        }

        if (has_nlpbs && has_lpb01)
            throw parameter_error("{0}s[01] and {0}s are incompatible\n", key);

        if (!has_nlpbs && !has_lpb01)
            return {};

        return lpb_arg;
    }

    public:

    /* This returns a vector that is either empty or has size exactly n,
     * based on the value of the parameters ${key0}${i} for i = 0 ...
     * n-1, or the comma separated list given by parameter ${key0}s.
     * Missing specifications are copied from the value with the
     * preceding index.
     */
    template <typename T>
    bool parse_per_side(std::string const & key0, std::vector<T> & r, size_t n, copy_previous_side const &)
    {
        r = parse_per_side_base<T>(key0, n);
        if (r.empty())
            return false;

        if (r.size() > n)
            throw parameter_error(
                    "Number of values for parameter {} exceeds the maximum {}.",
                    key0, n);
        for(size_t i = r.size() ; i < n ; i++)
            r.push_back(r.back());
        return true;
    }

    /* Same as above, but missing values are filled with the passed
     * default.
     */
    template <typename T>
    bool parse_per_side(std::string const & key0, std::vector<T> & r, size_t n, T const & v)
    {
        r = parse_per_side_base<T>(key0, n);
        if (r.empty())
            return false;

        if (r.size() > n)
            throw parameter_error(
                    "Number of values for parameter {} exceeds the maximum {}.",
                    key0, n);

        for(size_t i = r.size() ; i < n ; i++)
            r.push_back(v);
        return true;
    }

    /* accepts no defaults */
    template <typename T>
    bool parse_per_side(std::string const & key0, std::vector<T> & v, size_t n)
    {
        v = parse_per_side_base<T>(key0, n);
        if (v.size() != n)
            v.clear();
        return v.size() == n;
    }
    template <typename T, size_t N, typename... Args>
    bool parse_per_side(std::string const & key0, std::array<T, N> & v, Args && ...args)
    {
        std::vector<T> tv;
        if (!parse_per_side(key0, tv, N, std::forward<Args>(args)...))
            return false;
        std::copy(tv.begin(), tv.end(), v.begin());
        return true;
    }
    /* do we want to accept defaults here? */
    template <typename T, size_t N, typename ...Args>
    bool parse_per_side(std::string const & key0, T (&v)[N], Args &&...args)
    {
        std::vector<T> tv;
        if (!parse_per_side(key0, tv, N, std::forward<Args>(args)...))
            return false;
        std::copy(tv.begin(), tv.end(), v);
        return true;
    }

    template <typename... Args>
    bool parse_per_side_mandatory(std::string const & key0, Args &&...args)
    {
        if (!parse_per_side(key0, std::forward<Args>(args)...))
            fail("missing mandatory argument {0}0, {0}1, ... (or {0}s)", key0);
        return true;
    }
    /* this one is convenient in ctors */
    template <typename T, typename... Args>
    T parse_per_side_mandatory(std::string const & key0, Args &&...args)
    {
        T res;
        parse_per_side_mandatory(key0, res, std::forward<Args>(args)...);
        return res;
    }


    private:
    std::string const * get_assoc_ptr(std::string const & key0, bool stealth = false, bool * seen = nullptr);
    bool get_assoc(std::string const & key, std::string & value, bool stealth = false, bool * seen = nullptr);
    mutable std::mutex mutex;
    static std::string drop_one_or_two_leading_dashes(std::string s)
    {
        size_t i = 0;
        for( ; i < 2 && i < s.size() && s[i] == '-' ; i++) ;
        if (i) s = s.substr(i);
        return s;
    }

};

/* this one is really just a proxy around join() in utils_cxx.hpp */
extern std::string collect_command_line(int argc, char const * argv[]);

struct flags {
    using type = uint8_t;
    static constexpr type optional = 0;
    static constexpr type mandatory = 1;
    static constexpr type toggle = 2;
};

template <flags::type fl, typename T, cado::string_literal key,
          cado::string_literal documentation, cado::string_literal default_value>
struct parameter_base {
    static constexpr bool has_default = !decltype(default_value)::empty();
    using is_parameter_base = std::true_type;
    static void declare_usage(cxx_param_list & pl)
    {
        if constexpr (has_default) {
            auto f =
                fmt::format("{} (default {})", documentation, default_value);
            pl.declare_usage(key.value, f);
        } else {
            pl.declare_usage(key.value, documentation.value);
        }
    }
    static void configure(cxx_param_list & pl)
    {
        declare_usage(pl);
        if (fl & flags::toggle)
            pl.configure_switch(key.value);
    }
    static T get_default_value() {
        T x {};
        std::string s(default_value.begin(), default_value.end());
        if (s.empty() || !cado::params::parse(s, x))
            return {};
        return x;
    }
    T value = get_default_value();
    explicit parameter_base(cxx_param_list & pl)
    {
        bool parsed = pl.parse(key.value, value);
        if (!parsed && (fl & flags::mandatory))
            pl.fail("missing mandatory parameter {}", key);
    }
    parameter_base() = default;
    /* Implicit conversions are really, really desired */
    // NOLINTBEGIN(google-explicit-constructor,hicpp-explicit-conversions)
    operator T const &() const { return value; }
    operator T &() { return value; }
    // NOLINTEND(google-explicit-constructor,hicpp-explicit-conversions)

    /* when implicit conversions don't work, we have parameter_value() */
    T const & parameter_value() const { return value; }
    T & parameter_value() { return value; }
    T const & operator()() const { return value; }
    T & operator()() { return value; }
    bool is_default() const { return value == get_default_value(); }
    bool is_provided() const { return !is_default(); }
    using parameter_type = T;
};
} /* namespace cado::params */

/* The "with_default" options should be templates with default args,
 * obviously. Alas, this does not seem to work with g++-10 (debian-11)
 */
template <typename T, cado::string_literal key,
          cado::string_literal documentation>
using parameter =
    cado::params::parameter_base<cado::params::flags::optional, T, key,
                                 documentation, "">;

template <typename T, cado::string_literal key,
          cado::string_literal documentation,
          cado::string_literal default_value>
using parameter_with_default =
    cado::params::parameter_base<cado::params::flags::optional, T, key,
                                 documentation, default_value>;

/* a mandatory parameter does not have a default value. */
template <typename T, cado::string_literal key,
          cado::string_literal documentation>
using parameter_mandatory =
    cado::params::parameter_base<cado::params::flags::mandatory, T, key,
                                 documentation, "">;

/* in most cases, a switch does not have a default value, but
 * we can have a switch that can only be toggled *down*
 */
template <cado::string_literal key, cado::string_literal documentation>
using parameter_switch =
    cado::params::parameter_base<cado::params::flags::toggle, bool, key,
                                 documentation, "">;

template <cado::string_literal key,
         cado::string_literal documentation,
         cado::string_literal default_value>
using parameter_switch_with_default =
    cado::params::parameter_base<cado::params::flags::toggle, bool, key,
                                 documentation, default_value>;

using cxx_param_list = cado::params::cxx_param_list;

/***********************************************************************/
/* compatibility calls */
static inline void param_list_decl_usage_header(cxx_param_list & pl, const char * a)
{
    pl.declare_usage_header(a);
}
static inline void param_list_decl_usage(cxx_param_list & pl, const char * key, const char * msg)
{
    pl.declare_usage(key, msg);
}
static inline void param_list_print_usage(cxx_param_list & pl, const char *, FILE * f)
{
    pl.print_usage(f);
}
static inline void param_list_print_command_line(FILE * f, cxx_param_list & pl)
{
    pl.print_command_line(f);
}
static inline int param_list_read_stream(cxx_param_list & pl, FILE * f, int stop_on_empty_line)
{
    return pl.read(f, stop_on_empty_line);
}
static inline void param_list_configure_alias(cxx_param_list & pl, const char * key, const char * alias)
{
    pl.configure_alias(key, alias);
}
static inline void param_list_configure_switch(cxx_param_list & pl, const char * key, int * v)
{
    pl.configure_switch_old(key, v);
}
static inline int param_list_parse_double(cxx_param_list & pl, const char * key, double * d)
{
    return pl.parse(key, *d);
}
static inline int param_list_parse_uint(cxx_param_list & pl, const char * key, unsigned int * d)
{
    return pl.parse(key, *d);
}
static inline int param_list_parse_int(cxx_param_list & pl, const char * key, int * d)
{
    return pl.parse(key, *d);
}
static inline int param_list_parse_uint64(cxx_param_list & pl, const char * key, uint64_t * d)
{
    return pl.parse(key, *d);
}
static inline int param_list_parse_int64(cxx_param_list & pl, const char * key, int64_t * d)
{
    return pl.parse(key, *d);
}
static inline int param_list_parse_ulong(cxx_param_list & pl, const char * key, unsigned long * d)
{
    return pl.parse(key, *d);
}
static inline int param_list_parse_long(cxx_param_list & pl, const char * key, long * d)
{
    return pl.parse(key, *d);
}
static inline int param_list_parse_switch(cxx_param_list & pl, const char * key)
{
    /* for the needs of cumulative switches, it is better to return the
     * switch as an int.
     */
    return pl.parse<int>(key);
}
extern int param_list_parse_mpz(cxx_param_list & pl, const char * key, mpz_ptr x);
static inline int param_list_parse_uint_and_uint(cxx_param_list & pl, const char * key, unsigned int * x, const char * sep)
{
    std::array<unsigned int, 2> a {};
    int const r = pl.parse(key, a, sep);
    if (r)
        std::ranges::copy(a, x);
    return r;
}
static inline int param_list_parse_int_and_int(cxx_param_list & pl, const char * key, int * x, const char * sep)
{
    std::array<int, 2> a {};
    int const r = pl.parse(key, a, sep);
    if (r)
        std::ranges::copy(a, x);
    return r;
}
static inline int param_list_parse_intxint(cxx_param_list & pl, const char * key, int * x)
{
    return param_list_parse_int_and_int(pl, key, x, "x");
}
static inline int param_list_parse_int_list(cxx_param_list & pl, const char * key, int * x, size_t n, const char * sep)
{
    /* XXX param_list_parse_*_list gives a _maximum_ size, and expects
     * the number of parsed items as return value */
    std::vector<int> a;
    if (!pl.parse(key, a, sep))
        return 0;
    ASSERT_ALWAYS(a.size() <= n);
    std::ranges::copy(a, x);
    return static_cast<int>(a.size());
}
static inline int param_list_parse_uint_list(cxx_param_list & pl, const char * key, unsigned int * x, size_t n, const char * sep)
{
    /* XXX param_list_parse_*_list gives a _maximum_ size, and expects
     * the number of parsed items as return value */
    std::vector<unsigned int> a;
    if (!pl.parse(key, a, sep))
        return 0;
    ASSERT_ALWAYS(a.size() <= n);
    std::ranges::copy(a, x);
    return static_cast<int>(a.size());
}
template<typename T>
T param_list_parse(cxx_param_list & pl, std::string const & key)
{
    T r {};
    pl.parse(key, r);
    return r;
}
template<typename T>
int param_list_parse(cxx_param_list & pl, std::string const & key, T & r)
{
    return pl.parse(key, r);
}
template<typename T>
T param_list_parse_mandatory(cxx_param_list & pl, std::string const & key)
{
    return pl.parse_mandatory<T>(key);
}
static inline void param_list_generic_failure(cxx_param_list & pl, const char *missing)
{
    pl.fail("missing or invalid parameter {}", missing);
}

static inline const char * param_list_lookup_string(cxx_param_list & pl, const char * key)
{
    auto const * p = pl.has(key);
    return p ? p->c_str() : nullptr;
}

enum args_per_side_policy_t {
    ARGS_PER_SIDE_DEFAULT_AS_IS,
    ARGS_PER_SIDE_DEFAULT_COPY_PREVIOUS,
};
template<typename T>
int param_list_parse_per_side(cxx_param_list & pl, std::string const & key0, T * lpb_arg, int n, enum args_per_side_policy_t policy)
{
    std::vector<T> v;
    int r;
    if (policy == ARGS_PER_SIDE_DEFAULT_COPY_PREVIOUS) {
        r = pl.parse_per_side(key0, v, n, cado::params::copy_previous_side());
    } else {
        r = pl.parse_per_side(key0, v, n, lpb_arg[n-1]);
    }
    if (r)
        std::copy(v.begin(), v.end(), lpb_arg);
    /* XXX param_list_parse_*_per_side expects a size (but it is
     * particularly badly defined, so we should really not use it IMHO) */
    return r * n;
}

static inline int param_list_parse_int_args_per_side(cxx_param_list & pl, std::string const & key0, int * lpb_arg, int n, enum args_per_side_policy_t policy)
{
    return param_list_parse_per_side(pl, key0, lpb_arg, n, policy);
}
static inline int param_list_parse_uint_args_per_side(cxx_param_list & pl, std::string const & key0, unsigned int * lpb_arg, int n, enum args_per_side_policy_t policy)
{
    return param_list_parse_per_side(pl, key0, lpb_arg, n, policy);
}
enum parameter_origin { PARAMETER_FROM_FILE, PARAMETER_FROM_CMDLINE };
static inline void param_list_add_key(cxx_param_list & pl,
        const char * key, const char * value, enum parameter_origin o)
{
    pl.add_key(key, value, o == PARAMETER_FROM_FILE ? cado::params::origin::FROM_FILE : cado::params::origin::FROM_CMDLINE);
}
static inline void param_list_remove_key(cxx_param_list & pl,
        const char * key)
{
    pl.remove_key(key);
}
static inline void param_list_process_command_line_and_extra_parameter_files(cxx_param_list & pl,
        int * p_argc, char const *** p_argv)
{
    pl.process_command_line_and_extra_parameter_files(*p_argc, *p_argv);
}

static inline void param_list_process_command_line(cxx_param_list & pl,
        int * p_argc, char const *** p_argv, int accept_trailing_args)
{
    pl.process_command_line(*p_argc, *p_argv, accept_trailing_args);
}

static inline void param_list_display(cxx_param_list & pl,
        FILE * f)
{
    pl.display_debug(f);
}
static inline bool param_list_empty(cxx_param_list const & pl)
{
    return pl.empty();
}
static inline int param_list_warn_unused(cxx_param_list const & pl)
{
    return pl.warn_unused();
}
static inline int param_list_update_cmdline(cxx_param_list & pl,
    int * argc, char const *** argv)
{
    return pl.update_cmdline(*argc, *argv);
}
static inline size_t param_list_get_list_count(cxx_param_list & pl, const char * key)
{
    return pl.get_list_count(key);
}
using cado::params::collect_command_line;
using cado::params::parameter_error;





#endif /* CADO_PARAMS_HPP */
