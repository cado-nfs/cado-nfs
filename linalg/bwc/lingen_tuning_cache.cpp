#include "cado.h" // IWYU pragma: keep

#include <cctype>
#include <cerrno>
#include <cstdio>

#include <array>
#include <ios>
#include <istream>
#include <list>
#include <ostream>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>

#include "lingen_tuning_cache.hpp"
#include "macros.h"

template <size_t n, typename... T>
static std::ostream& print_tuple(std::ostream& os, const std::tuple<T...>&)
    requires (n >= sizeof...(T))
{ return os; }

template <size_t n, typename... T>
static std::ostream& print_tuple(std::ostream& os, const std::tuple<T...>& tup)
    requires (n < sizeof...(T))
{
    if (n)
        os << ";";
    os << std::get<n>(tup);
    return print_tuple<n+1>(os, tup);
}
template <typename... T>
static std::ostream& operator<<(std::ostream& os, const std::tuple<T...>& tup)
{
    return print_tuple<0>(os, tup);
}

template <size_t n, typename... T>
static std::istream& parse_tuple(std::istream& is, std::tuple<T...>&)
    requires (n >= sizeof...(T))
{ return is; }

template <size_t n, typename... T>
static std::istream& parse_tuple(std::istream& is, std::tuple<T...>& tup)
    requires (n < sizeof...(T))
{
    if (!is) return is;
    if (n) {
        int const c = is.get();
        if (c != ';' && c != ' ') {
            is.setstate(std::ios_base::failbit);
            return is;
        }
    }
    is >> std::get<n>(tup);
    return parse_tuple<n+1>(is, tup);
}
template <typename... T>
static std::istream& operator>>(std::istream& is, std::tuple<T...>& tup)
{
    return parse_tuple<0>(is, tup);
}

template<typename T, typename U>
static std::istream& operator>>(std::istream& is, std::pair<T,U> & x)
{
    is >> x.first;
    int const c = is.get();
    if (!is || c != ',') {
        /* since c++11, unget clears eofbit */
        is.setstate(std::ios_base::failbit);
        return is;
    }
    is >> x.second;
    return is;
}

template<typename T>
static std::istream& operator>>(std::istream& is, std::list<T> & L)
{
    /* Reads a ***COMMA*** separatd list */
    L.clear();
    for(unsigned int item = 0 ; ; item++) {
        if (item) {
            int const c = is.get();
            if (!is || c != ',') {
                /* since c++11, unget clears eofbit */
                is.unget();
                return is;
            }
        }
        T x;
        is >> x;
        L.emplace_back(std::move(x));
    }
}

template<typename T, typename U>
static std::ostream& operator<<(std::ostream& os, std::pair<T,U> const & x)
{
    os << x.first << ',' << x.second;
    return os;
}


template<typename T>
static std::ostream& operator<<(std::ostream& os, std::list<T> const & L)
{
    unsigned int item=0;
    for(auto const & x : L) {
        if (item++) os << ',';
        os << x;
    }
    return os;
}

template <typename T, size_t n>
static std::ostream& operator<<(std::ostream& os, const std::array<T, n>& arr)
{
    for(size_t i = 0 ; i < n ; i++) {
        if (i) os << ";";
        os << arr[i];
    }
    return os;
}
template <typename T, size_t n>
static std::istream& operator>>(std::istream& is, std::array<T, n>& arr)
{
    for(size_t i = 0 ; i < n ; i++) {
        if (i) {
            char c;
            std::ios_base::fmtflags const ff = is.flags();
            is.flags(ff & ~std::ios_base::skipws);
            is >> c;
            is.flags(ff);
            if (c != ';' && c != ' ') {
                is.setstate(std::ios_base::failbit);
                return is;
            }
        }
        is >> arr[i];
    }
    return is;
}


void lingen_tuning_cache::load(const char * timing_cache_filename)/*{{{*/
{
    if (timing_cache_filename == nullptr) return;

    FILE * f = fopen(timing_cache_filename, "r");

    if (f == nullptr && errno == ENOENT) {
        printf("# %s: no cache file found\n", timing_cache_filename);
        return;
    }

    printf("# reading timing cache from %s\n", timing_cache_filename);

    ASSERT_ALWAYS(f);
    for(;;) {
        char line[2048];
        if (!fgets(line, sizeof(line), f)) break;

        if (line[0] == '#') continue;
        char * q = line;
        for( ; *q && isspace(*q) ; q++);
        if (!*q) continue;

        for(char * p = q; *p ; p++)
            if (*p == ';') *p=' ';

        std::istringstream is(q);

        std::string step;

        is >> step;

        if (step == "basecase") {
            basecase_key K;
            basecase_value V;
            if (is >> K.np >> K.m >> K.n >> K.length >> K.nthreads >> V)
                basecase_cache[K] = V;
        } else {
            mul_or_mp_key K;
            bool good = false;
            if (step == op_mul_or_mp_base::op_name(op_mul_or_mp_base::OP_MP)) {
                K.op_type = op_mul_or_mp_base::OP_MP;
                good = true;
            } else if (step == op_mul_or_mp_base::op_name(op_mul_or_mp_base::OP_MUL)) {
                K.op_type = op_mul_or_mp_base::OP_MUL;
                good = true;
            } else {
                is.setstate(std::ios_base::failbit);
            }
            if (good) {
                lingen_substep_schedule::fft_type_unserialize(is, K.fft_type);
                mul_or_mp_value V;
                if (is >> K.np >> K.na >> K.nb >> V) {
                    mul_or_mp_cache[K] = V;
                }
            }
        }
        if (!is) {
            fprintf(stderr, "parse error in %s\nwhile reading line:\n%s\n",
                    timing_cache_filename, q);
        }
    }
    fclose(f);
}/*}}}*/

void lingen_tuning_cache::save(const char * timing_cache_filename)/*{{{*/
{
    if (timing_cache_filename == nullptr) return;

    FILE * f = fopen(timing_cache_filename, "w");
    ASSERT_ALWAYS(f);
    for(auto const & e : basecase_cache) {
        std::ostringstream os;
        os << "basecase"
            << ";" << e.first.np
            << ";" << e.first.m
            << ";" << e.first.n
            << ";" << e.first.length
            << ";" << e.first.nthreads
            << ";" << e.second;
        fprintf(f, "%s\n", os.str().c_str());
    }
    for(auto const & e : mul_or_mp_cache) {
        std::ostringstream os;
        os << op_mul_or_mp_base::op_name(e.first.op_type)
            << ";" << lingen_substep_schedule::fft_name(e.first.fft_type)
            << ";" << e.first.np
            << ";" << e.first.na
            << ";" << e.first.nb
            << ";" << e.second;
        fprintf(f, "%s\n", os.str().c_str());
    }
    fclose(f);
}/*}}}*/
