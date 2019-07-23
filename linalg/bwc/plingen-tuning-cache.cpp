#include "cado.h"
#include <sstream>
#include <cstdio>

#include "plingen-tuning-cache.hpp"
#include "macros.h"

#include <tuple>
#include <iostream>
#include <type_traits>

template <size_t n, typename... T>
typename std::enable_if<(n >= sizeof...(T)), std::ostream&>::type
    print_tuple(std::ostream& os, const std::tuple<T...>&)
{ return os; }

template <size_t n, typename... T>
typename std::enable_if<(n < sizeof...(T)), std::ostream&>::type
    print_tuple(std::ostream& os, const std::tuple<T...>& tup)
{
    if (n)
        os << ";";
    os << std::get<n>(tup);
    return print_tuple<n+1>(os, tup);
}
template <typename... T>
std::ostream& operator<<(std::ostream& os, const std::tuple<T...>& tup)
{
    return print_tuple<0>(os, tup);
}

template <size_t n, typename... T>
typename std::enable_if<(n >= sizeof...(T)), std::istream&>::type
    parse_tuple(std::istream& is, std::tuple<T...>&)
{ return is; }

template <size_t n, typename... T>
typename std::enable_if<(n < sizeof...(T)), std::istream&>::type
    parse_tuple(std::istream& is, std::tuple<T...>& tup)
{
    if (n) {
        char c = is.get();
        if (c != ';' && c != ' ') {
            is.setstate(std::ios_base::failbit);
            return is;
        }
    }
    is >> std::get<n>(tup);
    return parse_tuple<n+1>(is, tup);
}
template <typename... T>
std::istream& operator>>(std::istream& is, std::tuple<T...>& tup)
{
    return parse_tuple<0>(is, tup);
}

template<typename T, typename U>
std::istream& operator>>(std::istream& is, std::pair<T,U> & x)
{
    is >> x.first;
    char c = is.get();
    if (!is || c != ',') {
        /* since c++11, unget clears eofbit */
        is.setstate(std::ios_base::failbit);
        return is;
    }
    is >> x.second;
    return is;
}

template<typename T>
std::istream& operator>>(std::istream& is, std::list<T> & L)
{
    /* Reads a ***COMMA*** separatd list */
    L.clear();
    for(unsigned int item = 0 ; ; item++) {
        if (item) {
            char c = is.get();
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
std::ostream& operator<<(std::ostream& os, std::pair<T,U> const & x)
{
    os << x.first << ',' << x.second;
    return os;
}


template<typename T>
std::ostream& operator<<(std::ostream& os, std::list<T> const & L)
{
    unsigned int item=0;
    for(auto const & x : L) {
        if (item++) os << ',';
        os << x;
    }
    return os;
}

template <typename T, size_t n>
std::ostream& operator<<(std::ostream& os, const std::array<T, n>& arr)
{
    for(size_t i = 0 ; i < n ; i++) {
        if (i) os << ";";
        os << arr[i];
    }
    return os;
}
template <typename T, size_t n>
std::istream& operator>>(std::istream& is, std::array<T, n>& arr)
{
    for(size_t i = 0 ; i < n ; i++) {
        if (i) {
            char c;
            std::ios_base::fmtflags ff = is.flags();
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


std::istream& operator>>(std::istream& is, plingen_tuning_cache::basecase_key & K) {
    return is >> (plingen_tuning_cache::basecase_key::super&) K;
}
std::ostream& operator<<(std::ostream& os, plingen_tuning_cache::basecase_key const & K) {
    return os << (plingen_tuning_cache::basecase_key::super const&) K;
}
std::istream& operator>>(std::istream& is, plingen_tuning_cache::mul_key & K) {
    return is >> (plingen_tuning_cache::mul_key::super&) K;
}
std::ostream& operator<<(std::ostream& os, plingen_tuning_cache::mul_key const & K) {
    return os << (plingen_tuning_cache::mul_key::super const&) K;
}
std::istream& operator>>(std::istream& is, plingen_tuning_cache::mp_key & K) {
    return is >> (plingen_tuning_cache::mp_key::super&) K;
}
std::ostream& operator<<(std::ostream& os, plingen_tuning_cache::mp_key const & K) {
    return os << (plingen_tuning_cache::mp_key::super const&) K;
}

void plingen_tuning_cache::load(const char * timing_cache_filename)/*{{{*/
{
    if (timing_cache_filename == NULL) return;

    FILE * f = fopen(timing_cache_filename, "r");

    if (f == NULL && errno == ENOENT) {
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
            if (is >> K >> V)
                basecase_cache[K] = V;
        } else if (step == "MUL") {
            mul_key K;
            mul_value V;
            if (is >> K >> V)
                mul_cache[K] = V;
        } else if (step == "MP") {
            mp_key K;
            mp_value V;
            if (is >> K >> V)
                mp_cache[K] = V;
        } else {
            fprintf(stderr, "parse error in %s\nwhile reading line:\n%s\n",
                    timing_cache_filename, q);
        }
        if (!is) {
            fprintf(stderr, "parse error in %s\nwhile reading line:\n%s\n",
                    timing_cache_filename, q);
        }
    }
    fclose(f);
}/*}}}*/

void plingen_tuning_cache::save(const char * timing_cache_filename)/*{{{*/
{
    if (timing_cache_filename == NULL) return;

    FILE * f = fopen(timing_cache_filename, "w");
    ASSERT_ALWAYS(f);
    for(auto const & e : basecase_cache) {
        std::ostringstream os;
        os << "basecase" << ";" << e.first << ";" << e.second;
        fprintf(f, "%s\n", os.str().c_str());
    }
    for(auto const & e : mul_cache) {
        std::ostringstream os;
        os << "MUL" << ";" << e.first << ";" << e.second;
        fprintf(f, "%s\n", os.str().c_str());
    }
    for(auto const & e : mp_cache) {
        std::ostringstream os;
        os << "MP" << ";" << e.first << ";" << e.second;
        fprintf(f, "%s\n", os.str().c_str());
    }
    fclose(f);
}/*}}}*/
