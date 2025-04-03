#ifndef CADO_UTILS_ISTREAM_MATCHER_HPP
#define CADO_UTILS_ISTREAM_MATCHER_HPP

#include <cstddef>

#include <algorithm>
#include <string>
#include <istream>
#include <ios>
#include <iterator>

/* this is a hacky wrapper around istreams, so that we can do a very
 * simple scanf replacement for many of the fairly common use cases
 */
template<typename CharT, typename Traits = std::char_traits<CharT>>
struct basic_istream_matcher {
    std::basic_istream<CharT, Traits> & r;
    explicit basic_istream_matcher(std::basic_istream<CharT, Traits> & r) : r(r) {}
    operator bool() const { return (bool) r; }
    /* we can't conveniently overload >> std::ws for this type, probably
     * because it doesn't derive from std::ws as it should
     */
    basic_istream_matcher& ws() { r >> std::ws; return *this; }
    bool eof() const { return r.eof(); }
    bool good() const { return r.good(); }
    bool bad() const { return r.bad(); }
    bool fail() const { return r.fail(); }
    int peek() const { return r.peek(); }
    int get() const { return r.get(); }
    void setstate(std::ios_base::iostate state) { r.setstate(state); }
};

typedef basic_istream_matcher<char> istream_matcher;

template<typename CharT, typename Traits, typename T>
inline basic_istream_matcher<CharT, Traits> &operator >>(basic_istream_matcher<CharT, Traits> &is, T & x) {
    is.r >> x;
    return is;
}


template<typename CharT, typename Traits, std::size_t N>
basic_istream_matcher<CharT, Traits> &operator >>(basic_istream_matcher<CharT, Traits> &is, char const (& literal)[N])
{
    auto it = std::end(literal);
    --it;
    if(!std::equal(std::begin(literal), it, std::istreambuf_iterator<CharT>{is.r})) {
        is.r.setstate(is.r.rdstate() | std::ios::failbit);
    }
    return is;
}

template<typename CharT, typename Traits>
basic_istream_matcher<CharT, Traits> &operator >>(basic_istream_matcher<CharT, Traits> &is, std::string const &literal)
{
    auto it = std::end(literal);
    --it;
    if(!std::equal(std::begin(literal), it, std::istreambuf_iterator<CharT>{is.r})) {
        is.r.setstate(is.r.rdstate() | std::ios::failbit);
    }
    return is;
}

/* This is an alternative approach. Maybe it's better */
struct expect_stream_separator {
    std::string s;
    explicit expect_stream_separator(char c)
        : s(1, c)
    {}
    explicit expect_stream_separator(std::string s)
        : s(std::move(s))
    {}
};

static inline std::istream& operator>>(std::istream & is, expect_stream_separator const & e)
{
    if(!std::equal(std::begin(e.s), std::end(e.s), std::istreambuf_iterator<char>{is})) {
        is.setstate(is.rdstate() | std::ios::failbit);
    }
    return is;
}


#endif	/* CADO_UTILS_ISTREAM_MATCHER_HPP */
