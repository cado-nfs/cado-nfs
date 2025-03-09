#ifndef UTILS_ISTREAM_MATCHER_HPP_
#define UTILS_ISTREAM_MATCHER_HPP_

#include <string>
#include <istream>
#include <ios>
#include <iterator>

/* this is a hacky wrapper around istreams, so that we can do a very
 * simple scanf replacement for many of the fairly common use cases
 */
template<typename CharT, typename Traits = std::char_traits<CharT>>
struct istream_matcher {
    std::basic_istream<CharT, Traits> & r;
    explicit istream_matcher(std::basic_istream<CharT, Traits> & r) : r(r) {}
    operator bool() const { return (bool) r; }
    /* we can't conveniently overload >> std::ws for this type, probably
     * because it doesn't derive from std::ws as it should
     */
    istream_matcher& ws() { r >> std::ws; return *this; }
    bool eof() const { return r.eof(); }
    bool good() const { return r.good(); }
    bool bad() const { return r.bad(); }
    bool fail() const { return r.fail(); }
};

template<typename CharT, typename Traits, typename T>
inline istream_matcher<CharT, Traits> &operator >>(istream_matcher<CharT, Traits> &is, T & x) {
    is.r >> x;
    return is;
}


template<typename CharT, typename Traits, std::size_t N>
istream_matcher<CharT, Traits> &operator >>(istream_matcher<CharT, Traits> &is, char const (& literal)[N])
{
    auto it = std::end(literal);
    --it;
    if(!std::equal(std::begin(literal), it, std::istreambuf_iterator<CharT>{is.r})) {
        is.r.setstate(is.r.rdstate() | std::ios::failbit);
    }
    return is;
}

template<typename CharT, typename Traits>
istream_matcher<CharT, Traits> &operator >>(istream_matcher<CharT, Traits> &is, std::string const &literal)
{
    auto it = std::end(literal);
    --it;
    if(!std::equal(std::begin(literal), it, std::istreambuf_iterator<CharT>{is.r})) {
        is.r.setstate(is.r.rdstate() | std::ios::failbit);
    }
    return is;
}

#endif	/* UTILS_ISTREAM_MATCHER_HPP_ */
