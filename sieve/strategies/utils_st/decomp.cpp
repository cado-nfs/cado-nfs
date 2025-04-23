#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include <ostream>
#include <istream>
#include <algorithm>

#include "decomp.hpp"
#include "utils_cxx.hpp"
#include "istream_matcher.hpp"

std::ostream & operator<<(std::ostream & o, decomp const & D)
{
    return o << "[ " << join(D, " ") << " ]: " << D.nb_elem;
}

std::istream & operator>>(std::istream & is, decomp & D)
{
    D.clear();
    is >> expect_stream_separator('[');
    for(unsigned int v ; is.good() && is >> std::ws && is.peek() != ']' ; ) {
        is >> v;
        D.push_back(v);
    }
    return is >> expect_stream_separator("]:") >> D.nb_elem;
}

bool decomp::operator<(decomp const & o) const
{
    if (!std::is_sorted(begin(), end())) {
        decomp a = *this;
        std::sort(a.begin(), a.end());
        return a < o;
    }
    if (!std::is_sorted(o.begin(), o.end())) {
        decomp b = o;
        std::sort(b.begin(), b.end());
        return *this < b;
    }
    bool c = std::lexicographical_compare(begin(), end(), o.begin(), o.end());
    if (size() != o.size())
        return c;
    if (!std::equal(begin(), end(), o.begin()))
        return c;
    return nb_elem < o.nb_elem - 0.1;
}
