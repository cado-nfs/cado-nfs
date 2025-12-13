#include "cado.h" // IWYU pragma: keep

#include <ostream>
#include <istream>
#include <algorithm>

#include "decomp.hpp"
#include "utils_cxx.hpp"

std::ostream & operator<<(std::ostream & o, decomp const & D)
{
    return o << "[ " << join(D, " ") << " ]: " << D.nb_elem;
}

std::istream & operator>>(std::istream & is, decomp & D)
{
    D.clear();
    is >> expect("[");
    for(unsigned int v ; is.good() && is >> std::ws && is.peek() != ']' ; ) {
        is >> v;
        D.push_back(v);
    }
    return is >> expect("]:") >> D.nb_elem;
}

bool decomp::operator<(decomp const & o) const
{
    if (!std::ranges::is_sorted(*this)) {
        decomp a = *this;
        std::ranges::sort(a);
        return a < o;
    }
    if (!std::ranges::is_sorted(o)) {
        decomp b = o;
        std::ranges::sort(b);
        return *this < b;
    }
    bool c = std::ranges::lexicographical_compare(*this, o);
    if (size() != o.size())
        return c;
    if (!std::ranges::equal(*this, o))
        return c;
    return nb_elem < o.nb_elem - 0.1;
}
