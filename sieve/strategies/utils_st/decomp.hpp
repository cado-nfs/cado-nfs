#ifndef DECOMP_HPP
#define DECOMP_HPP

#include <cstdio>
#include <ostream>
#include <vector>

#include "fmt/base.h"
#include "fmt/ostream.h"

struct decomp : public std::vector<unsigned int> {
    double nb_elem = 0; // number of elements which satisfy this decomposition
    template <typename... Args>
    explicit decomp(double n, Args &&... args)
        : std::vector<unsigned int>(std::forward<Args>(args)...)
        , nb_elem(n)
    {
    }
    decomp() = default;
    ~decomp() = default;
    decomp(decomp const &) = default;
    decomp(decomp &&) = default;
    decomp & operator=(decomp const &) = default;
    decomp & operator=(decomp &&) = default;
};

static inline unsigned int
is_good_decomp(decomp const & D, unsigned int len_p_min, unsigned int len_p_max)
{
    for (auto x: D)
        if (x > len_p_max || x < len_p_min)
            return false;
    return true;
}

std::istream & operator>>(std::istream & is, decomp &);
std::ostream & operator<<(std::ostream & o, decomp const & D);

namespace fmt
{
template <> struct formatter<decomp> : ostream_formatter {
};
} // namespace fmt

#endif /* DECOMP_HPP */
