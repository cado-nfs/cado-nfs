#ifndef DECOMP_HPP
#define DECOMP_HPP

#include <cstdio>
#include <istream>
#include <ostream>
#include <vector>

#include "fmt/base.h"
#include "fmt/ostream.h"

struct decomp : public std::vector<unsigned int> {
    // nb_elem is the number of elements which satisfy this decomposition
    // I'm not sure it's a good idea, since this number is attached to a
    // target size, which we don't store here (it's _not_ directly
    // inferred by the summands in the std::vector<> parent).
    double nb_elem = 0;
    template <typename... Args>
    explicit decomp(double n, Args &&... args)
        : std::vector<unsigned int>(std::forward<Args>(args)...)
        , nb_elem(n)
    {
    }
    explicit decomp(double n, std::vector<unsigned int> v)
        : std::vector<unsigned int>(std::move(v))
        , nb_elem(n)
    {
    }
    decomp() = default;
    ~decomp() = default;
    decomp(decomp const &) = default;
    decomp(decomp &&) = default;
    decomp & operator=(decomp const &) = default;
    decomp & operator=(decomp &&) = default;

    bool operator<(decomp const &) const;
    bool operator==(decomp const & o) const { return !operator!=(o); }
    bool operator!=(decomp const & o) const { return *this < o || o < *this; }
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
