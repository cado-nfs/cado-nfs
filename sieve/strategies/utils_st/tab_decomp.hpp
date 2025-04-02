#ifndef TAB_DECOMP_HPP
#define TAB_DECOMP_HPP

#include <cstdio>

#include <istream>
#include <ostream>
#include <vector>

#include "fmt/base.h"

#include "decomp.hpp"

typedef std::vector<decomp> tabular_decomp;

std::ostream & operator<<(std::ostream & is, tabular_decomp const &);
std::istream & operator>>(std::istream & is, tabular_decomp &);

namespace fmt
{
template <> struct formatter<tabular_decomp> : ostream_formatter {
};
} // namespace fmt

#endif /* TAB_DECOMP_HPP */
