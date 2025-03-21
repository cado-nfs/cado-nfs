#ifndef TAB_DECOMP_HPP
#define TAB_DECOMP_HPP

#include <cstdio>

#include <vector>
#include <istream>

#include "fmt/base.h"

#include "decomp.hpp"

typedef std::vector<decomp> tabular_decomp;


std::ostream& operator<<(std::ostream& is, tabular_decomp const &);
std::istream& operator>>(std::istream& is, tabular_decomp &);

namespace fmt {
    template<>
    struct formatter<tabular_decomp>: ostream_formatter {};
}


#endif				/* TAB_DECOMP_HPP */
