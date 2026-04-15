#include "cado.h"       // IWYU pragma: keep

#include <string>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <iostream>

#include "cado_subrange_streambuf.hpp"

int main()
{
    std::string s = "The quick brown fox jumps over the lazy dog";

    std::istringstream is(s);

    cado::subrange_streambuf X(is.rdbuf(), 10, 19);

    std::istream isx(&X);

    std::string ex;

    isx >> std::noskipws;
    std::copy(std::istream_iterator<char>(isx), std::istream_iterator<char>(), std::back_inserter(ex));

    std::cout << ex << "\n";
    std::cout << ex.size() << "\n";

    return ex == "brown fox" ? EXIT_SUCCESS : EXIT_FAILURE;
}

