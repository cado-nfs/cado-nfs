#include "cado.h"       // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include <ios>
#include <string>
#include <istream>
#include <sstream>
#include <iostream>

#include "fmt/base.h"
#include "number_literal.hpp"

static bool test_tokenize(std::string const & s)
{
    std::istringstream is(s);
    cado::number_literal N;
    if (!(is >> N))
        return false;
    is >> std::ws;
    return is.eof();
}

int main()
{
    setvbuf(stderr, nullptr, _IONBF, 0);
    setvbuf(stdout, nullptr, _IONBF, 0);

    bool ok = true;

    for(auto const & s : {
            "1",
            "0",
            "-1",
            "0.0",
            "-inf",
            "nan", "NaN",
            "-1.23",
            "-1.23e2",
            "-0x1.deadbeefp+10",
            "-1i",
            "nan i",
            }) {
        const bool b = test_tokenize(s);
        fmt::print(stdout, "{} -> {}: {}\n", s, b, b ? "ok" : "NOK NOK");
        ok = ok && b;
    }
    for(auto const & s : {
            "1p",
            "1.1deadbeef",
            }) {
        const bool b = test_tokenize(s);
        fmt::print(stdout, "{} -> {}: {}\n", s, b, b ? "NOK NOK" : "ok");
        ok = ok && !b;
    }
    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
