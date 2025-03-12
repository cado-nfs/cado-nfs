#include "cado.h" // IWYU pragma: keep

#include <string>
#include <stdexcept>

#include "utils_cxx.hpp"
#include "macros.h"

decomposed_path::decomposed_path(std::string const & s)
{
    if (s.empty())
        throw std::invalid_argument("empty path");

    auto it = s.begin();
    if (*it == '/')
        emplace_back();
    for( ; it != s.end() && *it == '/' ; ++it);

    for( ; it != s.end() ; ) {
        auto it0 = it;
        for( ; it != s.end() && *it != '/' ; ++it);
        emplace_back(it0, it);
        for( ; it != s.end() && *it == '/' ; ++it);
    }
    ASSERT_ALWAYS(!empty());
}

std::string decomposed_path::dirname() const
{
    auto b = begin();
    auto e = end();
    --e;
    auto s = join(b, e, "/");
    if (s.empty())
        return ".";
    return s;
}

bool decomposed_path::is_relative() const {
    return !front().empty();
}

std::string decomposed_path::extension() const {
    auto const & s = back();
    if (s == "." || s == "..")
        return {};
    auto v = split(s, ".");
    if (v.size() == 2 && v.front().empty())
        return {};
    if (v.size() == 1)
        return {};
    return std::string(".").append(v.back());
}

decomposed_path::operator std::string() const {
    if (is_absolute() && size() == 1)
        return "/";
    return join(*this, "/");
}
