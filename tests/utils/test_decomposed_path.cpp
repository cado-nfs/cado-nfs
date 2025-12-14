#include "cado.h" // IWYU pragma: keep

#include <string>

#include "fmt/base.h"
#include "fmt/format.h"

#include "utils_cxx.hpp"
#include "macros.h"

int main()
{
    for(std::string s : {
            "/", "/etc/fstab", "./foo",
            "/usr/", "/usr", "/usr/bin/vim",
            "file.cpp"
            })
    {
        const decomposed_path D(s);
        auto quote = [](std::string const & s) {
            return fmt::format("\"{}\"", s);
        };
        fmt::print("path: \"{}\" -> {}, components: {}\n",
                s,
                D.is_relative() ? "relative" : "absolute",
                join(D, ",", quote));
        /* We don't make a distinction between paths based on whether
         * they're null-terminated or not.
         */
        ASSERT_ALWAYS(std::string(D) == s || (s.back() == '/'));
    }

    fmt::print("\n");

    for(const char * s : {
            "/foo/bar.txt",
            "/foo/bar.",
            "/foo/bar",
            "/foo/bar.txt/bar.cc",
            "/foo/bar.txt/bar.",
            "/foo/bar.txt/bar",
            "/foo/.",
            "/foo/..",
            "/foo/.hidden",
            "/foo/..bar",
            "/foo/..",
            "/foo/.profile.bak",
            })
    {
        const decomposed_path D(s);
        fmt::print("path: \"{}\", extension: \"{}\"\n", s, D.extension());
    }

    fmt::print("\n");

    for (const char * p : {
            "/foo/bar.txt", "/foo/.bar", "foo.bar.baz.tar",
            "/one/two/three/four", "five/six/seven/eight"
            }) {
        const decomposed_path D(p);
        fmt::print("path: \"{}\", dirname: \"{}\", basename: \"{}\"\n",
                std::string(D), D.dirname(), D.basename());
    }

    fmt::print("\n");

    return 0;
}
