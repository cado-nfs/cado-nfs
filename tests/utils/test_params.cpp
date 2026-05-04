#include "cado.h" // IWYU pragma: keep

#include <cstdio>

#include <stdexcept>

#include "params.hpp"

int main(int argc, char const * argv[])
{
    /* 1. and 2. */
    std::string foo = "forty-two";
    int bar;
    bool add_newlines = false;
    cxx_param_list pl;

    /* 3. */
    /* These two calls are here so that the output of --help looks nicer */
    pl.declare_usage_header("This is an example program\n");
    pl.declare_usage_section("General operational flags");

    pl.declare_usage("foo", "give the h2g2 number (default forty-two)");
    pl.declare_usage("bar", "number of repeats");
    pl.declare_usage("nl", "add newlines");

    /* 4. */
    pl.configure_switch("nl");

    /* most of the stuff happens here */
    pl.process_command_line(argc, argv);

    /* 5. */
    pl.parse("foo", foo);
    pl.parse_mandatory("bar", bar);
    pl.parse("nl", add_newlines);

    /* this will trap errors */
    if (pl.warn_unused())
        throw std::runtime_error("not all parameters were parsed");

    /* 6. */
    if (bar <= 0)
        throw std::runtime_error("must have bar>0");

    for(int i = 0 ; i < bar ; i++) {
        printf("%s\n", foo.c_str());
        if (add_newlines)
            printf("\n");
    }

    return 0;
}
