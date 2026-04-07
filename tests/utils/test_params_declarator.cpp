#include "cado.h" // IWYU pragma: keep

#include <cstdio>

#include <stdexcept>

#include "params.hpp"

struct puppet {
    /* 1. and 2. */
    parameter<std::string, "foo", "give the h2g2 number", "forty-two"> foo;
    parameter_mandatory<int, "bar", "number of repeats"> bar;
    parameter_switch<"nl", "add newlines"> add_newlines;

    /* 3. and 4. */
    static void configure(cxx_param_list & pl)
    {
        pl.declare_usage_section("General operational flags");

        decltype(puppet::foo)::configure(pl);
        decltype(puppet::bar)::configure(pl);
        decltype(puppet::add_newlines)::configure(pl);
    }

    /* 5. and 6. */
    explicit puppet(cxx_param_list & pl)
        : foo(pl)
        , bar(pl)
        , add_newlines(pl)
    {
        if (bar <= 0)
            throw std::runtime_error("must have bar>0");
    }

    void doit() const
    {
        for(int i = 0 ; i < bar ; i++) {
            printf("%d\n", foo.parameter_value());
            if (add_newlines)
                printf("\n");
        }
    }
};

int main(int argc, char const * argv[])
{
    cxx_param_list pl;

    /* These two calls are here so that the output of --help looks nicer */
    pl.declare_usage_header("This is an example program\n");

    /* 3. and 4. */
    puppet::configure(pl);

    /* most of the stuff happens here */
    pl.process_command_line(argc, argv);

    const puppet p(pl);

    /* this will trap errors */
    if (pl.warn_unused())
        throw std::runtime_error("not all parameters were parsed");


    p.doit();

    return 0;
}
