#include "cado.h" // IWYU pragma: keep

#include <cstdio>

#include <stdexcept>

#include "params.h"

struct puppet {
    /* 1. and 2. */
    std::string foo = "forty-two";
    int bar;
    bool add_newlines = false;

    /* 3. */
    static void declare_usage(cxx_param_list & pl)
    {
        pl.declare_usage_section("General operational flags");

        pl.declare_usage("foo", "give the h2g2 number (default forty-two)");
        pl.declare_usage("bar", "number of repeats");
        pl.declare_usage("nl", "add newlines");
    }

    /* 4. */
    static void configure_switches(cxx_param_list & pl)
    {
        pl.configure_switch("nl");
    }

    /* 5. and 6. */
    explicit puppet(cxx_param_list & pl)
        : foo(pl.parse_optional<std::string>("foo", "forty-two"))
        , bar(pl.parse_mandatory<int>("bar"))
        , add_newlines(pl.parse<bool>("nl"))
    {
        /* as an alternative to using initializers, we could also do:
        pl.parse("foo", foo);
        pl.parse_mandatory("bar", bar);
         */
        if (bar <= 0)
            throw std::runtime_error("must have bar>0");
    }

    void doit() const
    {
        for(int i = 0 ; i < bar ; i++) {
            printf("%s\n", foo.c_str());
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
    puppet::declare_usage(pl);
    puppet::configure_switches(pl);

    /* most of the stuff happens here */
    pl.process_command_line(argc, argv);

    const puppet p(pl);

    /* this will trap errors */
    if (pl.warn_unused())
        throw std::runtime_error("not all parameters were parsed");


    p.doit();

    return 0;
}
