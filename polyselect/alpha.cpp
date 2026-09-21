#include "cado.h" // IWYU pragma: keep

#include <cstdlib>

#include <string>

#include "fmt/base.h"

#include "cado_poly.hpp"
#include "macros.h"
#include "params.hpp"
#include "polyselect_alpha.h"
#include "cado_main.hpp"

struct alpha_process {
    parameter_mandatory<std::string, "poly", "input polynomial file">
        polyfilename;
    parameter_with_default<unsigned long , "B",
        "bound on primes used for alpha computation",
        CADO_STRINGIZE(ALPHA_BOUND)>
        B;
    cxx_cado_poly cpoly;

    static void configure(cxx_param_list & pl)
    {
        pl.declare_usage_header("Compute alpha of the given polynomial\n");
        decltype(polyfilename)::configure(pl);
        decltype(B)::configure(pl);
    }

    explicit alpha_process(cxx_param_list & pl)
        : polyfilename(pl)
        , B(pl)
    {
        if (!cpoly.read(polyfilename.parameter_value().c_str()))
        {
            pl.fail("Error reading polynomial file\n");
        }
    }

    void run()
    {
        for (size_t i = 0; auto const & poly: cpoly) {
            double alpha = get_alpha(poly, B);
            double alpha_proj = get_alpha_projective(poly, B);
            fmt::print("alpha(poly{}) = {:.2f} (proj: {:.2f})\n",
                       i, alpha, alpha_proj);
            ++i;
        }
    }
};

static int main_(int argc, char const * argv[]);

int main(int argc, char const * argv[])
{
    return cado::main_wrapper(main_, argc, argv);
}

static int main_(int argc, char const * argv[])
{
    cxx_param_list pl;

    alpha_process::configure(pl);
    pl.process_command_line(argc, argv, false);

    pl.print_command_line(stdout);
    fflush(stdout);

    alpha_process alpha(pl);

    if (pl.warn_unused())
        pl.fail("Error, unused parameters are given\n");

    alpha.run();

    return EXIT_SUCCESS;
}
