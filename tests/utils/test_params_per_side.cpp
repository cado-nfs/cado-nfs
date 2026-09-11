#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <string>
#include <vector>

#include "params.hpp"

/* Regression test for per-side parameters that are specified for only
 * some of the sides.
 *
 * -fb0 without -fb1 used to be dropped silently: parse_per_side_base()
 * recorded only whether the *last* side had been given, so a value
 * supplied for an earlier side was discarded and replaced by the default,
 * and no "unused parameter" warning was issued either, because the
 * argument had been consumed. This broke, among others, the documented
 * DLP-240 sieving command line, which sieves on side 0 only.
 */

static int failures = 0;

static void check(std::string const & what,
        std::string const & got, std::string const & expected)
{
    if (got != expected) {
        fprintf(stderr, "FAIL %s: got \"%s\", expected \"%s\"\n",
                what.c_str(), got.c_str(), expected.c_str());
        failures++;
    }
}

int main(int argc, char const * argv[])
{
    cxx_param_list pl;

    pl.declare_usage("fb", "factor base files per side");
    pl.declare_usage("batch", "batch limits per side");
    pl.declare_usage("lim", "sieving bounds per side");

    pl.process_command_line(argc, argv);

    /* default-filled: sides that were not given get the default */
    {
        std::vector<std::string> v;
        bool const r = pl.parse_per_side("fb", v, 2, std::string(""));
        check("fb parsed", r ? "yes" : "no", "yes");
        if (r) {
            check("fb[0]", v[0], "zero");
            check("fb[1]", v[1], "");
        }
    }

    {
        std::vector<std::string> v;
        bool const r = pl.parse_per_side("batch", v, 2, std::string(""));
        check("batch parsed", r ? "yes" : "no", "yes");
        if (r) {
            check("batch[0]", v[0], "");
            check("batch[1]", v[1], "zero");
        }
    }

    /* copy-previous-side: side 1 inherits side 0 */
    {
        std::vector<unsigned long> v;
        bool const r = pl.parse_per_side("lim", v, 2,
                cado::params::copy_previous_side());
        check("lim parsed", r ? "yes" : "no", "yes");
        if (r) {
            check("lim[0]", std::to_string(v[0]), "17");
            check("lim[1]", std::to_string(v[1]), "17");
        }
    }

    if (failures)
        fprintf(stderr, "%d check(s) failed\n", failures);
    return failures ? 1 : 0;
}
