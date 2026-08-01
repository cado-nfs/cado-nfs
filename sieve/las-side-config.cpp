#include "cado.h" // IWYU pragma: keep

#include <climits>

#include <algorithm>
#include <stdexcept>
#include <map>
#include <vector>
#include <string>

#include "params.hpp"
#include "las-side-config.hpp"

void siever_side_config::declare_usage(cxx_param_list & pl)
{
    pl.declare_usage("lim", "sieving bounds per side");
    pl.declare_usage("lpb", "large prime bounds per side, in bits");
    pl.declare_usage("mfb", "post-sieve cofactor bounds per side, in bits");
    pl.declare_usage("lambda", "post-sieve qualification multipliers per side");
    pl.declare_usage("powlim", "limits on powers sieved per side");
    pl.declare_usage("ncurves", "cofactoring effort (number of curves) per side");
    pl.declare_usage("fb", "factor base files per side");
}

void batch_side_config::declare_usage(cxx_param_list & pl)
{
    pl.declare_usage("batchfile", "prime product tree files per side");
    pl.declare_usage("batchmfb", "cofactor bounds, per side, to be considered after batch cofactorization. After primes below 2^batchlpbX have been extracted, cofactors below this bound will go through ecm. Defaults to lpbX.");
    pl.declare_usage("batchlpb", "large prime bounds, per side, to be considered by batch cofactorization. Primes between limX and 2^batchlpbX will be extracted by product trees. Defaults to lpbX.");
}

#define DISPATCH_PARAMETER_COPY(type, name, argname, dfl) do {              \
        std::vector<type> t;                                                \
        bool const r = pl.parse_per_side(argname, t, nb_polys,              \
                            cado::params::copy_previous_side());            \
        if (!r) t.assign(nb_polys, dfl);                                    \
        res[#name] = r ? nb_polys : 0;                                      \
        if (!r && std::ranges::find(mo, argname) != mo.end())               \
            pl.fail("Parameter " argname " is mandatory (for all sides)");  \
        for(int i = 0 ; i < nb_polys ; i++)                                 \
            v[i].name = t[i];                                               \
    } while (0)

#define DISPATCH_PARAMETER_DEFAULT(type, name, argname, dfl) do {           \
        std::vector<type> t;                                                \
        bool const r = pl.parse_per_side(argname, t, nb_polys, type(dfl));  \
        if (!r) t.assign(nb_polys, dfl);                                    \
        res[#name] = r ? nb_polys : 0;                                      \
        if (!r && std::ranges::find(mo, argname) != mo.end())               \
            pl.fail("Parameter " argname " is mandatory (for all sides)");  \
        for(int i = 0 ; i < nb_polys ; i++)                                 \
            v[i].name = t[i];                                               \
    } while (0)

std::map<std::string, int> siever_side_config::parse(cxx_param_list & pl, std::vector<siever_side_config> & v, int nb_polys,
    std::vector<std::string> const & mandatory)
{
    std::map<std::string, int> res;
    v.assign(nb_polys, {});
    auto const & mo = mandatory;

    DISPATCH_PARAMETER_COPY(unsigned long, lim, "lim", ULONG_MAX);
    DISPATCH_PARAMETER_COPY(unsigned long, powlim, "powlim", ULONG_MAX);
    DISPATCH_PARAMETER_COPY(unsigned int, lpb, "lpb", UINT_MAX);
    DISPATCH_PARAMETER_COPY(unsigned int, mfb, "mfb", UINT_MAX);
    DISPATCH_PARAMETER_COPY(int, ncurves, "ncurves", -1);
    DISPATCH_PARAMETER_COPY(double, lambda, "lambda", 0);
    DISPATCH_PARAMETER_DEFAULT(std::string, fbfilename, "fb", "");
    return res;
}

void siever_side_config::lookup_parameters(cxx_param_list & pl, int nsides)
{
    std::vector<siever_side_config> v;
    siever_side_config::parse(pl, v, nsides);
}

std::map<std::string, int> batch_side_config::parse(cxx_param_list & pl, std::vector<batch_side_config> & v, int nb_polys,
    std::vector<std::string> const & mandatory)
{
    std::map<std::string, int> res;
    v.assign(nb_polys, {});
    auto const & mo = mandatory;

    // for(auto & s: v) s.batchmfb = s.batchlpb = s.lpb;
    DISPATCH_PARAMETER_COPY(unsigned int, batchlpb, "batchlpb", UINT_MAX);
    DISPATCH_PARAMETER_COPY(unsigned int, batchmfb, "batchmfb", UINT_MAX);
    DISPATCH_PARAMETER_DEFAULT(std::string, batchfilename, "batchfile", "");
    return res;
}

void batch_side_config::lookup_parameters(cxx_param_list & pl, int nsides)
{
    std::vector<batch_side_config> v;
    batch_side_config::parse(pl, v, nsides);
}

#undef DISPATCH_PARAMETER
