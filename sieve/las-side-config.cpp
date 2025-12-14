#include "cado.h" // IWYU pragma: keep

#include <climits>

#include <algorithm>
#include <stdexcept>
#include <map>
#include <vector>
#include <string>

#include "params.h"
#include "las-side-config.hpp"

void siever_side_config::declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "lim", "sieving bounds per side");
    param_list_decl_usage(pl, "lpb", "large prime bounds per side, in bits");
    param_list_decl_usage(pl, "mfb", "post-sieve cofactor bounds per side, in bits");
    param_list_decl_usage(pl, "lambda", "post-sieve qualification multipliers per side");
    param_list_decl_usage(pl, "powlim", "limits on powers sieved per side");
    param_list_decl_usage(pl, "ncurves", "cofactoring effort (number of curves) per side");
    param_list_decl_usage(pl, "fb", "factor base files per side");
}

void batch_side_config::declare_usage(cxx_param_list & pl)
{
    param_list_decl_usage(pl, "batchfile", "prime product tree files per side");
    param_list_decl_usage(pl, "batchmfb", "cofactor bounds, per side, to be considered after batch cofactorization. After primes below 2^batchlpbX have been extracted, cofactors below this bound will go through ecm. Defaults to lpbX.");
    param_list_decl_usage(pl, "batchlpb", "large prime bounds, per side, to be considered by batch cofactorization. Primes between limX and 2^batchlpbX will be extracted by product trees. Defaults to lpbX.");
}

#define DISPATCH_PARAMETER2(type, name, argname, dfl, policy) do {	\
        std::vector<type> t(nb_polys, dfl);		                \
        const int r = param_list_parse_per_side<type>(pl,		\
                    argname, t.data(), nb_polys, policy);		\
        res[#name] = r;                                                 \
        if (!r && std::ranges::find(mo, argname) != mo.end())	\
            throw std::runtime_error("Parameter " argname		\
                    " is mandatory (for all sides)");			\
        for(int i = 0 ; i < nb_polys ; i++)				\
            v[i].name = t[i];						\
    } while (0)

#define DISPATCH_PARAMETER(type, name, dfl, policy)     \
            DISPATCH_PARAMETER2(type, name, #name, dfl, policy)

std::map<std::string, int> siever_side_config::parse(cxx_param_list & pl, std::vector<siever_side_config> & v, int nb_polys,
    std::vector<std::string> const & mandatory)
{
    std::map<std::string, int> res;
    v.assign(nb_polys, {});
    auto const & mo = mandatory;


    DISPATCH_PARAMETER(unsigned long, lim, ULONG_MAX, ARGS_PER_SIDE_DEFAULT_COPY_PREVIOUS);
    DISPATCH_PARAMETER(unsigned long, powlim, ULONG_MAX, ARGS_PER_SIDE_DEFAULT_COPY_PREVIOUS);
    DISPATCH_PARAMETER(unsigned int, lpb, UINT_MAX, ARGS_PER_SIDE_DEFAULT_COPY_PREVIOUS);
    DISPATCH_PARAMETER(unsigned int, mfb, UINT_MAX, ARGS_PER_SIDE_DEFAULT_COPY_PREVIOUS);
    DISPATCH_PARAMETER(int, ncurves, -1, ARGS_PER_SIDE_DEFAULT_COPY_PREVIOUS);
    DISPATCH_PARAMETER(double, lambda, 0, ARGS_PER_SIDE_DEFAULT_COPY_PREVIOUS);
    DISPATCH_PARAMETER2(std::string, fbfilename, "fb", "", ARGS_PER_SIDE_DEFAULT_AS_IS);
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
    DISPATCH_PARAMETER(unsigned int, batchlpb, UINT_MAX, ARGS_PER_SIDE_DEFAULT_COPY_PREVIOUS);
    DISPATCH_PARAMETER(unsigned int, batchmfb, UINT_MAX, ARGS_PER_SIDE_DEFAULT_COPY_PREVIOUS);
    DISPATCH_PARAMETER2(std::string, batchfilename, "batchfile", "", ARGS_PER_SIDE_DEFAULT_AS_IS);
    return res;
}

void batch_side_config::lookup_parameters(cxx_param_list & pl, int nsides)
{
    std::vector<batch_side_config> v;
    batch_side_config::parse(pl, v, nsides);
}

#undef DISPATCH_PARAMETER
