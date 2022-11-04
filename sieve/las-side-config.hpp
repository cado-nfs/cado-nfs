#ifndef LAS_SIDE_CONFIG_HPP_
#define LAS_SIDE_CONFIG_HPP_

#include "params.h"
#include <string>
#include <vector>
#include <map>

#define SIDE_CONFIG_COLLECTOR(C, lpb)				\
    static std::vector<decltype(C::lpb)>                        \
    CPP_PAD(collect_, lpb)                                      \
    (std::vector<C> const & A) {                                \
        std::vector<decltype(C::lpb)> res;		        \
        res.reserve(A.size());					\
        for(auto const & s : A) res.emplace_back(s.lpb);	\
        return res;						\
    }

struct siever_side_config {
    unsigned long lim;    /* factor base bound */
    unsigned long powlim; /* bound on powers in the factor base */

    unsigned int lpb;     /* large prime bound is 2^lpb */
    unsigned int mfb;     /* bound for residuals is 2^mfb */
    int ncurves;          /* number of cofactorization curves */
    double lambda;        /* lambda sieve parameter */

    std::string fbfilename;

    static void declare_usage(cxx_param_list &);
    static void lookup_parameters(cxx_param_list &, int nsides = 2);

    /* This parses a command line into a vector of siever_side_configs.
     * The returned value gives, per config key, the number of sides for
     * which config data was found */
    static std::map<std::string, int> parse(cxx_param_list &, std::vector<siever_side_config> &, int, std::vector<std::string> const & mandatory = std::vector<std::string>());

    SIDE_CONFIG_COLLECTOR(siever_side_config, lim)
    SIDE_CONFIG_COLLECTOR(siever_side_config, powlim)
    SIDE_CONFIG_COLLECTOR(siever_side_config, lpb)
    SIDE_CONFIG_COLLECTOR(siever_side_config, mfb)
    SIDE_CONFIG_COLLECTOR(siever_side_config, ncurves)
    SIDE_CONFIG_COLLECTOR(siever_side_config, lambda)
    // SIDE_CONFIG_COLLECTOR(siever_side_config, fbfilename)
};

struct batch_side_config {
    std::string batchfilename;
    unsigned int batchlpb;
    unsigned int batchmfb;

    static void declare_usage(cxx_param_list &);
    static void lookup_parameters(cxx_param_list &, int nsides = 2);

    /* This parses a command line into a vector of batch_side_configs.
     * The returned value gives, per config key, the number of sides for
     * which config data was found */
    static std::map<std::string, int> parse(cxx_param_list &, std::vector<batch_side_config> &, int, std::vector<std::string> const & mandatory = std::vector<std::string>());

    SIDE_CONFIG_COLLECTOR(batch_side_config, batchlpb)
    SIDE_CONFIG_COLLECTOR(batch_side_config, batchmfb)
    SIDE_CONFIG_COLLECTOR(batch_side_config, batchfilename)
};

#undef SIDE_CONFIG_COLLECTOR


#endif	/* LAS_SIDE_CONFIG_HPP_ */
