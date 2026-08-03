#ifndef SIEVE_RELATION_CACHE_HPP_
#define SIEVE_RELATION_CACHE_HPP_

#include <string>
#include <vector>

#include "params.hpp"
#include "cxx_mpz.hpp"

class relation_cache {
    parameter<std::string,
        "relation_cache",
        "Directory with cache of collected relation for sampling within a known data set. Useful only with --random-sample">
        path;
    std::vector<unsigned long> splits;
    std::string subdir_name(std::vector<unsigned long> const & split_q) const;
public:
    bool active() const { return !path().empty(); }
    static void configure(cxx_param_list & pl)
    {
        decltype(path)::configure(pl);
    }
    explicit relation_cache(cxx_param_list & pl);
    std::string find_filepath(cxx_mpz const & q) const;
};

#endif	/* SIEVE_RELATION_CACHE_HPP_ */
