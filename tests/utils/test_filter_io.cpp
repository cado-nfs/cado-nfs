#include "cado.h"       // IWYU pragma: keep

#include <string>

#include "fmt/base.h"

#include "filter_io.hpp"
#include "utils_cxx.hpp"

int main()
{
    for(std::string const s :
            { "11,2a:3,4/-,5,8/-s+2ss,b,10/s^4",
              "11,2a:3,4/-1,5,8/-s+2*s*s,b,10/s^2*s^2\n",
              "11,2a:3/,4/-s^0,5,8/-1*s+s*s,8/s*s,b,10/s^2*s^2",
              })
    {
        using R = tnfs_indexed_relation;
        auto it = s.begin();
        R rel;
        int c = rel.parse(it);
        fmt::print("parser input: {}\n", s);
        fmt::print("parser output: c='0x{:02x}', a={}, b={}, #primes={}\n",
                c, rel.a, rel.b, rel.primes.size());
        for(auto const & pe : rel.primes) {
            fmt::print(" p={}, e=[{}]\n",
                    pe.p_or_h(),
                    join(pe.e.begin(), pe.e.end(), ", "));
        }
        fmt::print("{}\n", rel);
    }

    return 0;
}

