#include "cado.h"
#include <iostream>
#include <iterator>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include "indexed_relation.hpp"
#include "omp_proxy.h"

/* This binary can be used on files in indexed form (e.g. out of dup2 or
 * out of purge) to recover the table of per-ideal density. This is just
 * a handy tool to have, but we don't make any particular use of it at
 * the moment.
 */

/* Note that this is slow, because indexed_relation parsing is slow -- it
 * was never meant to be optimized in the first place, since the fast
 * parsing is done by filter_io instead. (I'm curious as to the reason of
 * the slow processing, though).
 *
 * having the thing openmp'ed alleviates the concerns a little bit.
 *
 */
int main()
{
    std::vector<unsigned int> density;
#pragma omp parallel
    {
#pragma omp single
        {
            for(std::string s ; std::getline(std::cin, s) ; ) {
#pragma omp task
                {
                    std::istringstream is(s);
                    indexed_relation rel;
                    is >> rel;
                    rel.sort();
                    for(auto it = rel.data.begin() ; it != rel.data.end() ; ) {
                        int n = 0;
                        auto itx = it;
                        for(; itx != rel.data.end() && *it == *itx ; ++itx, n++);
                        index_t p = *it;
#pragma omp critical
                        {
                            if (p >= density.size())
                                std::fill_n(std::back_inserter(density), p + 1 - density.size(), 0);
                            density[p]++;
                        }
                        it = itx;
                    }
                }
            }
#pragma omp taskwait
        }
    }
    for(index_t p = 0 ; p < density.size() ; p++)
        std::cout << p << " " << density[p] << "\n";

    return 0;
}
