#include "cado.h"
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include "indexed_relation.hpp"

int main()
{
    std::vector<unsigned int> density;
    for(indexed_relation rel ; std::cin >> rel ; ) {
        rel.sort();
        for(auto it = rel.data.begin() ; it != rel.data.end() ; ) {
            int n = 0;
            auto itx = it;
            for(; itx != rel.data.end() && *it == *itx ; ++itx, n++);
            index_t p = *it;
            if (p >= density.size())
                std::fill_n(std::back_inserter(density), p + 1 - density.size(), 0);
            density[p]++;
            it = itx;
        }
    }
    for(index_t p = 0 ; p < density.size() ; p++)
        std::cout << p << " " << density[p] << "\n";

    return 0;
}
