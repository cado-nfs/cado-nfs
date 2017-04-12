#include "cado.h"
#include "get_successive_minima.hpp"
#include <iostream>
#include <iterator>

int main()
{
    std::vector<int> foo;
    for(int i = 0 ; i < 200 ; i++) foo.push_back(rand() % 97);
    std::vector<size_t> mins = get_successive_minima(foo.begin(), foo.end(), 10);
    std::vector<size_t> maxs = get_successive_minima(foo.begin(), foo.end(), 10, std::greater<int>());

    std::cout << "10 minima:";
    for(size_t i = 0 ; i < mins.size() ; i++)
        std::cout << " " << mins[i] << ":" << foo[mins[i]];
    std::cout << std::endl;

    std::cout << "10 maxima: ";
    for(size_t i = 0 ; i < maxs.size() ; i++)
        std::cout << " " << maxs[i] << ":" << foo[maxs[i]];
    std::cout << std::endl;
}

