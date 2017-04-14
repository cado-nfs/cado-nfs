#include "cado.h"
#include "minimum_spanning_tree.hpp"
#include <iostream>
#include <climits>
#include "macros.h"

using namespace std;

int main()
{
    vector<int> A;
    int _ = INT_MAX;
    // this is the graph from
    // https://en.wikipedia.org/wiki/File:Minimum_spanning_tree.svg
    for(auto x : {
        _,  4,  1,  4,  _,  _,  _,  _,  _,  _,
        4,  _,  5,  _,  9,  9,  _,  7,  _,  _,
        1,  5,  _,  3,  _,  _,  _,  9,  _,  _,
        4,  _,  3,  _,  _,  _,  _, 10,  _, 18,
        _,  9,  _,  _,  _,  2,  4,  _,  6,  _,
        _,  9,  _,  _,  2,  _,  2,  8,  _,  _,
        _,  _,  _,  _,  4,  2,  _,  9,  3,  9,
        _,  7,  9, 10,  _,  8,  9,  _,  _,  8,
        _,  _,  _,  _,  6,  _,  3,  _,  _,  9,
        _,  _,  _, 18,  _,  _,  9,  8,  9,  _,
        }) A.push_back(x);
    auto mst = minimalSpanningTreePrimNaive(A, 10);
    cout << "MST has weight " << mst.first << "\n";
    for(auto const & x : mst.second)
        cout << " " << x.first << "--" << x.second << "\n";
    ASSERT_ALWAYS(mst.first == 38);
    return 0;
}

/*
 * we expect an MST of weight 38, which the code would be welcome to
 * print nicely as follows
 0--2--3
 |
 \--1--7--9
       |
       \--5--6--8
          |
          \--4
 */

