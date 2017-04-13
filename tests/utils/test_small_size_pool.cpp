#include "cado.h"
#include <vector>
#include <iostream>
#define DEBUG_SMALL_SIZE_POOL
#include "small_size_pool.hpp"

using namespace std;

int main()
{
    small_size_pool<int, unsigned int, 4> pool;

    vector<pair<unsigned int, unsigned int>> values;

    for(int i = 0 ; i < 100 ; i++) {
        /* pick a new size, add a new something */
        unsigned int w = rand() % 40 + 1;
        unsigned int v = pool.alloc(w);
        values.push_back(make_pair(w,v));
        int * blah = pool(w, v);
        for(unsigned int k = 0 ; k < w ; k++) {
            blah[k] = rand();
        }
    }
    /* free one half of them */
    for(int i = 0 ; i < 50 ; i++) {
        if (rand() & 4) continue;
        pool.free(values[i]);
    }
    /* now reallocate everyone */
    for(int i = 0 ; i < 100 ; i++) {
        unsigned int w = rand() % 40 + 1;
        pool.realloc(values[i], w);
        int * blah = pool(values[i]);
        for(unsigned int k = 0 ; k < w ; k++) {
            blah[k] = rand();
        }
    }
    cout << pool;
}

