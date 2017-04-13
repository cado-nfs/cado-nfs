#include "cado.h"
#include <iostream>
#define DEBUG_INDEXED_PRIORITY_QUEUE
#include "indexed_priority_queue.hpp"
#include "macros.h"

int main() {
    std::vector<int> foo(20);
    for(int i = 0, j = 42 ; i < 20 ; i++, j = (j * 17) % 97) {
        foo[i] = j;
    }

    indexed_priority_queue<int, int, std::greater<int> > q(foo.begin(), foo.end());

    // std::cout << q;
    // std::cout << q.is_heap() << std::endl;
    ASSERT_ALWAYS(q.is_heap());

    for(int i = 0 ; i < 20 ; i++) {
        int j = rand() % q.size();
        int v = rand() % 97;
        std::cout << "== upate value for col " << j << " to " << v << std::endl;
        q.update(j,v);
        // std::cout << q;
    }
    
    for( ; !q.empty() ; q.pop()) {
        auto v = q.top();
        std::cout << "col " << v.first << " with value " << v.second << std::endl;
    }
}
