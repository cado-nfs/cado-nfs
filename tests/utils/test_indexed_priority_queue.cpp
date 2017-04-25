#include "cado.h"
#include <iostream>
#define DEBUG_INDEXED_PRIORITY_QUEUE
#include "indexed_priority_queue.hpp"
#include "macros.h"
#include <map>
#include <string>

int test1() {
    std::vector<int> foo(20);
    for(int i = 0, j = 42 ; i < 20 ; i++, j = (j * 17) % 97) {
        foo[i] = j;
    }

    indexed_priority_queue<int, int, std::greater<std::pair<int,int>>> q(foo.begin(), foo.end());

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

    std::cout << q;
    
    for( ; !q.empty() ; q.pop()) {
        auto v = q.top();
        std::cout << "col " << v.first << " with value " << v.second << std::endl;
    }
    return 0;
}

int test2() {
    std::vector<std::pair<int, int>> foo;
    // indexed_priority_queue<int, int, std::greater<int>, int, std::map<int, int>> q;
    sparse_indexed_priority_queue<int, int, std::greater<std::pair<int,int>>> q;

    for(int i = 0, j = 42 ; i < 20 ; i++, j = (j * 17) % 97) {
        foo.push_back(std::make_pair(i * 17 + rand() % 17, j));
        q.push(foo.back());
    }

    std::cout << q;
    // std::cout << q.is_heap() << std::endl;
    ASSERT_ALWAYS(q.is_heap());

    for(size_t k = 0 ; k < foo.size() ; k++) {
        int j = foo[rand() % foo.size()].first;
        int v = rand() % 97;
        if (v < 25) {
            std::cout << "== remove col " << j << std::endl;
            q.remove(j);
        } else {
            std::cout << "== upate value for col " << j << " to " << v << std::endl;
            q.update(j,v);
        }
        // std::cout << q;
    }

    std::cout << q;
    
    for( ; !q.empty() ; q.pop()) {
        auto v = q.top();
        std::cout << "col " << v.first << " with value " << v.second << std::endl;
    }
    return 0;
}

int test3()
{
    high_score_table<std::string, int> HS(5);
    HS.push("mkz", 1000000);
    HS.push("jmz", 123700);
    HS.push("eth", 1234500);
    HS.push("jml", 238500);
    HS.push("mwz", 691000);
    HS.push("chf", 224000);
    HS.push("dav", 847000);
    HS.remove("jml");
    HS.remove("mkz");
    HS.push("mrs", 610000);
    HS.push("sam", 780000);
    HS.set_depth(3);
    HS.push("xyz", 810000);
    HS.push("waw", 123400);
    HS.push("zaz", 582900);
    HS.set_depth(5);
    HS.push("bob", 1225000);
    HS.push("alf", 795200);

    std::cout << HS;
    for( ; !HS.empty() ; HS.pop()) {
        auto v = HS.top();
        std::cout << v.first << " " << v.second << std::endl;
    }
    return 0;
}

int main()
{
    test1();
    test2();
    test3();
}
