#ifndef TEST_BBLAS_LEVEL5_HPP_
#define TEST_BBLAS_LEVEL5_HPP_

#include <stdio.h>              // for printf
#include <set>                  // for set
#include <string>               // for string, operator<
#include <vector>               // for vector
#include "test_bblas_base.hpp"

struct test_bblas_level5 : public test_bblas_base
{
    test_bblas_level5(unsigned int n) : test_bblas_base(n) {}

    static tags_t polmul_tags;
    void polmul();

    static tags_t polblockmul_tags;
    void polblockmul();

    static tags_t matpolmul_tags;
    void matpolmul();

    static tags_t matpolscale_tags;
    void matpolscale();

    void banner()
    {
        printf("-- level-5 tests (polynomials) --\n");
    }
    void operator()(std::vector<std::string> const & tests, std::set<std::string> & seen)
    {
        int has_banner = 0;
        for(auto const & s : tests) {
            bool match = false;
            if (matches(s, polmul_tags, match)) {
                if (!has_banner++) banner();
                polmul();
            }
            if (matches(s, polblockmul_tags, match)) {
                if (!has_banner++) banner();
                polblockmul();
            }
            if (matches(s, matpolmul_tags, match)) {
                if (!has_banner++) banner();
                matpolmul();
            }
            if (matches(s, matpolscale_tags, match)) {
                if (!has_banner++) banner();
                matpolscale();
            }
            if (match) seen.insert(s);
        }
    }
};

#endif /* TEST_BBLAS_LEVEL5_HPP_ */
