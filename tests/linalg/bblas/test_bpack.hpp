#ifndef TEST_BPACK_HPP_
#define TEST_BPACK_HPP_

#include <stdio.h>              // for printf
#include <set>                  // for set
#include <string>               // for string, operator<
#include <vector>               // for vector
#include "test_bblas_base.hpp"

struct test_bpack : public test_bblas_base
{
    test_bpack(unsigned int n) : test_bblas_base(n) {}

    /* in test_bpack.cpp */
    /* don't call this "bpack", as this would clash with the bpack
     * template ! */
    static tags_t do_bpack_tags;
    void do_bpack();

    template<typename T>
    void meta_bpack();

    /* internal tests */
    template<typename T>
    int test_invert_triangular(unsigned int m, unsigned int n);

    void banner()
    {
        printf("-- bpack tests (blocks of matrices) --\n");
    }
    void operator()(std::vector<std::string> const & tests, std::set<std::string> & seen)
    {
        int has_banner = 0;
        for(auto const & s : tests) {
            bool match = false;
            if (matches(s, do_bpack_tags, match)) {
                if (!has_banner++) banner();
                do_bpack();
            }
            if (match) seen.insert(s);
        }
    }
};

#endif	/* TEST_BPACK_HPP_ */

