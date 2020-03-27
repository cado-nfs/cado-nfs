#ifndef TEST_BBLAS_LEVEL2_HPP_
#define TEST_BBLAS_LEVEL2_HPP_

#include "test_bblas_base.hpp"

struct test_bblas_level2 : public test_bblas_base
{
    test_bblas_level2(unsigned int n) : test_bblas_base(n) {}

    static tags_t level2a_tags;
    void level2a();

    static tags_t level2_tags;
    void level2();

    void operator()(std::vector<std::string> const & tests, std::set<std::string> & seen)
    {
        printf("-- level-2 tests (inputs + outputs = 2 vectors + 1 matrix) --\n");
        for(auto const & s : tests) {
            bool match = false;
            if (matches(s, level2a_tags, match)) level2a();
            if (matches(s, level2_tags, match)) level2();
            if (match) seen.insert(s);
        }
    }
};

#endif	/* TEST_BBLAS_LEVEL2_HPP_ */
