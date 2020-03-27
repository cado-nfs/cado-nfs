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

    bool operator()(std::string const & s)
    {
            bool match = false;
            if (matches(s, level2a_tags, match)) level2a();
            if (matches(s, level2_tags, match)) level2();
            return match;
    }
};

#endif	/* TEST_BBLAS_LEVEL2_HPP_ */
