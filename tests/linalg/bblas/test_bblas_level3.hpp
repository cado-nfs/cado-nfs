#ifndef TEST_BBLAS_LEVEL3_HPP_
#define TEST_BBLAS_LEVEL3_HPP_

#include "test_bblas_base.hpp"

struct test_bblas_level3 : public test_bblas_base
{
    test_bblas_level3(unsigned int n) : test_bblas_base(n) {}

    static tags_t level3a_tags;
    void level3a();

    static tags_t transpose_tags;
    void transpose();

    static tags_t matpoly_polmat_tags;
    void matpoly_polmat();

    static tags_t matmul_tags;
    void matmul();

    static tags_t rank_n_update_tags;
    void rank_n_update();

    static tags_t level3c_tags;
    void level3c_list();
    void level3c();

    static tags_t trsm_tags;
    void trsm();

    void operator()(std::vector<std::string> const & tests, std::set<std::string> & seen)
    {
        printf("-- level-3 tests (operations on matrices) --\n");
        for(auto const & s : tests) {
            bool match = false;
            if (matches(s, level3a_tags, match)) level3a();
            if (matches(s, transpose_tags, match)) transpose();
            if (matches(s, matpoly_polmat_tags, match)) matpoly_polmat();
            if (matches(s, matmul_tags, match)) matmul();
            if (matches(s, rank_n_update_tags, match)) rank_n_update();
            if (matches(s, level3c_tags, match)) level3c();
            if (matches(s, trsm_tags, match)) trsm();
            if (match) seen.insert(s);
        }
    }
};

#endif	/* TEST_BBLAS_LEVEL3_HPP_ */
