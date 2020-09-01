#ifndef TEST_BBLAS_LEVEL4_HPP_
#define TEST_BBLAS_LEVEL4_HPP_

#include <set>
#include <string>
#include <vector>

#include "test_bblas_base.hpp"

struct test_bblas_level4 : public test_bblas_base
{
    test_bblas_level4(unsigned int n) : test_bblas_base(n) {}

    /* in test_bblas_level4_pluq.cpp */
    static tags_t pluq_tags;
    void pluq();

    /* in test_bblas_level4_gauss.cpp */
    static tags_t gauss_tags;
    void gauss();

    /* in test_bblas_level4_ple.cpp */
    static tags_t ple_tags;
    void ple();

    template<typename T>
    void meta_ple();
    /* unit tests for the internal ops of PLE */
    template<typename T>
    int test_PLE_find_pivot(unsigned int m, unsigned int n);
    template<typename T>
    int test_PLE_propagate_pivot(unsigned int m, unsigned int n);
    template<typename T>
    int test_PLE_propagate_row_permutations(unsigned int m, unsigned int n);
    template<typename T>
    int test_PLE_move_L_fragments(unsigned int m, unsigned int n);
    template<typename T>
    int test_PLE(unsigned int m, unsigned int n);

#ifdef  HAVE_M4RI
    void m4ri_plu_tests(int n);
#endif

    void banner()
    {
        printf("-- level-4 tests (reductions / factorizations of matrices) --\n");
    }
    void operator()(std::vector<std::string> const & tests, std::set<std::string> & seen)
    {
        int has_banner = 0;
        for(auto const & s : tests) {
            bool match = false;
            if (matches(s, pluq_tags, match)) {
                if (!has_banner++) banner();
                pluq();
            }
            if (matches(s, gauss_tags, match)) {
                if (!has_banner++) banner();
                gauss();
            }
            if (matches(s, ple_tags, match)) {
                if (!has_banner++) banner();
                ple();
            }
            if (match) seen.insert(s);
        }
    }
};

#endif	/* TEST_BBLAS_LEVEL4_HPP_ */
