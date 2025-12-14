#include "cado.h" // IWYU pragma: keep

/**
 * @file MinMaxHeap.hpp
 * @author John Sullivan (jsull003 at ucr.edu)
 * @date May 16, 2012
 *
 * @brief tests for @ref MinMaxHeap.
 **/

#include <cstdlib>

#include <vector>
#include <algorithm>

#include <gmp.h>

#include "macros.h"
#include "min_max_heap.hpp"
#include "tests_common.h"

struct min_max_heap_test_case {
    std::vector<int> unsortedValues;

    /// The values in unsorted values in ascending order.
    std::vector<int> sortedValues;

    explicit min_max_heap_test_case(size_t nvalues)
    {
        unsortedValues.resize(nvalues);
        sortedValues.resize(nvalues);

        for (unsigned int i = 0; i < nvalues; ++i)
            unsortedValues[i] = sortedValues[i] = (int) gmp_urandomb_ui(state, ULONG_BITS);

        std::ranges::sort(sortedValues);

        for (unsigned int i = 0; i < nvalues; ++i)
            heap.push(unsortedValues[i]);
    }

    min_max_heap<int> heap;
};

static void test_SatisfiesMaxHeap(size_t nvalues)
{
    min_max_heap_test_case T(nvalues);
    // Get all the values out of the heap and store them in ascending order
    std::vector<int> values(nvalues);
    for (int i = (int)nvalues - 1; i >= 0; --i)
        values[i] = T.heap.popMax();

    FATAL_ERROR_CHECK(
            !std::ranges::equal(T.sortedValues, values), "Values were not popped in the correct order.");
}

static void test_SatisfiesMinHeap(size_t nvalues)
{
    min_max_heap_test_case T(nvalues);
    // Get all the values out of the heap and store them in ascending order
    std::vector<int> values(nvalues);
    for (unsigned int i = 0; i < nvalues; ++i)
        values[i] = T.heap.popMin();

    FATAL_ERROR_CHECK(
            !std::ranges::equal(T.sortedValues, values), "Values were not popped in the correct order.");
}

static void test_SatisfiesMinMaxHeap_HalfAndHalf(size_t nvalues)
{
    min_max_heap_test_case T(nvalues);
    // Get all the values out of the heap and store them in ascending order
    std::vector<int> values(nvalues);
    for (unsigned int i = 0; i < nvalues / 2; ++i)
        values[i] = T.heap.popMin();

    for (int i = (int)nvalues - 1; i >= (int)nvalues / 2; --i)
        values[i] = T.heap.popMax();

    FATAL_ERROR_CHECK(
            !std::ranges::equal(T.sortedValues, values),
            "Values were not popped in the correct order.");
}

static void test_SatisfiesMinMaxHeap_Alternating(size_t nvalues)
{
    min_max_heap_test_case T(nvalues);
    // Get all the values out of the heap and store them in ascending order
    std::vector<int> values(nvalues);
    for (unsigned int i = 0; i < nvalues / 2; ++i)
    {
        values[i] = T.heap.popMin();
        values[nvalues - i - 1] = T.heap.popMax();
    }

    if (nvalues % 2 == 1)
        values[nvalues / 2 + 1] = T.heap.popMin();

    /*
    for (unsigned int i = 0; i < nvalues; ++i)
        std::cout << values[i] << " ?= " << T.sortedValues[i] << "\n";
        */

    FATAL_ERROR_CHECK(
            !std::ranges::equal(T.sortedValues, values),
            "Values were not popped in the correct order.");
}

int main(int argc, char const * argv[])
{
    unsigned long iter = 20000;
    tests_common_cmdline (&argc, &argv, PARSE_SEED | PARSE_ITER);
    tests_common_get_iter (&iter);
    test_SatisfiesMaxHeap(iter);
    test_SatisfiesMinHeap(iter);
    test_SatisfiesMinMaxHeap_HalfAndHalf(iter);
    test_SatisfiesMinMaxHeap_Alternating(iter);
    tests_common_clear ();
    exit (EXIT_SUCCESS);
}
