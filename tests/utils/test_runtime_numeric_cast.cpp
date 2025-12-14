#include "cado.h" // IWYU pragma: keep

#include <limits>
#include <stdexcept>
#include <iostream>

#include <cstdlib>
#include <cstdint>

#include "runtime_numeric_cast.hpp"

// NOLINTBEGIN(bugprone-empty-catch)
template<typename FROM_TYPE, typename TO_TYPE>
static bool test_pass_narrowing_min()
{
    FROM_TYPE v_from = std::numeric_limits<TO_TYPE>::min();
    TO_TYPE v_to = runtime_numeric_cast<TO_TYPE>(v_from);
    return v_to && v_from;
}

template<typename FROM_TYPE, typename TO_TYPE>
static bool test_pass_narrowing_max()
{
    FROM_TYPE v_from = std::numeric_limits<TO_TYPE>::max();
    TO_TYPE v_to = runtime_numeric_cast<TO_TYPE>(v_from);
    return v_to && v_from;
}

template<typename FROM_TYPE, typename TO_TYPE>
static bool test_fail_narrowing_min()
{
    FROM_TYPE v_from = std::numeric_limits<TO_TYPE>::min();
    TO_TYPE v_to = 0;
    try {
        v_from--;
        v_to = runtime_numeric_cast<TO_TYPE>(v_from);
        throw std::runtime_error("we should not reach here");
    } catch (runtime_numeric_cast_details::underflow const & e) {
    }
    return v_to && v_from;
}

template<typename FROM_TYPE, typename TO_TYPE>
static bool test_fail_narrowing_max()
{
    FROM_TYPE v_from = std::numeric_limits<TO_TYPE>::max();
    TO_TYPE v_to = 0;
    try {
        v_from++;
        v_to = runtime_numeric_cast<TO_TYPE>(v_from);
        throw std::runtime_error("we should not reach here");
    } catch (runtime_numeric_cast_details::overflow const & e) {
    }
    return v_to && v_from;
}

template<typename FROM_TYPE, typename TO_TYPE>
static bool test_pass_strict_widening_max()
{
    FROM_TYPE v_from = std::numeric_limits<FROM_TYPE>::max();
    TO_TYPE v_to = runtime_numeric_cast<TO_TYPE>(v_from);
    return v_to && v_from;
}

template<typename FROM_TYPE, typename TO_TYPE>
static bool test_pass_widening_min()
{
    FROM_TYPE v_from = std::numeric_limits<FROM_TYPE>::min();
    TO_TYPE v_to = runtime_numeric_cast<TO_TYPE>(v_from);
    return v_to && v_from;
}

template<typename FROM_TYPE, typename TO_TYPE>
static bool test_pass_widening_max()
{
    FROM_TYPE v_from = std::numeric_limits<FROM_TYPE>::max();
    TO_TYPE v_to = runtime_numeric_cast<TO_TYPE>(v_from);
    return v_to && v_from;
}

template<typename FROM_TYPE, typename TO_TYPE>
static bool test_fail_widening_min()
{
    FROM_TYPE v_from = std::numeric_limits<FROM_TYPE>::min();
    TO_TYPE v_to = 0;
    try {
        v_to = runtime_numeric_cast<TO_TYPE>(v_from);
        throw std::runtime_error("we should not reach here");
    } catch (runtime_numeric_cast_details::underflow const & e) {
    }
    return v_to && v_from;
}
// NOLINTEND(bugprone-empty-catch)


int main()
{
    try {
        int r = 0;

        /* both signed, narrowing, pass */
        r += test_pass_narrowing_min<int64_t, int32_t>();
        r += test_pass_narrowing_max<int64_t, int32_t>();
        /* both signed, narrowing, underflow */
        r += test_fail_narrowing_min<int64_t, int32_t>();
        /* both signed, narrowing, overflow */
        r += test_fail_narrowing_max<int64_t, int32_t>();

        /* both signed, same size, pass */
        r += test_pass_widening_min<int32_t, int32_t>();
        r += test_pass_widening_max<int32_t, int32_t>();
        r += test_pass_widening_min<int64_t, int64_t>();
        r += test_pass_widening_max<int64_t, int64_t>();

        /* both signed, strict widening, pass */
        r += test_pass_widening_min<int32_t, int64_t>();
        r += test_pass_widening_max<int32_t, int64_t>();

        /* both unsigned, narrowing, pass */
        r += test_pass_narrowing_min<uint64_t, uint32_t>();
        r += test_pass_narrowing_max<uint64_t, uint32_t>();
        /* both unsigned, narrowing, overflow */
        r += test_fail_narrowing_max<uint64_t, uint32_t>();
        /* both unsigned, same size, pass */
        r += test_pass_widening_min<uint32_t, uint32_t>();
        r += test_pass_widening_max<uint32_t, uint32_t>();
        r += test_pass_widening_min<uint64_t, uint64_t>();
        r += test_pass_widening_max<uint64_t, uint64_t>();
        /* both unsigned, strict widening, pass */
        r += test_pass_widening_min<uint32_t, uint64_t>();
        r += test_pass_widening_max<uint32_t, uint64_t>();

        /* to signed, narrowing, pass */
        r += test_pass_narrowing_max<uint64_t, int32_t>();
        /* to signed, narrowing, overflow */
        r += test_fail_narrowing_max<uint64_t, int32_t>();
        /* to signed, same size, pass */
        r += test_pass_narrowing_max<uint64_t, int64_t>();
        /* to signed, same size, overflow */
        r += test_fail_narrowing_max<uint64_t, int64_t>();
        /* to signed, strict widening, pass */
        r += test_pass_strict_widening_max<uint32_t, int64_t>();

        /* to unsigned, narrowing, pass */
        r += test_pass_narrowing_max<int64_t, uint32_t>();
        /* to unsigned, narrowing, underflow */
        r += test_fail_narrowing_min<int64_t, uint32_t>();
        /* to unsigned, same size, pass */
        r += test_pass_widening_max<int64_t, uint64_t>();
        r += test_pass_widening_max<int32_t, uint32_t>();
        /* to unsigned, same size, underflow */
        r += test_fail_widening_min<int64_t, uint64_t>();
        r += test_fail_widening_min<int32_t, uint32_t>();
        /* to unsigned, strict widening, pass */
        r += test_pass_widening_max<int32_t, uint64_t>();
        /* to unsigned, strict widening, underflow */
        r += test_fail_widening_min<int32_t, uint64_t>();


        return r ? EXIT_SUCCESS : EXIT_FAILURE;
    } catch (runtime_numeric_cast_details::failure const & e) {
        std::cerr << e.what() << "\n";
        return EXIT_FAILURE;
    } catch (std::runtime_error const & e) {
        std::cerr << e.what() << "\n";
        return EXIT_FAILURE;
    }
}
