#include "cado.h" // IWYU pragma: keep
#include "smallset.hpp" // IWYU pragma: keep
#include "macros.h"

#if defined(HAVE_SSE2) && GNUC_VERSION_ATLEAST(4,7,0)
/* see utils/smallset.h */

template <int SIZE, typename ELEMENTTYPE>
struct
test_smallset
{
    static_assert(SIZE != 0, "use size-0 specialization instead");
    static constexpr const size_t nr_items = smallset<SIZE, ELEMENTTYPE>::nr_items;
    ELEMENTTYPE items[nr_items];

    test_smallset() {
        for (size_t i = 0; i < nr_items; i++)
            items[i] = i;
    }

    void operator()() {
        for (size_t i = 1; i <= nr_items; i++) {
            smallset<SIZE, ELEMENTTYPE> set(items, i);
            ASSERT_ALWAYS(set.contains(0));
            ASSERT_ALWAYS(set.contains(i-1));
            ASSERT_ALWAYS(!set.contains(i));
            ASSERT_ALWAYS(set.contains(2) == (i > 2));
        }
    }
};

template <typename ELEMENTTYPE>
struct
test_smallset<0, ELEMENTTYPE>
{
    void operator()() {
        smallset<0, ELEMENTTYPE> set(NULL, 0);
        ASSERT_ALWAYS(!set.contains(0));
    }
};

template <typename ELEMENTTYPE>
void
test_one_type()
{
  test_smallset<0, ELEMENTTYPE>()();
  test_smallset<1, ELEMENTTYPE>()();
  test_smallset<2, ELEMENTTYPE>()();
  test_smallset<3, ELEMENTTYPE>()();
  test_smallset<4, ELEMENTTYPE>()();
  test_smallset<5, ELEMENTTYPE>()();
}

// coverity[root_function]
int main()
{
  test_one_type<unsigned char>();
  test_one_type<unsigned short>();
  test_one_type<unsigned int>();
  return(0);
}

#else

int main()
{
  return 0;
}

#endif
