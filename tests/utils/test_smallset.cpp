#include "cado.h" // IWYU pragma: keep
#include "smallset.hpp" // IWYU pragma: keep
#include "macros.h"

#if defined(HAVE_SSE2) && GNUC_VERSION_ATLEAST(4,7,0)
/* see utils/smallset.h */

template <int SIZE, typename ELEMENTTYPE>
void
test_smallset()
{
  ASSERT_ALWAYS(SIZE != 0);

  const size_t nr_items = smallset<SIZE, ELEMENTTYPE>::nr_items;
  ELEMENTTYPE items[nr_items];
  for (size_t i = 0; i < nr_items; i++)
    items[i] = i;

  for (size_t i = 1; i <= nr_items; i++) {
    smallset<SIZE, ELEMENTTYPE> set(items, i);
    ASSERT(set.contains(0));
    ASSERT(set.contains(i-1));
    ASSERT(!set.contains(i));
    ASSERT(set.contains(2) == (i > 2));
  }
}

template <typename ELEMENTTYPE>
void
test_smallset<0, ELEMENTTYPE>()
{
    smallset<SIZE, ELEMENTTYPE> set(items, 0);
    ASSERT(!set.contains(0));
}

template <typename ELEMENTTYPE>
void
test_one_type()
{
  test_smallset<0, ELEMENTTYPE>();
  test_smallset<1, ELEMENTTYPE>();
  test_smallset<2, ELEMENTTYPE>();
  test_smallset<3, ELEMENTTYPE>();
  test_smallset<4, ELEMENTTYPE>();
  test_smallset<5, ELEMENTTYPE>();
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
