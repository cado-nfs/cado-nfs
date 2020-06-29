#include "cado.h" // IWYU pragma: keep
#include <stdlib.h>
#include "getprime.h"
#include "macros.h"

void
test_getprime ()
{
  unsigned long p;
  int i;
  prime_info pi;

  /* there is no really something random we can test here, since
     getprime can start only at 2 */
  prime_info_init (pi);
  for (p = 2, i = 0; i < 1000000; p = getprime_mt (pi), i++);
  ASSERT_ALWAYS (p == 15485867);
  prime_info_clear (pi);

  prime_info_init (pi);
  for (p = 2, i = 0; i < 1000000; p = getprime_mt (pi), i++);
  ASSERT_ALWAYS (p == 15485867);
  prime_info_clear (pi);
}

int
main ()
{
  test_getprime ();
  exit (EXIT_SUCCESS);
}
