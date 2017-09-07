#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
int main (int argc, char * argv[]) {
  uint64_t o = argc, y, z;
  __sync_fetch_and_add(&o, 1);
  /* according to http://trac.wxwidgets.org/ticket/4542, one should use the
     return value for a complete configure test */
  y = __sync_sub_and_fetch(&o, 1);
  z = __sync_add_and_fetch(&o, y);

  /* make sure we use the return value */
  printf ("%lx", z);
  exit (0);
}
