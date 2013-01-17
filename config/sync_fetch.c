#include <stdlib.h>
int main () {
  int o, *p;
  __sync_fetch_and_add(&o, 1);
  /* according to http://trac.wxwidgets.org/ticket/4542, one should use the
     return value for a complete configure test */
  p = __sync_sub_and_fetch(&o, 1);
  /* make sure we use the return value */
  printf ("%lx", p);
  exit (0);
}
