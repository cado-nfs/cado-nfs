#include <stdlib.h>
/* 
  This program is a test case for GCC bug 58805, see
  http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58805
  This bug affects the cofactorization routines.
  If the gcc version is unaffected, the compiled executable returns with
  exit status 0, otherwise it aborts.
*/

static inline void bar(unsigned long *r)
{
  unsigned long t;
  __asm__ (
    "movq $42, %[t]\n\t"
    "movq %[t], %[r]\n\t"
    : [t] "=&r" (t), [r] "=r" (*r)
  );
}

static void __attribute__((__noinline__))
foo(int n, unsigned long *x, unsigned long *y)
{
  if (n == 0)
    bar(x);
  else
    bar(y);
}

int na = 0, nb = 1;
unsigned long xa = 17, ya = 17;
unsigned long xb = 17, yb = 17;

int main (void)
{
  foo(na, &xa, &ya);
  if (xa != 42 || ya != 17)
    __builtin_abort ();
  foo(nb, &xb, &yb);
  if (xb != 17 || yb != 42)
    __builtin_abort ();
  return 0;
}
