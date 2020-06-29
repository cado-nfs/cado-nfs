#include <stdio.h>

static inline void
ularith_div_2ul_ul_ul_r (unsigned long *r, unsigned long a1,
                 const unsigned long a2, const unsigned long b)
{
  __asm__(
            "# ularith_div_2ul_ul_ul_r: divq %0 %1 %2 %3\n\t"
            "divq %3"
            : "+a" (a1), "=d" (*r)
            : "1" (a2), "rm" (b)
            : "cc");
}

int
modredcul_initmod_ul (unsigned long *m, const unsigned long s)
{
  if (s == 1UL) {
    *m = 0UL;
  }
  else
    {
      *m = 1UL;
      ularith_div_2ul_ul_ul_r (m, 0UL, 1UL, s);
    }

  // fprintf(stderr, "m = %lu\n", *m);

  if (sizeof(unsigned long) == 8)
    if (*m != 8241071186)
      return 1;
  else
    return 0;
}

int main(void)
{
    unsigned long m;
    unsigned long pp = 10003800361UL;
    int re = modredcul_initmod_ul (&m,  pp);
    return re;
}
