
static void 
mul(unsigned long *rlo, unsigned long *rhi, const unsigned long a, const unsigned long b)
{
  __asm__ volatile(
   "    umull   %[rlo], %[rhi], %[a], %[b]\n\t"
  : [rhi] "=r" (*rhi), [rlo] "=r" (*rlo)
  : [a] "r" (a), [b] "r" (b)
  );
}

int main()
{
    unsigned long z[2] = {1, 2};

    mul(&z[0], &z[1], 17, 42);
    
    return (z[0] == 17*42 && z[1] == 0) ? 0 : 1;
}
