#include "cado.h" // IWYU pragma: keep

#include <cstdio>
/* This file implements some methods that work more or less the same
   with Modulus64 and ModulusREDC64. E.g., div3() and gcd() work
   unchanged with plain and Montgomery representation (so we can work on
   the stored residue directly, whatever its representation is);
   jacobi() converts to plain Integer first, the others use
   only Modulus methods.
   Speed-critical functions need to be rewritten in assembly for REDC,
   but this is a start.
*/

// scan-headers: stop here


#if 1


#else

/*
#!/usr/bin/env python3
# Python program to create mult.h

def rate(k,b):
  # The 0.25 magic constant here tries to estimate the ratio m/x, 
  # to minimize (x+c*m)/2^b
  r = (abs(k)*0.25 + 1.)/2**b
  # print ("rate(" + str(k) + ", " + str(b) + ") = " + str(r))
  return(r)

def bestb(k):
  best_b = 0
  best_r = 1
  best_c = 0
  for b in range(1, 8):
    c = k % (2**b)
    r = rate(c, b) 
    if r < best_r: 
      best_r=r 
      best_b = b 
      best_c = c
    c = - (2**b - c)
    r = rate(c, b)
    if r < best_r: 
      best_r=r 
      best_b = b 
      best_c = c
  # print ("bestb(" + str(k) + ") = " + str(best_b))
  return([k, best_b, best_c, 1/best_r])


r = [str(-bestb(2*i+1)[2]) for i in range(0, 128) ]
print("static char mult[128] = {" + ", ".join(r) + "};")
*/

#include "mult.h"
#include "macros.h"
static unsigned char invmod8[256] = {
0, 1, 0, 171, 0, 205, 0, 183, 0, 57, 0, 163, 0, 197, 0, 239, 0, 241, 0, 27, 0, 61, 0, 167, 0, 41, 0, 19, 0, 53, 0, 223, 0, 225, 0, 139, 0, 173, 0, 151, 0, 25, 0, 131, 0, 165, 0, 207, 0, 209, 0, 251, 0, 29, 0, 135, 0, 9, 0, 243, 0, 21, 0, 191, 0, 193, 0, 107, 0, 141, 0, 119, 0, 249, 0, 99, 0, 133, 0, 175, 0, 177, 0, 219, 0, 253, 0, 103, 0, 233, 0, 211, 0, 245, 0, 159, 0, 161, 0, 75, 0, 109, 0, 87, 0, 217, 0, 67, 0, 101, 0, 143, 0, 145, 0, 187, 0, 221, 0, 71, 0, 201, 0, 179, 0, 213, 0, 127, 0, 129, 0, 43, 0, 77, 0, 55, 0, 185, 0, 35, 0, 69, 0, 111, 0, 113, 0, 155, 0, 189, 0, 39, 0, 169, 0, 147, 0, 181, 0, 95, 0, 97, 0, 11, 0, 45, 0, 23, 0, 153, 0, 3, 0, 37, 0, 79, 0, 81, 0, 123, 0, 157, 0, 7, 0, 137, 0, 115, 0, 149, 0, 63, 0, 65, 0, 235, 0, 13, 0, 247, 0, 121, 0, 227, 0, 5, 0, 47, 0, 49, 0, 91, 0, 125, 0, 231, 0, 105, 0, 83, 0, 117, 0, 31, 0, 33, 0, 203, 0, 237, 0, 215, 0, 89, 0, 195, 0, 229, 0, 15, 0, 17, 0, 59, 0, 93, 0, 199, 0, 73, 0, 51, 0, 85, 0, 255
};

static inline int
s_val(unsigned int s) {
  return ((s & 2) == 0) ? 1 : -1;
}

static int
Modulus::jacobi1 (const Residue &a_par)
{
  uint64_t x, m;
  unsigned int s, j;

  x = mod_get_ul (a_par, m_par);
  m = mod_getmod_ul (m_par);
  ASSERT (x < m);
  ASSERT(m % 2 == 1);

  j = u64arith_ctz(x);
  x = x >> j;

  s = ((j<<1) & (m ^ (m>>1)));

  while (x > 1) {
    uint64_t t;
    unsigned char inv;

    // printf ("kronecker(%lu, %lu) == %d * kronecker(%lu, %lu)\n",
    //        mod_get_ul(a_par, m_par), mod_getmod_ul (m_par), s_val(s), x, m);
    /* Here x < m. Swap to make x > m */
    t = m;
    m = x; 
    x = t;
    s = s ^ (x&m);

    /* Now x > m */
    inv = invmod8[(unsigned char)m];
    do {
      /* Do a REDC-like step. We want a multiplier k such that the low 
         8 bits of x+k*m are all zero. 
         That is, we want k = -x/m (mod 2^8). */
      unsigned char k;
      uint64_t t1;
      long int c, t2;
      // const uint64_t old_x = x;
      
      k = inv * (unsigned char)x;
      // ASSERT_ALWAYS((k & 1) == 1);
      c = mult[k / 2];
      
      /* Compute x+cm */
      long tmp = c >> 63;
      u64arith_mul_1_1_2 (&t1, (uint64_t *)&t2, (c ^ tmp) - tmp, m);
      t2 ^= tmp;
      t1 ^= tmp;
      u64arith_add_1_2 (&t1, (uint64_t *)&t2, x-tmp);
      tmp = ((long) t2) >> 63;
      
      t2 ^= tmp;
      t1 ^= tmp;
      s ^= m & tmp;
      u64arith_add_1_2 (&t1, (uint64_t *)&t2, -tmp);
      // ASSERT_ALWAYS(t2 >= 0);

      if (t1 == 0) {
        if (t2 == 0) {
          x = 0;
          break;
        }
        t1 = t2;
        /* Divided by 2^64 which is square, so no adjustment to s */
        t2 = 0;
      }

      j = u64arith_ctz(t1);
      u64arith_shrd (&t1, t2, j);
      // ASSERT_ALWAYS((t2 >> j) == 0);
      x = t1;
      s ^= ((j<<1) & (m ^ (m>>1)));
      // ASSERT_ALWAYS(x < old_x);
      // printf ("%f\n", (double)old_x / (double)x);
    } while (x >= m);
  }

  if (x == 0)
    return 0;
  return s_val(s);
}

int
Modulus::jacobi (const Residue &a)
{
  uint64_t x, m;
  unsigned int s, j;

  x = mod_get_ul (a_par, m_par);
  m = mod_getmod_ul (m_par);
  ASSERT (x < m);
  ASSERT(m % 2 == 1);

  if ((LONG_MAX - x) / 50 < m)
    return mod_jacobi1 (a_par, m_par);

  j = u64arith_ctz(x);
  x = x >> j;

  s = ((j<<1) & (m ^ (m>>1)));

  while (x > 1) {
    uint64_t t;
    unsigned char inv;

    // printf ("kronecker(%lu, %lu) == %d * kronecker(%lu, %lu)\n",
    //        mod_get_ul(a_par, m_par), mod_getmod_ul (m_par), s_val(s), x, m);
    /* Here x < m. Swap to make x > m */
    t = m;
    m = x; 
    x = t;
    s = s ^ (x&m);

    /* Now x > m */
    inv = invmod8[(unsigned char)m];
    do {
      /* Do a REDC-like step. We want a multiplier k such that the low 
         8 bits of x+k*m are all zero. 
         That is, we want k = -x/m (mod 2^8). */
      unsigned char k;
      long int c;
      // const uint64_t old_x = x;
      
      k = inv * x;
      // ASSERT_ALWAYS((k & 1) == 1);
      c = mult[k / 2];
      
      c = x + c*m;
      x = c;
      c >>= 63;
      x = (x ^ c) - c;
      
      if (x == 0) {
        break;
      }
      s ^= m & c;

      j = u64arith_ctz(x);
      x >>= j;
      s ^= ((j<<1) & (m ^ (m>>1)));
      // printf ("%f\n", (double)old_x / (double)x);
    } while (x >= m);
  }

  if (x == 0)
    return 0;
  return s_val(s);
}

#endif
