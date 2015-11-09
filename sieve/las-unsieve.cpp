#include "cado.h"
#include <stdlib.h>
#include <string.h>

#include "las-unsieve.h"
#include "ularith.h"
#include "las-norms.h"
#include "las-debug.h"
#include "gcd.h"
#include "memory.h"

static const int verify_gcd = 0; /* Enable slow but thorough test */

/* Set every stride-th byte, starting at index 0, to 255 in an array of
   stride unsieve_pattern_t's, and set all other bytes to 0. */
static void
minisieve(unsieve_pattern_t * const array, const size_t stride)
{
    memset (array, 0, stride * sizeof(unsieve_pattern_t));
    for (size_t i = 0; i < stride * sizeof(unsieve_pattern_t); i += stride)
      ((unsigned char *) array)[i] = 255;
}

unsieve_aux_data_srcptr
init_unsieve_data(uint32_t I)
{
  unsieve_aux_data_ptr us;
  us = (unsieve_aux_data_ptr) malloc_aligned(sizeof(unsieve_aux_data), 16);
  /* Store largest prime factor of k in us->lpf[k], 0 for k=0, 1 for k=1 */
  us->entries = (unsieve_entry_t *) malloc (sizeof (unsieve_entry_t) * I);
  FATAL_ERROR_CHECK(us->entries == NULL, "malloc failed");
  us->entries[0].lpf = 0U;
  us->entries[0].cof = 0U;
  us->entries[0].start = 0U;
  us->entries[1].lpf = 1U;
  us->entries[1].cof = 1U;
  us->entries[1].start = 0U;
  for (unsigned int k = 2U; k < I; k++)
    {
      unsigned int p, c = k;
      for (p = 2U; p * p <= c; p += 1U + p % 2U)
        {
          while (c % p == 0U)
            c /= p;
          if (c == 1U)
            break;
        }
      p = (c == 1U) ? p : c;
      c = k; do {c /= p;} while (c % p == 0);
      us->entries[k].lpf = p;
      us->entries[k].cof = c;
      us->entries[k].start = (I / 2U) % p;
    }

    /* Create pattern for sieving 3 */
    minisieve(us->pattern3, 3);
    /* Create pattern for sieving 5 */
    minisieve(us->pattern5, 5);
    /* Create pattern for sieving 7 */
    minisieve(us->pattern7, 7);

  return us;
}

void clear_unsieve_data(unsieve_aux_data_srcptr us)
{
  free (us->entries);
  free_aligned (us);
}

static inline void
unsieve_one_prime (unsigned char *line_start, const unsigned int p, 
                   const unsigned int j, const unsigned int start_idx,
                   const unsigned int I)
{
  unsigned int x, np = p; /* if 2|j, np=2p, else np=p */

  x = start_idx;
  if (j % 2U == 0U)
    {
      np += p;
      if (x % 2U == 0U)
        x += p;
    }
  for ( ; x < I; x += np)
    line_start[x] = 255;
}


static inline void
unsieve_3(unsigned char *line_start, const unsigned int start_idx,
          const unsigned int I, unsieve_aux_data_srcptr us)
{
  const unsigned int I_upt  = I / sizeof (unsieve_pattern_t);
  unsigned int i, pattern_idx;
  unsieve_pattern_t p0, p1, p2;
  unsieve_pattern_t * ul_line_start = (unsieve_pattern_t *) line_start;

  if (sizeof(unsieve_pattern_t) == 4) {
    /* -4^(-1) == 2 (mod 3) */
    pattern_idx = (2 * start_idx) % 3;
  } else if (sizeof(unsieve_pattern_t) == 8) {
    /* -8^(-1) == 1 (mod 3) */
    pattern_idx = start_idx;
  } else if (sizeof(unsieve_pattern_t) == 16) {
    /* -16^(-1) == 2 (mod 3) */
    pattern_idx = (2 * start_idx) % 3;
  } else
    abort();
  
  p0 = us->pattern3[pattern_idx];
  p1 = us->pattern3[(pattern_idx + 1) % 3];
  p2 = us->pattern3[(pattern_idx + 2) % 3];

  ASSERT_ALWAYS(((unsigned char *)&p0)[start_idx] == 255);
  
  /* Apply pattern to array */
  for (i = 0U; i < I_upt - 2U; i += 3U)
    {
      UNSIEVE_OR(ul_line_start[i], p0);
      UNSIEVE_OR(ul_line_start[i + 1], p1);
      UNSIEVE_OR(ul_line_start[i + 2], p2);
    }
  if (i < I_upt)
    UNSIEVE_OR(ul_line_start[i], p0);
  if (i + 1 < I_upt)
    UNSIEVE_OR(ul_line_start[i + 1], p1);
}


static inline void
unsieve_5(unsigned char *line_start, const unsigned int start_idx,
          const unsigned int I, unsieve_aux_data_srcptr us)
{
  const unsigned int I_upt  = I / sizeof (unsieve_pattern_t);
  unsigned int i;
  unsieve_pattern_t p0, p1, p2, p3, p4;
  unsieve_pattern_t * ul_line_start = (unsieve_pattern_t *) line_start;
  size_t pattern_idx;

  if (sizeof(unsieve_pattern_t) == 4) {
    /* -4^(-1) == 1 (mod 5) */
    pattern_idx = start_idx;
  } else if (sizeof(unsieve_pattern_t) == 8) {
    /* -8^(-1) == 3 (mod 5) */
    pattern_idx = (3 * start_idx) % 5;
  } else if (sizeof(unsieve_pattern_t) == 16) {
    /* -16^(-1) == -1 (mod 5) */
    pattern_idx = (5 - start_idx) % 5;
  } else
    abort();
  
  p0 = us->pattern5[pattern_idx];
  p1 = us->pattern5[(pattern_idx + 1) % 5];
  p2 = us->pattern5[(pattern_idx + 2) % 5];
  p3 = us->pattern5[(pattern_idx + 3) % 5];
  p4 = us->pattern5[(pattern_idx + 4) % 5];

  if (start_idx < sizeof(p0)) {
      ASSERT_ALWAYS(((unsigned char *)&p0)[start_idx] == 255);
  } else {
      ASSERT_ALWAYS(start_idx < 2 * sizeof(p0));
      ASSERT_ALWAYS(((unsigned char *)&p1)[start_idx - sizeof(p0)] == 255);
  }

  /* Apply pattern to array */
  for (i = 0U; i < I_upt - 4U; i += 5U)
    {
      UNSIEVE_OR(ul_line_start[i], p0);
      UNSIEVE_OR(ul_line_start[i + 1], p1);
      UNSIEVE_OR(ul_line_start[i + 2], p2);
      UNSIEVE_OR(ul_line_start[i + 3], p3);
      UNSIEVE_OR(ul_line_start[i + 4], p4);
    }
  if (i < I_upt)
    UNSIEVE_OR(ul_line_start[i], p0);
  if (i + 1 < I_upt)
    UNSIEVE_OR(ul_line_start[i + 1], p1);
  if (i + 2 < I_upt)
    UNSIEVE_OR(ul_line_start[i + 2], p2);
  if (i + 3 < I_upt)
    UNSIEVE_OR(ul_line_start[i + 3], p3);
}

static inline void
unsieve_7(unsigned char *line_start, const unsigned int start_idx,
          const unsigned int I, unsieve_aux_data_srcptr us)
{
  const unsigned int I_upt  = I / sizeof (unsieve_pattern_t);
  unsigned int i;
  unsieve_pattern_t p0, p1, p2, p3, p4, p5, p6;
  unsieve_pattern_t * ul_line_start = (unsieve_pattern_t *) line_start;
  size_t pattern_idx;

  if (sizeof(unsieve_pattern_t) == 4) {
    /* -4^(-1) == 5 (mod 7) */
    pattern_idx = (5 * start_idx) % 7;
  } else if (sizeof(unsieve_pattern_t) == 8) {
    /* -8^(-1) == -1 (mod 7) */
    pattern_idx = (7 - start_idx) % 7;
  } else if (sizeof(unsieve_pattern_t) == 16) {
    /* -16^(-1) == 3 (mod 7) */
    pattern_idx = (3 * start_idx) % 7;
  } else
    abort();
  
  p0 = us->pattern7[pattern_idx];
  p1 = us->pattern7[(pattern_idx + 1) % 7];
  p2 = us->pattern7[(pattern_idx + 2) % 7];
  p3 = us->pattern7[(pattern_idx + 3) % 7];
  p4 = us->pattern7[(pattern_idx + 4) % 7];
  p5 = us->pattern7[(pattern_idx + 5) % 7];
  p6 = us->pattern7[(pattern_idx + 6) % 7];

  if (start_idx < sizeof(p0)) {
      ASSERT_ALWAYS(((unsigned char *)&p0)[start_idx] == 255);
  } else {
      ASSERT_ALWAYS(start_idx < 2 * sizeof(p0));
      ASSERT_ALWAYS(((unsigned char *)&p1)[start_idx - sizeof(p0)] == 255);
  }

  /* Apply pattern to array */
  for (i = 0U; i < I_upt - 6U; i += 7U)
    {
      UNSIEVE_OR(ul_line_start[i], p0);
      UNSIEVE_OR(ul_line_start[i + 1], p1);
      UNSIEVE_OR(ul_line_start[i + 2], p2);
      UNSIEVE_OR(ul_line_start[i + 3], p3);
      UNSIEVE_OR(ul_line_start[i + 4], p4);
      UNSIEVE_OR(ul_line_start[i + 5], p5);
      UNSIEVE_OR(ul_line_start[i + 6], p6);
    }
  if (i < I_upt)
    UNSIEVE_OR(ul_line_start[i], p0);
  if (i + 1 < I_upt)
    UNSIEVE_OR(ul_line_start[i + 1], p1);
  if (i + 2 < I_upt)
    UNSIEVE_OR(ul_line_start[i + 2], p2);
  if (i + 3 < I_upt)
    UNSIEVE_OR(ul_line_start[i + 3], p3);
  if (i + 4 < I_upt)
    UNSIEVE_OR(ul_line_start[i + 4], p4);
  if (i + 5 < I_upt)
    UNSIEVE_OR(ul_line_start[i + 5], p5);
}


static void
unsieve_not_coprime_line(unsigned char * line_start,
                         const unsigned int j, const unsigned int min_p,
                         const unsigned int I, unsieve_aux_data_srcptr us)
{
  unsigned int p, start_idx, c = j;

  if (j == 0)
    return;

  while (c % 2U == 0U) 
    c >>= 1;

  while (1)
    {
      p = us->entries[c].lpf; /* set p to largest prime factor of c */
      start_idx = us->entries[c].start;
      c = us->entries[c].cof;
      if (p < min_p)
        return;
      if (p <= 7)
        break;
      unsieve_one_prime (line_start, p, j, start_idx, I);
    }
  
  if (p == 7U)
    {
      unsieve_7(line_start, start_idx, I, us);
      p = us->entries[c].lpf;
      start_idx = us->entries[c].start;
      c = us->entries[c].cof;
    }

  if (p < min_p)
    return;
  
  if (p == 5U)
    {
      unsieve_5(line_start, start_idx, I, us);
      p = us->entries[c].lpf;
      start_idx = us->entries[c].start;
      c = us->entries[c].cof;
    }

  if (p < min_p)
    return;
  
  if (p == 3U)
    {
      unsieve_3(line_start, start_idx, I, us);
    }
  ASSERT_ALWAYS(c <= 1);
}


j_div_srcptr
init_j_div(uint32_t J)
{
  /* Store largest prime factor of k in j_div[k].p, for 1 < k < J,
     and store 0 for k=0, 1 for k=1 */
  j_div_ptr j_div = (j_div_ptr) malloc (sizeof (struct j_div_s) * J);
  FATAL_ERROR_CHECK(j_div == NULL, "malloc failed");
  j_div[0].p   = 0U;
  j_div[0].cof = 0U;
  j_div[0].inv = 0U;
  j_div[0].bound = 0U;
  j_div[1].p   = 1U;
  j_div[1].cof = 1U;
  j_div[1].inv = 1U;
  j_div[1].inv = UINT_MAX;
  for (unsigned int k = 2U; k < J; k++) {
    /* Find largest prime factor of k */
    unsigned int p, c = k;
    for (p = 2U; p * p <= c; p += 1U + p % 2U)
      {
        while (c % p == 0U)
          c /= p;
        if (c == 1U)
          break;
      }
    p = (c == 1U) ? p : c;
    c = k; do {c /= p;} while (c % p == 0);
    j_div[k].p = p;
    j_div[k].cof = c;
    j_div[k].inv = p == 2 ? 0 : (unsigned int)ularith_invmod(p);
    j_div[k].bound = UINT_MAX / p;
  }
  return (j_div_srcptr) j_div;
}

void
clear_j_div(j_div_srcptr j_div)
{
  free((void *)j_div);
}

static inline unsigned int 
extract_j_div(unsigned int (*div)[2], const unsigned int j, j_div_srcptr j_div, 
              const unsigned int pmin, const unsigned int pmax)
{
    unsigned int c, nr_div = 0;
    /* For each distict odd prime factor p of j, if pmin <= p <= pmax,
       store the inverse and bound in array */
    c = j;
    c >>= ularith_ctz(c);
    while (c > 1) {
      unsigned int p = j_div[c].p;
      if (p < pmin)
          break;
      if (p <= pmax) {
          div[nr_div][0] = j_div[c].inv;
          div[nr_div++][1] = j_div[c].bound;
      }
      c = j_div[c].cof;
    }
    return nr_div;
}

#ifndef HAVE_SSE2
static inline int
sieve_info_test_lognorm (const unsigned char C1, const unsigned char C2,
                         const unsigned char S1, const unsigned char S2)
{
  return S1 <= C1 && S2 <= C2;
}
#endif

#ifdef HAVE_SSE2
static inline int 
sieve_info_test_lognorm_sse2(__m128i * S0, const __m128i pattern0,
                             const __m128i *S1, const __m128i pattern1)
{
    const __m128i zero = _mm_set1_epi8(0);
    const __m128i ff = _mm_set1_epi8(0xff);
    const __m128i sign_conversion = _mm_set1_epi8(-128);
    __m128i a = *S0;
    __m128i r = *S1;
    __m128i m1, m2;
    /* _mm_cmpgt_epi8() performs signed comparison, but we have unsigned
       bytes. We can switch to signed in a way that preserves ordering by
       flipping the MSB, e.g., 0xFF (255 unsigned) becomes 0x7F (+127 signed), 
       and 0x00 (0 unsigned) becomes 0x80 (-128 signed).

       If a byte in the first operand is greater than the corresponding byte in
       the second operand, the corresponding byte in the result is set to all 1s
       (i.e., to 0xFF); otherwise, it is set to all 0s.
       
       Normally, S[x] <= bound means a sieve survivor. However, for skipping over
       locations where 3 | gcd(i,j), we set the bound to 0 in the pattern at
       those locations. We then need a comparison that never produces a survivor
       in those locations, even when S[x] is in fact 0. Thus we initialise the
       pattern to bound + 1, then set the bound to 0 where 3 | gcd(i,j), and
       change the comparison to S[x] < bound, which is guaranteed not to let any
       survivors through where the pattern byte is 0. */
    m1 = _mm_cmpgt_epi8 (pattern0, _mm_xor_si128(a, sign_conversion));
    m2 = _mm_cmpgt_epi8 (pattern1, _mm_xor_si128(r, sign_conversion));
    /* m1 is 0xFF where pattern[x] > S0[x], i.e., where it survived.
       Same for S1. */
    /* Logically AND the two masks: survivor only where both sided survived */
    m1 = _mm_and_si128(m1, m2);

    /* m1 is 0xFF in those locations where there the sieve entry survived on
       both sides */
    /* For the OR mask we need the bit complement, via m1 XOR 0xFF */
    m2 = _mm_xor_si128(m1, ff);
    *S0 = _mm_or_si128(a, m2);
    /* Do we want to update this one? */
    // *S1 = _mm_or_si128(r, m2);

    /* Compute number of non-zero bytes. We want 1 is those bytes that
    survived, and 0 in the others. m1 has 0xFF in those bytes that
    survived, 0 in the others. First sign flip: 0xFF -> 0x1 */
    m1 = _mm_sub_epi8(zero, m1);
    /* Using Sum of Absolute Differences with 0, which for us gives the
    number of non-zero bytes. This Sum of Absolute Differences uses
    unsigned arithmetic, thus we needed the sign flip first */
    m1 = _mm_sad_epu8(m1, zero);
    /* Sum is stored in two parts */
    int nr_set = _mm_extract_epi16(m1, 0) + _mm_extract_epi16(m1, 4);
    /* Return number of bytes that were not set to 255 */
    return nr_set;
}
#endif


/* In SS[2][x_start] ... SS[2][x_start * x_step - 1], look for survivors.
   If SSE is available, the bounds check was already done. If no SSE2 is
   available, we still have to do it. In both cases, we still have to test
   divisibility of the resulting i value by the trial-divided primes.
   Return the number of survivors found. */
static inline int
search_single_survivors(unsigned char * const SS[2],
        const unsigned char bound[2] MAYBE_UNUSED, const unsigned int log_I,
        const unsigned int j, const int N MAYBE_UNUSED,
        const int x_start, const int x_step, const unsigned int nr_div,
        unsigned int (*div)[2])
{
  int survivors = 0;
  for (int x = x_start; x < x_start + x_step; x++) {
#ifndef HAVE_SSE2
      if (!sieve_info_test_lognorm(bound[0], bound[1], SS[0][x], SS[1][x]))
      {
          SS[0][x] = 255;
          continue;
      }
#else
      if (SS[0][x] == 255)
          continue;
#endif
      survivors++;

      /* The very small prime used in the bound pattern, and unsieving larger
         primes have not identified this as gcd(i,j) > 1. It remains to check
         the trial-divided primes. */
      const unsigned int i = abs (x - (1 << (log_I - 1)));
      int divides = 0;
      switch (nr_div) {
          case 6: divides |= (i * div[5][0] <= div[5][1]);
          case 5: divides |= (i * div[4][0] <= div[4][1]);
          case 4: divides |= (i * div[3][0] <= div[3][1]);
          case 3: divides |= (i * div[2][0] <= div[2][1]);
          case 2: divides |= (i * div[1][0] <= div[1][1]);
          case 1: divides |= (i * div[0][0] <= div[0][1]);
          case 0: while(0){};
      }

      if (divides)
      {
          if (verify_gcd)
              ASSERT_ALWAYS(bin_gcd_int64_safe (i, j) != 1);
#ifdef TRACE_K
          if (trace_on_spot_Nx(N, x)) {
              verbose_output_print(TRACE_CHANNEL, 0, "# Slot [%u] in bucket %u has non coprime (i,j)=(%d,%u)\n",
                      x, N, i, j);
          }
#endif
          SS[0][x] = 255;
      } else {
          if (verify_gcd)
              ASSERT_ALWAYS(bin_gcd_int64_safe (i, j) == 1);
#ifdef TRACE_K
          if (trace_on_spot_Nx(N, x)) {
              verbose_output_print(TRACE_CHANNEL, 0, "# Slot [%u] in bucket %u is survivor with coprime (i,j)\n",
                      x, N);
          }
#endif
      }
  }
  return survivors;
}


/* This function works for all j */
static int
search_survivors_in_line1(unsigned char * const SS[2],
        const unsigned char bound[2], const unsigned int log_I,
        const unsigned int j, const int N MAYBE_UNUSED, j_div_srcptr j_div,
        const unsigned int td_max)
{
    unsigned int div[6][2], nr_div;

    nr_div = extract_j_div(div, j, j_div, 3, td_max);
    ASSERT_ALWAYS(nr_div <= 6);

#ifdef HAVE_SSE2
    const __m128i sse2_sign_conversion = _mm_set1_epi8(-128);
    /* The reason for the bound+1 here is documented in
       sieve_info_test_lognorm_sse2() */
    __m128i patterns[2] = {
        _mm_xor_si128(_mm_set1_epi8(bound[0] + 1), sse2_sign_conversion),
        _mm_xor_si128(_mm_set1_epi8(bound[1] + 1), sse2_sign_conversion)
    };
    const int x_step = sizeof(__m128i);

    /* If j is even, set all the even entries on the bound pattern to
       unsigned 0 -> signed -128 */
    if (j % 2 == 0)
        for (size_t i = 0; i < sizeof(__m128i); i += 2)
            ((unsigned char *)&patterns[0])[i] = 0x80;
#else
    const int x_step = 1;
#endif
    int survivors = 0;

    for (int x_start = 0; x_start < (1 << log_I); x_start += x_step)
    {
#ifdef HAVE_SSE2
        /* Do bounds check using SSE pattern, set non-survivors in SS[0] array
           to 255 */
        int sse_surv = 
            sieve_info_test_lognorm_sse2((__m128i*) (SS[0] + x_start), patterns[0],
                                         (__m128i*) (SS[1] + x_start), patterns[1]);
        if (sse_surv == 0)
            continue;
#endif
        int surv = search_single_survivors(SS, bound, log_I, j, N, x_start,
            x_step, nr_div, div);
#ifdef HAVE_SSE2
        ASSERT(sse_surv == surv);
#endif
        survivors += surv;
    }
    return survivors;
}


/* This function assumes j % 3 == 0. It uses an SSE bound pattern where 
   i-coordinates with i % 3 == 0 are set to a bound of 0. */
#ifdef HAVE_SSE2
int
search_survivors_in_line3(unsigned char * const SS[2], 
        const unsigned char bound[2], const unsigned int log_I,
        const unsigned int j, const int N MAYBE_UNUSED, j_div_srcptr j_div,
        const unsigned int td_max)
{
    const __m128i sse2_sign_conversion = _mm_set1_epi8(-128);
    __m128i patterns[2][3];
    const int x_step = sizeof(__m128i);
    const int pmin = 5;
    int next_pattern = 0;
    int survivors = 0;

    /* We know that 3 does not divide j in the code of this function, and we
       don't store 2. Allowing 5 distinct odd prime factors >3 thus handles
       all j < 1616615, which is the smallest integer with 6 such factors */
    unsigned int div[5][2];
    unsigned int nr_div;

    nr_div = extract_j_div(div, j, j_div, pmin, td_max);
    ASSERT_ALWAYS(nr_div <= 5);

    patterns[0][0] = patterns[0][1] = patterns[0][2] = 
        _mm_xor_si128(_mm_set1_epi8(bound[0] + 1), sse2_sign_conversion);
    patterns[1][0] = patterns[1][1] = patterns[1][2] = 
        _mm_xor_si128(_mm_set1_epi8(bound[1] + 1), sse2_sign_conversion);

    if (j % 2 == 0)
        for (size_t i = 0; i < 3 * sizeof(__m128i); i += 2)
            ((unsigned char *)&patterns[0][0])[i] = 0x80;

    /* Those locations in patterns[0] that correspond to i being a multiple
       of 3 are set to 0. Byte 0 of patterns[0][0] corresponds to i = -I/2.
       We want d s.t. -I/2 + d == 0 (mod 3), or d == I/2 (mod 3). With
       I = 2^log_I and 2 == -1 (mod 3), we have d == -1^(log_I-1) (mod 3),
       or d = 2 if log_I is even and d = 1 if log_I is odd.
       We use the sign conversion trick (i.e., XOR 0x80), so to get an
       effective bound of unsigned 0, we need to set the byte to 0x80. */
    size_t d = 2 - log_I % 2;
    for (size_t i = 0; i < sizeof(__m128i); i++)
        ((unsigned char *)&patterns[0][0])[3*i + d] = 0x80;

    for (int x_start = 0; x_start < (1 << log_I); x_start += x_step)
    {
        int sse_surv = 
            sieve_info_test_lognorm_sse2((__m128i*) (SS[0] + x_start), patterns[0][next_pattern],
                                         (__m128i*) (SS[1] + x_start), patterns[1][next_pattern]);
        if (++next_pattern == 3)
            next_pattern = 0;
        if (sse_surv == 0)
            continue;
        int surv = search_single_survivors(SS, bound, log_I, j, N, x_start,
            x_step, nr_div, div);
        ASSERT(sse_surv == surv);
        survivors += surv;
    }
    return survivors;
}
#endif


#ifdef HAVE_SSE2
/* This function assumes j % 3 != 0 and j % 5 == 0. It uses an SSE bound 
   pattern where i-coordinates with i % 5 == 0 are set to a bound of 0,
   and trial divies only by primes > 5. */
int
search_survivors_in_line5(unsigned char * const SS[2], 
        const unsigned char bound[2], const unsigned int log_I,
        const unsigned int j, const int N MAYBE_UNUSED, 
        j_div_srcptr j_div, const unsigned int td_max)
{
    const __m128i sse2_sign_conversion = _mm_set1_epi8(-128);
    const int nr_patterns = 5;
    __m128i patterns[2][nr_patterns];
    const int x_step = sizeof(__m128i);
    const int pmin = 7;
    int next_pattern = 0;
    int survivors = 0;

    /* We know that 3 and 5 do not divide j in the code of this function, and
       we don't store 2. Allowing 5 distinct odd prime factors >5 thus handles
       all j < 7436429 */
    unsigned int div[5][2];
    unsigned int nr_div;

    nr_div = extract_j_div(div, j, j_div, pmin, td_max);
    ASSERT_ALWAYS(nr_div <= 5);

    for (int i = 0; i < nr_patterns; i++) {
        patterns[0][i] = _mm_xor_si128(_mm_set1_epi8(bound[0] + 1), sse2_sign_conversion);
        patterns[1][i] = _mm_xor_si128(_mm_set1_epi8(bound[1] + 1), sse2_sign_conversion);
    }

    if (j % 2 == 0)
        for (size_t i = 0; i < nr_patterns * sizeof(__m128i); i += 2)
            ((unsigned char *)&patterns[0][0])[i] = 0x80;

    /* Those locations in patterns[0] that correspond to i being a multiple
       of 5 are set to 0. Byte 0 of patterns[0][0] corresponds to i = -I/2.
       We want d s.t. -I/2 + d == 0 (mod 5), or d == I/2 (mod 5). With
       I = 2^log_I and ord_5(2) == 4 (mod 5), we have d == 2^((log_I-1)%4)
       (mod 5), so we want a function: 0->3, 1->1, 2->2, 3->4.
       We use the sign conversion trick (i.e., XOR 0x80), so to get an
       effective bound of unsigned 0, we need to set the byte to 0x80. */
    static const unsigned char d_lut[] = {3,1,2,4};
    size_t d = d_lut[log_I % 4];
    for (size_t i = 0; i < sizeof(__m128i); i++)
        ((unsigned char *)&patterns[0][0])[nr_patterns*i + d] = 0x80;

    for (int x_start = 0; x_start < (1 << log_I); x_start += x_step)
    {
        int sse_surv = 
            sieve_info_test_lognorm_sse2((__m128i*) (SS[0] + x_start), patterns[0][next_pattern],
                                         (__m128i*) (SS[1] + x_start), patterns[1][next_pattern]);
        if (++next_pattern == nr_patterns)
            next_pattern = 0;
        if (sse_surv == 0)
            continue;
        int surv = search_single_survivors(SS, bound, log_I, j, N, x_start,
            x_step, nr_div, div);
        ASSERT(sse_surv == surv);
        survivors += surv;
    }
    return survivors;
}
#endif


#define USE_PATTERN_3 1
#define USE_PATTERN_5 1

int
search_survivors_in_line(unsigned char * const SS[2], 
        const unsigned char bound[2], const unsigned int log_I,
        const unsigned int j, const int N, j_div_srcptr j_div,
        const unsigned int td_max, unsieve_aux_data_srcptr us)
{
    /* In line j = 0, only the coordinate (i, j) = (-1, 0) may survive */
    if (j == 0) {
        const size_t I = (size_t) 1 << log_I;
        const unsigned char s0 = SS[0][I / 2 - 1], s1 = SS[1][I / 2 - 1];
        memset(SS[0], 255, I);
        if (s0 <= bound[0] && s1 <= bound[1]) {
            SS[0][I / 2 - 1] = s0;
            SS[1][I / 2 - 1] = s1;
            return 1;
        } else {
            return 0;
        }
    }

    unsieve_not_coprime_line(SS[0], j, td_max + 1, 1U<<log_I, us);

#if defined(HAVE_SSE2) && USE_PATTERN_3
    if (j % 3 == 0)
      return search_survivors_in_line3(SS, bound, log_I, j, N, j_div, td_max);
#if defined(HAVE_SSE2) && USE_PATTERN_5
    else if (j % 5 == 0)
      return search_survivors_in_line5(SS, bound, log_I, j, N, j_div, td_max);
#endif
    else
#endif
      return search_survivors_in_line1(SS, bound, log_I, j, N, j_div, td_max);
}
