#include "cado.h" // IWYU pragma: keep

/* This compilation units reacts to TRACK_CODE_PATH and uses macros
 * such as WHERE_AM_I_UPDATE.
 * This compilation unit _must_ produce different object files depending
 * on the value of TRACK_CODE_PATH.
 * The WHERE_AM_I_UPDATE macro itself is defined in las-where-am-i.hpp
 */

#ifndef HAVE_SSE2
#error "This file assumes SSE-2 support!"
#endif /* HAVE_SSE2 */

#include <cstdint>
#include <cstdlib>

#include <vector>

#include <emmintrin.h>

#ifdef TRACE_K
#include "las-where-am-i.hpp"
#include "las-output.hpp"
#include "verbose.h"
#endif

#include "las-unsieve.hpp"
#include "macros.h"
#include "arith/ularith.h"
#include "gcd.h"

static const int verify_gcd = 0; /* Enable slow but thorough test */
static const __m128i sign_conversion = _mm_set1_epi8(-128);
static const __m128i even_masks[2] = {
  _mm_set_epi8(0xFF, 0, 0xFF, 0, 0xFF, 0, 0xFF, 0,
               0xFF, 0, 0xFF, 0, 0xFF, 0, 0xFF, 0),
  _mm_set1_epi8(0xff)};
static const __m128i ff = _mm_set1_epi8(0xff);

static inline unsigned int
sieve_info_test_lognorm_sse2_mask(__m128i * S0, const __m128i pattern0,
                             const __m128i *S1, const __m128i pattern1)
{
    __m128i const a = *S0;
    __m128i const r = *S1;
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
    /* Logically AND the two masks: survivor only where both sides survived */
    m1 = _mm_and_si128(m1, m2);

    /* m1 is 0xFF in those locations where the sieve entry survived on
       both sides */
    /* For the OR mask we need the bit complement, via m1 XOR 0xFF */
    m2 = _mm_xor_si128(m1, ff);
    *S0 = _mm_or_si128(a, m2);
    /* Do we want to update this one? */
    // *S1 = _mm_or_si128(r, m2);

    /* Compute mask of non-zero bytes */
    return (unsigned int) _mm_movemask_epi8(m1);
}

static inline unsigned int
sieve_info_test_lognorm_sse2_mask_oneside(__m128i * S0, const __m128i pattern0)
{
    __m128i const a = *S0;
    __m128i const m1 = _mm_cmpgt_epi8 (pattern0, _mm_xor_si128(a, sign_conversion));
    *S0 = _mm_or_si128(a, _mm_xor_si128(m1, ff));
    return (unsigned int) _mm_movemask_epi8(m1);
}

/* Look for survivors as indicated by a bit mask.
   We still have to test divisibility of the resulting i value by the
   trial-divided primes. Return the number of survivors found. */
static inline void
search_single_survivors_mask(unsigned char * const SS,
        unsigned int j,
        int i0,
        int i1 MAYBE_UNUSED,
        int N MAYBE_UNUSED,
        int x_start,
        unsigned int nr_div,
        unsigned int (*div)[2],
        unsigned int bitmask,
        std::vector<uint32_t> &survivors)
{
  for (int x = x_start; UNLIKELY(bitmask != 0); x++) {
      const int tz = ularith_ctz(bitmask);
      x += tz;
      bitmask >>= tz + 1;

      /* The very small prime used in the bound pattern, and unsieving larger
         primes have not identified this as gcd(i,j) > 1. It remains to check
         the trial-divided primes. */
      const unsigned int i = abs (i0 + x);
      int divides = 0;
      switch (nr_div) {
            // coverity[unterminated_case]
          case 6: divides |= (i * div[5][0] <= div[5][1]); no_break();
            // coverity[unterminated_case]
          case 5: divides |= (i * div[4][0] <= div[4][1]); no_break();
            // coverity[unterminated_case]
          case 4: divides |= (i * div[3][0] <= div[3][1]); no_break();
            // coverity[unterminated_case]
          case 3: divides |= (i * div[2][0] <= div[2][1]); no_break();
            // coverity[unterminated_case]
          case 2: divides |= (i * div[1][0] <= div[1][1]); no_break();
            // coverity[unterminated_case]
          case 1: divides |= (i * div[0][0] <= div[0][1]); no_break();
          case 0: while(0){};
      }

      if (divides)
      {
          if (verify_gcd)
              ASSERT_ALWAYS(bin_gcd_int64_safe (i, j) != 1);
#ifdef TRACE_K
          if (trace_on_spot_Nx(N, x)) {
              verbose_fmt_print(TRACE_CHANNEL, 0, "# Slot [{}] in bucket {} has non coprime (i,j)=({},{})\n",
                      x, N, i, j);
          }
#endif
          SS[x] = 255;
      } else {
          survivors.push_back(x);
          if (verify_gcd)
              ASSERT_ALWAYS(bin_gcd_int64_safe (i, j) == 1);
#ifdef TRACE_K
          if (trace_on_spot_Nx(N, x)) {
              verbose_fmt_print(TRACE_CHANNEL, 0, "# Slot [{}] in bucket {} is survivor with coprime (i,j)\n",
                      x, N);
          }
#endif
      }
  }
}

/* This function works for all j. Uses SSE2. */
static void
search_survivors_in_line1_sse2(unsigned char * const SS[2],
        const unsigned char bound[2],
        unsigned int j,
        int i0, int i1,
        int N MAYBE_UNUSED,
        j_divisibility_helper const & j_div,
        unsigned int td_max,
        std::vector<uint32_t> &survivors)
{
    unsigned int div[6][2], nr_div;

    nr_div = extract_j_div(div, j, j_div, 3, td_max);
    ASSERT_ALWAYS(nr_div <= 6);

    /* If j is even, set all the even entries in the bound pattern to
       unsigned 0 */
    const __m128i even_mask = even_masks[j % 2];

    /* The reason for the bound+1 here is documented in
       sieve_info_test_lognorm_sse2_mask() */
    __m128i const patterns[2] = {
        _mm_xor_si128(_mm_and_si128(_mm_set1_epi8(bound[0] + 1), even_mask), sign_conversion),
        _mm_xor_si128(_mm_and_si128(_mm_set1_epi8(bound[1] + 1), even_mask), sign_conversion)
    };
    const int x_step = sizeof(__m128i);

    for (int x_start = 0; x_start < (i1 - i0); x_start += x_step)
    {
        /* Do bounds check using SSE pattern, set non-survivors in SS[0] array
           to 255 */
        const unsigned int mask = sieve_info_test_lognorm_sse2_mask(
                    (__m128i*) (SS[0] + x_start), patterns[0],
                    (__m128i*) (SS[1] + x_start), patterns[1]);
        search_single_survivors_mask(SS[0], j, i0, i1, N, x_start,
            nr_div, div, mask, survivors);
    }
}

static void
search_survivors_in_line1_sse2_oneside(unsigned char * SS,
        const unsigned char bound,
        unsigned int j,
        int i0, int i1,
        int N MAYBE_UNUSED,
        j_divisibility_helper const & j_div,
        unsigned int td_max,
        std::vector<uint32_t> &survivors)
{
    unsigned int div[6][2], nr_div;

    nr_div = extract_j_div(div, j, j_div, 3, td_max);
    ASSERT_ALWAYS(nr_div <= 6);

    /* If j is even, set all the even entries in the bound pattern to
       unsigned 0 */
    const __m128i even_mask = even_masks[j % 2];

    /* The reason for the bound+1 here is documented in
       sieve_info_test_lognorm_sse2_mask() */
    __m128i const pattern =
        _mm_xor_si128(_mm_and_si128(_mm_set1_epi8(bound + 1), even_mask), sign_conversion);
    const int x_step = sizeof(__m128i);

    for (int x_start = 0; x_start < (i1 - i0); x_start += x_step)
    {
        /* Do bounds check using SSE pattern, set non-survivors in SS[0] array
           to 255 */
        const unsigned int mask = sieve_info_test_lognorm_sse2_mask_oneside(
                    (__m128i*) (SS + x_start), pattern);
        search_single_survivors_mask(SS, j, i0, i1, N, x_start,
            nr_div, div, mask, survivors);
    }
}

/* This function assumes j % 3 == 0. It uses an SSE bound pattern where 
   i-coordinates with i % 3 == 0 are set to a bound of 0. */
static void
search_survivors_in_line3_sse2(unsigned char * const SS[2], 
        const unsigned char bound[2], 
        unsigned int j,
        int i0, int i1,
        int N MAYBE_UNUSED,
        j_divisibility_helper const & j_div,
        unsigned int td_max,
        std::vector<uint32_t> &survivors)
{
    __m128i patterns[2][3];
    const int x_step = sizeof(__m128i);
    const int pmin = 5;
    int next_pattern = 0;

    /* We know that 3 does not divide j in the code of this function, and we
       don't store 2. Allowing 5 distinct odd prime factors >3 thus handles
       all j < 1616615, which is the smallest integer with 6 such factors */
    unsigned int div[5][2];
    unsigned int nr_div;

    nr_div = extract_j_div(div, j, j_div, pmin, td_max);
    ASSERT_ALWAYS(nr_div <= 5);

    /* If j is even, set all the even entries in the bound pattern to
       unsigned 0 */
    const __m128i even_mask = even_masks[j % 2];

    patterns[0][0] = patterns[0][1] = patterns[0][2] = 
        _mm_xor_si128(_mm_and_si128(_mm_set1_epi8(bound[0] + 1), even_mask), sign_conversion);
    patterns[1][0] = patterns[1][1] = patterns[1][2] = 
        _mm_xor_si128(_mm_and_si128(_mm_set1_epi8(bound[1] + 1), even_mask), sign_conversion);

    /* Those locations in patterns[0] that correspond to i being a multiple
     * of 3 are set to 0. Byte 0 of patterns[0][0] corresponds to i = i0.
     * We want d s.t. i0 + d == 0 (mod 3), or d == -i0 (mod 3).
     */
    int d = (-i0) % 3;
    if (d < 0) d += 3;
       
    /*
     * Special hack for i0=-(I/2):
     * I = 2^logI and 2 == -1 (mod 3), we have d == -1^(logI-1) (mod 3),
     * or d = 2 if logI is even and d = 1 if logI is odd.

         size_t d = 2 - logI % 2;
     */

    /* We use the sign conversion trick (i.e., XOR 0x80), so to get an
     * effective bound of unsigned 0, we need to set the byte to 0x80.
     */
    for (size_t i = 0; i < sizeof(__m128i); i++)
        ((unsigned char *)&patterns[0][0])[3*i + d] = 0x80;

    for (int x_start = 0; x_start < (i1 - i0); x_start += x_step)
    {
        const unsigned int mask =
            sieve_info_test_lognorm_sse2_mask(
                    (__m128i*) (SS[0] + x_start), patterns[0][next_pattern],
                    (__m128i*) (SS[1] + x_start), patterns[1][next_pattern]);
        if (++next_pattern == 3)
            next_pattern = 0;
        search_single_survivors_mask(SS[0], j, i0, i1, N, x_start,
            nr_div, div, mask, survivors);
    }
}

/* This function assumes j % 3 == 0. It uses an SSE bound pattern where 
   i-coordinates with i % 3 == 0 are set to a bound of 0. */
static void
search_survivors_in_line3_sse2_oneside(unsigned char * const SS, 
        const unsigned char bound, 
        unsigned int j,
        int i0, int i1,
        int N MAYBE_UNUSED,
        j_divisibility_helper const & j_div,
        unsigned int td_max,
        std::vector<uint32_t> &survivors)
{
    __m128i patterns[3];
    const int x_step = sizeof(__m128i);
    const int pmin = 5;
    int next_pattern = 0;

    /* We know that 3 does not divide j in the code of this function, and we
       don't store 2. Allowing 5 distinct odd prime factors >3 thus handles
       all j < 1616615, which is the smallest integer with 6 such factors */
    unsigned int div[5][2];
    unsigned int nr_div;

    nr_div = extract_j_div(div, j, j_div, pmin, td_max);
    ASSERT_ALWAYS(nr_div <= 5);

    /* If j is even, set all the even entries in the bound pattern to
       unsigned 0 */
    const __m128i even_mask = even_masks[j % 2];

    patterns[0] = patterns[1] = patterns[2] = 
        _mm_xor_si128(_mm_and_si128(_mm_set1_epi8(bound + 1), even_mask), sign_conversion);

    /* Those locations in patterns[0] that correspond to i being a multiple
     * of 3 are set to 0. Byte 0 of patterns[0][0] corresponds to i = i0.
     * We want d s.t. i0 + d == 0 (mod 3), or d == -i0 (mod 3).
     */
    int d = (-i0) % 3;
    if (d < 0) d += 3;
       
    /*
     * Special hack for i0=-(I/2):
     * I = 2^logI and 2 == -1 (mod 3), we have d == -1^(logI-1) (mod 3),
     * or d = 2 if logI is even and d = 1 if logI is odd.

         size_t d = 2 - logI % 2;
     */

    /* We use the sign conversion trick (i.e., XOR 0x80), so to get an
     * effective bound of unsigned 0, we need to set the byte to 0x80.
     */
    for (size_t i = 0; i < sizeof(__m128i); i++)
        ((unsigned char *)&patterns[0])[3*i + d] = 0x80;

    for (int x_start = 0; x_start < (i1 - i0); x_start += x_step)
    {
        const unsigned int mask =
            sieve_info_test_lognorm_sse2_mask_oneside(
                    (__m128i*) (SS + x_start), patterns[next_pattern]);
        if (++next_pattern == 3)
            next_pattern = 0;
        search_single_survivors_mask(SS, j, i0, i1, N, x_start,
            nr_div, div, mask, survivors);
    }
}


/* This function assumes j % 3 != 0 and j % 5 == 0. It uses an SSE bound 
   pattern where i-coordinates with i % 5 == 0 are set to a bound of 0,
   and trial divides only by primes > 5. */
static void
search_survivors_in_line5_sse2(unsigned char * const SS[2], 
        const unsigned char bound[2],
        unsigned int j,
        int i0, int i1,
        int N MAYBE_UNUSED,
        j_divisibility_helper const & j_div,
        unsigned int td_max,
        std::vector<uint32_t> &survivors)
{
    const int nr_patterns = 5;
    __m128i patterns[2][nr_patterns];
    const int x_step = sizeof(__m128i);
    const int pmin = 7;
    int next_pattern = 0;

    /* We know that 3 and 5 do not divide j in the code of this function, and
       we don't store 2. Allowing 5 distinct odd prime factors >5 thus handles
       all j < 7436429 */
    unsigned int div[5][2];
    unsigned int nr_div;

    nr_div = extract_j_div(div, j, j_div, pmin, td_max);
    ASSERT_ALWAYS(nr_div <= 5);

    /* If j is even, set all the even entries in the bound pattern to
       unsigned 0 */
    const __m128i even_mask = even_masks[j % 2];

    for (int i = 0; i < nr_patterns; i++) {
        patterns[0][i] = _mm_xor_si128(_mm_and_si128(_mm_set1_epi8(bound[0] + 1), even_mask), sign_conversion);
        patterns[1][i] = _mm_xor_si128(_mm_and_si128(_mm_set1_epi8(bound[1] + 1), even_mask), sign_conversion);
    }

    if (j % 2 == 0)
        for (size_t i = 0; i < nr_patterns * sizeof(__m128i); i += 2)
            ((unsigned char *)&patterns[0][0])[i] = 0x80;

    /* Those locations in patterns[0] that correspond to i being a multiple
       of 5 are set to 0. Byte 0 of patterns[0][0] corresponds to i = i0.
       We want d s.t. i0 + d == 0 (mod 5), or d == i0 (mod 5).
     */
    int d = (-i0) % 5;
    if (d < 0) d += 5;

    /* Special trick for i0 = -(I/2) ; With
       I = 2^logI and ord_5(2) == 4 (mod 5), we have d == 2^((logI-1)%4)
       (mod 5), so we want a function: 0->3, 1->1, 2->2, 3->4.
    static const unsigned char d_lut[] = {3,1,2,4};
    size_t d = d_lut[logI % 4];
     */

    /* We use the sign conversion trick (i.e., XOR 0x80), so to get an
       effective bound of unsigned 0, we need to set the byte to 0x80. */
    for (size_t i = 0; i < sizeof(__m128i); i++)
        ((unsigned char *)&patterns[0][0])[nr_patterns*i + d] = 0x80;

    for (int x_start = 0; x_start < (i1 - i0); x_start += x_step)
    {
        const unsigned int mask = sieve_info_test_lognorm_sse2_mask(
                (__m128i*) (SS[0] + x_start), patterns[0][next_pattern],
                (__m128i*) (SS[1] + x_start), patterns[1][next_pattern]);
        if (++next_pattern == nr_patterns)
            next_pattern = 0;
        search_single_survivors_mask(SS[0], j, i0, i1, N, x_start,
            nr_div, div, mask, survivors);
    }
}

/* This function assumes j % 3 != 0 and j % 5 == 0. It uses an SSE bound 
   pattern where i-coordinates with i % 5 == 0 are set to a bound of 0,
   and trial divides only by primes > 5. */
static void
search_survivors_in_line5_sse2_oneside(unsigned char * const SS, 
        const unsigned char bound,
        unsigned int j,
        int i0, int i1,
        int N MAYBE_UNUSED,
        j_divisibility_helper const & j_div,
        unsigned int td_max,
        std::vector<uint32_t> &survivors)
{
    const int nr_patterns = 5;
    __m128i patterns[nr_patterns];
    const int x_step = sizeof(__m128i);
    const int pmin = 7;
    int next_pattern = 0;

    /* We know that 3 and 5 do not divide j in the code of this function, and
       we don't store 2. Allowing 5 distinct odd prime factors >5 thus handles
       all j < 7436429 */
    unsigned int div[5][2];
    unsigned int nr_div;

    nr_div = extract_j_div(div, j, j_div, pmin, td_max);
    ASSERT_ALWAYS(nr_div <= 5);

    /* If j is even, set all the even entries in the bound pattern to
       unsigned 0 */
    const __m128i even_mask = even_masks[j % 2];

    for (int i = 0; i < nr_patterns; i++) {
        patterns[i] = _mm_xor_si128(_mm_and_si128(_mm_set1_epi8(bound + 1), even_mask), sign_conversion);
    }

    if (j % 2 == 0)
        for (size_t i = 0; i < nr_patterns * sizeof(__m128i); i += 2)
            ((unsigned char *)&patterns[0])[i] = 0x80;

    /* Those locations in patterns[0] that correspond to i being a multiple
       of 5 are set to 0. Byte 0 of patterns[0][0] corresponds to i = i0.
       We want d s.t. i0 + d == 0 (mod 5), or d == i0 (mod 5).
     */
    int d = (-i0) % 5;
    if (d < 0) d += 5;

    /* Special trick for i0 = -(I/2) ; With
       I = 2^logI and ord_5(2) == 4 (mod 5), we have d == 2^((logI-1)%4)
       (mod 5), so we want a function: 0->3, 1->1, 2->2, 3->4.
    static const unsigned char d_lut[] = {3,1,2,4};
    size_t d = d_lut[logI % 4];
     */

    /* We use the sign conversion trick (i.e., XOR 0x80), so to get an
       effective bound of unsigned 0, we need to set the byte to 0x80. */
    for (size_t i = 0; i < sizeof(__m128i); i++)
        ((unsigned char *)&patterns[0])[nr_patterns*i + d] = 0x80;

    for (int x_start = 0; x_start < (i1 - i0); x_start += x_step)
    {
        const unsigned int mask = sieve_info_test_lognorm_sse2_mask_oneside(
                (__m128i*) (SS + x_start), patterns[next_pattern]);
        if (++next_pattern == nr_patterns)
            next_pattern = 0;
        search_single_survivors_mask(SS, j, i0, i1, N, x_start,
            nr_div, div, mask, survivors);
    }
}


#define USE_PATTERN_3 1
#define USE_PATTERN_5 1
#if USE_PATTERN_5 && ! USE_PATTERN_3
#error USE_PATTERN_5 requires USE_PATTERN_3
#endif

void
search_survivors_in_line_sse2(unsigned char * const SS[2], 
        const unsigned char bound[2],
        unsigned int j,
        int i0, int i1,
        int N,
        j_divisibility_helper const & j_div,
        const unsigned int td_max, std::vector<uint32_t> &survivors)
{
#if USE_PATTERN_3
    if (j % 3 == 0)
      search_survivors_in_line3_sse2(SS, bound, j, i0, i1, N, j_div,
              td_max, survivors);
#if USE_PATTERN_5
    else if (j % 5 == 0)
      search_survivors_in_line5_sse2(SS, bound, j, i0, i1, N, j_div,
              td_max, survivors);
#endif
    else
#endif
      search_survivors_in_line1_sse2(SS, bound, j, i0, i1, N, j_div,
              td_max, survivors);
}

void
search_survivors_in_line_sse2_oneside(unsigned char * const SS, 
        const unsigned char bound,
        unsigned int j,
        int i0, int i1,
        int N,
        j_divisibility_helper const & j_div,
        const unsigned int td_max, std::vector<uint32_t> &survivors)
{
#if USE_PATTERN_3
    if (j % 3 == 0)
      search_survivors_in_line3_sse2_oneside(SS, bound, j, i0, i1, N, j_div,
              td_max, survivors);
#if USE_PATTERN_5
    else if (j % 5 == 0)
      search_survivors_in_line5_sse2_oneside(SS, bound, j, i0, i1, N, j_div,
              td_max, survivors);
#endif
    else
#endif
      search_survivors_in_line1_sse2_oneside(SS, bound, j, i0, i1, N, j_div,
              td_max, survivors);
}
