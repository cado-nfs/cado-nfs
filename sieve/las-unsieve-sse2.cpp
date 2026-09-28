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
#include "verbose.hpp"
#endif

#include "las-unsieve.hpp"
#include "macros.h"
#include "arith/ularith.h"
#include "gcd.h"

static const int verify_gcd = 0; /* Enable slow but thorough test */
static const __m128i sign_conversion = _mm_set1_epi8(-128);
/* Masks used to force a bound of zero -- i.e. "never a survivor" -- on the
 * positions whose real abscissa ii is even, on rows whose real jj is even.
 * Index 0 kills the even x, index 1 kills nothing, index 2 kills the odd x.
 * Which one applies is decided by sublat_coords::parity_skip_class(); note
 * that index 2 only ever arises with an odd sublattice modulus and an odd
 * sublat.i0, and index 1 covers both "jj is odd" and "the modulus is even,
 * so ii and jj are never both even". */
static const __m128i even_masks[3] = {
  _mm_set_epi8(0xFF, 0, 0xFF, 0, 0xFF, 0, 0xFF, 0,
               0xFF, 0, 0xFF, 0, 0xFF, 0, 0xFF, 0),
  _mm_set1_epi8(0xff),
  _mm_set_epi8(0, 0xFF, 0, 0xFF, 0, 0xFF, 0, 0xFF,
               0, 0xFF, 0, 0xFF, 0, 0xFF, 0, 0xFF)};

/* Pick the mask for a line, given its sublattice row index. */
static inline __m128i parity_mask_for_line(sublat_runtime_t const & S,
                                           unsigned int j)
{
    if (S.jj(j) % 2 != 0)
        return even_masks[1];          /* jj odd: nothing to exclude */
    switch (S.parity_skip_class()) {
        case 1: return even_masks[0];  /* exclude even x */
        case 2: return even_masks[2];  /* exclude odd x */
        default: return even_masks[1]; /* even modulus: nothing to exclude */
    }
}
/* The pattern-P variants (P = 3 or 5) kill, in bulk, the positions whose
 * real abscissa ii is a multiple of P. With
 *
 *      ii = m*(i0 + x) + sublat.i0
 *
 * that is still a single residue class of x modulo P whenever P does not
 * divide m, only shifted: x == -i0 - sublat.i0/m (mod P). This returns that
 * residue, or -1 when there is nothing to kill.
 *
 * The latter happens exactly when P divides m. Then ii == sublat.i0 (mod P)
 * over the whole class, and the caller has already established P | jj,
 * which for P dividing m means P | sublat.j0. A class with P dividing both
 * sublat.i0 and sublat.j0 is degenerate and never sieved, so P divides ii
 * nowhere and the P-pattern has no work to do.
 *
 * P is a template parameter on purpose: this runs once per line, and with a
 * runtime P the "% P" become real divisions and the inverse a search loop.
 * m is one of 1, 2, 3, 6, so the inverse is a table lookup. */
template<unsigned int P>
static inline int pattern_kill_offset(sublat_runtime_t const & S, int i0)
{
    static_assert(P == 3 || P == 5);
    constexpr unsigned char inv_mod_P[5] = { 0, 1, (P == 3) ? 2 : 3,
                                             (P == 3) ? 0 : 2, 4 };
    if (S.m == 1) {
        int const d = (int) ((-(int64_t) i0) % (int64_t) P);
        return d < 0 ? d + (int) P : d;
    }
    unsigned int const mi = S.m % P;
    if (mi == 0) {
        ASSERT(S.i0 % P != 0);
        return -1;
    }
    int64_t t = -(int64_t) i0 - (int64_t) (S.i0 % P) * (int64_t) inv_mod_P[mi];
    t %= (int64_t) P;
    return (int) (t < 0 ? t + (int64_t) P : t);
}

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
        std::vector<uint32_t> &survivors,
        sublat_runtime_t const & S)
{
  for (int x = x_start; UNLIKELY(bitmask != 0); x++) {
      const int tz = ularith_ctz(bitmask);
      x += tz;
      bitmask >>= tz + 1;

      /* The very small prime used in the bound pattern, and unsieving larger
         primes have not identified this as gcd(i,j) > 1. It remains to check
         the trial-divided primes. Note that div[] was extracted from the
         real row jj, so the value tested here has to be the real abscissa
         too. */
      const unsigned int i = S.abs_ii(i0, x);
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
              ASSERT_ALWAYS(bin_gcd_int64_safe (i, S.jj(j)) != 1);
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
              ASSERT_ALWAYS(bin_gcd_int64_safe (i, S.jj(j)) == 1);
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
        std::vector<uint32_t> &survivors,
        sublat_runtime_t const & S)
{
    unsigned int div[6][2], nr_div;

    /* the primes to test against are those of the *real* row */
    nr_div = extract_j_div(div, S.jj(j), j_div, 3, td_max);
    ASSERT_ALWAYS(nr_div <= 6);

    const __m128i even_mask = parity_mask_for_line(S, j);

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
            nr_div, div, mask, survivors, S);
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
        std::vector<uint32_t> &survivors,
        sublat_runtime_t const & S)
{
    unsigned int div[6][2], nr_div;

    nr_div = extract_j_div(div, S.jj(j), j_div, 3, td_max);
    ASSERT_ALWAYS(nr_div <= 6);

    const __m128i even_mask = parity_mask_for_line(S, j);

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
            nr_div, div, mask, survivors, S);
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
        std::vector<uint32_t> &survivors,
        sublat_runtime_t const & S)
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

    nr_div = extract_j_div(div, S.jj(j), j_div, pmin, td_max);
    ASSERT_ALWAYS(nr_div <= 5);

    /* If j is even, set all the even entries in the bound pattern to
       unsigned 0 */
    const __m128i even_mask = parity_mask_for_line(S, j);

    patterns[0][0] = patterns[0][1] = patterns[0][2] = 
        _mm_xor_si128(_mm_and_si128(_mm_set1_epi8(bound[0] + 1), even_mask), sign_conversion);
    patterns[1][0] = patterns[1][1] = patterns[1][2] = 
        _mm_xor_si128(_mm_and_si128(_mm_set1_epi8(bound[1] + 1), even_mask), sign_conversion);

    /* Those locations in patterns[0] that correspond to i being a multiple
     * of 3 are set to 0. Byte 0 of patterns[0][0] corresponds to i = i0.
     * We want d s.t. i0 + d == 0 (mod 3), or d == -i0 (mod 3).
     */
    const int d = pattern_kill_offset<3>(S, i0);
       
    /*
     * Special hack for i0=-(I/2):
     * I = 2^logI and 2 == -1 (mod 3), we have d == -1^(logI-1) (mod 3),
     * or d = 2 if logI is even and d = 1 if logI is odd.

         size_t d = 2 - logI % 2;
     */

    /* We use the sign conversion trick (i.e., XOR 0x80), so to get an
     * effective bound of unsigned 0, we need to set the byte to 0x80.
     */
    if (d >= 0)
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
            nr_div, div, mask, survivors, S);
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
        std::vector<uint32_t> &survivors,
        sublat_runtime_t const & S)
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

    nr_div = extract_j_div(div, S.jj(j), j_div, pmin, td_max);
    ASSERT_ALWAYS(nr_div <= 5);

    /* If j is even, set all the even entries in the bound pattern to
       unsigned 0 */
    const __m128i even_mask = parity_mask_for_line(S, j);

    patterns[0] = patterns[1] = patterns[2] = 
        _mm_xor_si128(_mm_and_si128(_mm_set1_epi8(bound + 1), even_mask), sign_conversion);

    /* Those locations in patterns[0] that correspond to i being a multiple
     * of 3 are set to 0. Byte 0 of patterns[0][0] corresponds to i = i0.
     * We want d s.t. i0 + d == 0 (mod 3), or d == -i0 (mod 3).
     */
    const int d = pattern_kill_offset<3>(S, i0);
       
    /*
     * Special hack for i0=-(I/2):
     * I = 2^logI and 2 == -1 (mod 3), we have d == -1^(logI-1) (mod 3),
     * or d = 2 if logI is even and d = 1 if logI is odd.

         size_t d = 2 - logI % 2;
     */

    /* We use the sign conversion trick (i.e., XOR 0x80), so to get an
     * effective bound of unsigned 0, we need to set the byte to 0x80.
     */
    if (d >= 0)
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
            nr_div, div, mask, survivors, S);
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
        std::vector<uint32_t> &survivors,
        sublat_runtime_t const & S)
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

    nr_div = extract_j_div(div, S.jj(j), j_div, pmin, td_max);
    ASSERT_ALWAYS(nr_div <= 5);

    /* If j is even, set all the even entries in the bound pattern to
       unsigned 0 */
    const __m128i even_mask = parity_mask_for_line(S, j);

    for (int i = 0; i < nr_patterns; i++) {
        patterns[0][i] = _mm_xor_si128(_mm_and_si128(_mm_set1_epi8(bound[0] + 1), even_mask), sign_conversion);
        patterns[1][i] = _mm_xor_si128(_mm_and_si128(_mm_set1_epi8(bound[1] + 1), even_mask), sign_conversion);
    }

    /* the AND with even_mask above has already forced the excluded
     * parity to a bound of zero, whichever parity that is. */

    /* Those locations in patterns[0] that correspond to i being a multiple
       of 5 are set to 0. Byte 0 of patterns[0][0] corresponds to i = i0.
       We want d s.t. i0 + d == 0 (mod 5), or d == i0 (mod 5).
     */
    const int d = pattern_kill_offset<5>(S, i0);

    /* Special trick for i0 = -(I/2) ; With
       I = 2^logI and ord_5(2) == 4 (mod 5), we have d == 2^((logI-1)%4)
       (mod 5), so we want a function: 0->3, 1->1, 2->2, 3->4.
    static const unsigned char d_lut[] = {3,1,2,4};
    size_t d = d_lut[logI % 4];
     */

    /* We use the sign conversion trick (i.e., XOR 0x80), so to get an
       effective bound of unsigned 0, we need to set the byte to 0x80. */
    if (d >= 0)
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
            nr_div, div, mask, survivors, S);
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
        std::vector<uint32_t> &survivors,
        sublat_runtime_t const & S)
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

    nr_div = extract_j_div(div, S.jj(j), j_div, pmin, td_max);
    ASSERT_ALWAYS(nr_div <= 5);

    /* If j is even, set all the even entries in the bound pattern to
       unsigned 0 */
    const __m128i even_mask = parity_mask_for_line(S, j);

    for (int i = 0; i < nr_patterns; i++) {
        patterns[i] = _mm_xor_si128(_mm_and_si128(_mm_set1_epi8(bound + 1), even_mask), sign_conversion);
    }

    /* the AND with even_mask above has already forced the excluded
     * parity to a bound of zero, whichever parity that is. */

    /* Those locations in patterns[0] that correspond to i being a multiple
       of 5 are set to 0. Byte 0 of patterns[0][0] corresponds to i = i0.
       We want d s.t. i0 + d == 0 (mod 5), or d == i0 (mod 5).
     */
    const int d = pattern_kill_offset<5>(S, i0);

    /* Special trick for i0 = -(I/2) ; With
       I = 2^logI and ord_5(2) == 4 (mod 5), we have d == 2^((logI-1)%4)
       (mod 5), so we want a function: 0->3, 1->1, 2->2, 3->4.
    static const unsigned char d_lut[] = {3,1,2,4};
    size_t d = d_lut[logI % 4];
     */

    /* We use the sign conversion trick (i.e., XOR 0x80), so to get an
       effective bound of unsigned 0, we need to set the byte to 0x80. */
    if (d >= 0)
      for (size_t i = 0; i < sizeof(__m128i); i++)
        ((unsigned char *)&patterns[0])[nr_patterns*i + d] = 0x80;

    for (int x_start = 0; x_start < (i1 - i0); x_start += x_step)
    {
        const unsigned int mask = sieve_info_test_lognorm_sse2_mask_oneside(
                (__m128i*) (SS + x_start), patterns[next_pattern]);
        if (++next_pattern == nr_patterns)
            next_pattern = 0;
        search_single_survivors_mask(SS, j, i0, i1, N, x_start,
            nr_div, div, mask, survivors, S);
    }
}


/* Whether to use the pattern-3 and pattern-5 variants at all.
 *
 * The patterns cost about 1.5 extra cycles per 16 bytes scanned -- the bound
 * pattern rotates through three (or five) values held on the stack, where the
 * generic variant keeps a single one in a register. What they buy is one
 * fewer entry in div[], i.e. one multiply-compare saved for each position
 * that passes the bound test. That only pays off when such positions are
 * dense, and in las they are not: the survivor ratio is around 1e-4, and
 * tests/sieve/test-unsieve-patterns puts the break-even point two orders of
 * magnitude above that.
 *
 * Measured on whole RSA-768 runs, with -falign-loops=64 so that code
 * placement does not drown the signal, the patterns cost 17 to 18% of
 * search_survivors on Intel Skylake-SP and 14 to 18% on AMD Zen 4, in every
 * configuration tried. So the default is to leave them alone. None of this
 * has a dramatic impact overall -- search_survivors is well under one percent
 * of a run -- but there is no reason to pay for it.
 *
 * -DUNSIEVE_PATTERNS=0/1 overrides the default, which is how one acts on what
 * the tuner reports on a microarchitecture we have not measured.
 */
#if defined(UNSIEVE_PATTERNS)
static constexpr bool use_unsieve_patterns = UNSIEVE_PATTERNS != 0;
#else
static constexpr bool use_unsieve_patterns = false;
#endif

bool const search_survivors_uses_patterns = use_unsieve_patterns;

/* The two functions below are templates only so that the tuner in
 * tests/sieve can time both settings against the code that las really
 * runs. Production code always goes through the wrappers that follow,
 * which hard-wire use_unsieve_patterns.
 */
template<bool use_patterns>
void
search_survivors_in_line_sse2_choice(unsigned char * const SS[2], 
        const unsigned char bound[2],
        unsigned int j,
        int i0, int i1,
        int N,
        j_divisibility_helper const & j_div,
        const unsigned int td_max, std::vector<uint32_t> &survivors,
        sublat_runtime_t sublat)
{
    if constexpr (use_patterns) {
        /* The patterns are indexed by x but select the positions whose *real*
         * abscissa is a multiple of 3 or 5, and the branch is on the real row.
         * pattern_kill_offset() carries the sublattice into the offset. */
        const unsigned int jj = sublat.jj(j);
        if (jj % 3 == 0) {
            search_survivors_in_line3_sse2(SS, bound, j, i0, i1, N, j_div,
                    td_max, survivors, sublat);
            return;
        }
        if (jj % 5 == 0) {
            search_survivors_in_line5_sse2(SS, bound, j, i0, i1, N, j_div,
                    td_max, survivors, sublat);
            return;
        }
    }
    search_survivors_in_line1_sse2(SS, bound, j, i0, i1, N, j_div,
            td_max, survivors, sublat);
}

template<bool use_patterns>
void
search_survivors_in_line_sse2_oneside_choice(unsigned char * const SS, 
        const unsigned char bound,
        unsigned int j,
        int i0, int i1,
        int N,
        j_divisibility_helper const & j_div,
        const unsigned int td_max, std::vector<uint32_t> &survivors,
        sublat_runtime_t sublat)
{
    if constexpr (use_patterns) {
        /* see the comment in search_survivors_in_line_sse2() */
        const unsigned int jj = sublat.jj(j);
        if (jj % 3 == 0) {
            search_survivors_in_line3_sse2_oneside(SS, bound, j, i0, i1, N, j_div,
                    td_max, survivors, sublat);
            return;
        }
        if (jj % 5 == 0) {
            search_survivors_in_line5_sse2_oneside(SS, bound, j, i0, i1, N, j_div,
                    td_max, survivors, sublat);
            return;
        }
    }
    search_survivors_in_line1_sse2_oneside(SS, bound, j, i0, i1, N, j_div,
            td_max, survivors, sublat);
}

/* Instantiated only for the tuner, which compiles this file itself. las gets
 * the single instantiation that the wrappers below ask for, and no extra
 * symbols: this code is placement-sensitive enough that it is not worth
 * moving it around for the sake of a test.
 */
#ifdef UNSIEVE_PATTERN_TUNER
template void search_survivors_in_line_sse2_choice<false>(
        unsigned char * const SS[2], const unsigned char bound[2],
        unsigned int, int, int, int, j_divisibility_helper const &,
        unsigned int, std::vector<uint32_t> &, sublat_runtime_t);
template void search_survivors_in_line_sse2_choice<true>(
        unsigned char * const SS[2], const unsigned char bound[2],
        unsigned int, int, int, int, j_divisibility_helper const &,
        unsigned int, std::vector<uint32_t> &, sublat_runtime_t);
template void search_survivors_in_line_sse2_oneside_choice<false>(
        unsigned char * const SS, unsigned char,
        unsigned int, int, int, int, j_divisibility_helper const &,
        unsigned int, std::vector<uint32_t> &, sublat_runtime_t);
template void search_survivors_in_line_sse2_oneside_choice<true>(
        unsigned char * const SS, unsigned char,
        unsigned int, int, int, int, j_divisibility_helper const &,
        unsigned int, std::vector<uint32_t> &, sublat_runtime_t);
#endif

void
search_survivors_in_line_sse2(unsigned char * const SS[2],
        const unsigned char bound[2],
        unsigned int j,
        int i0, int i1,
        int N,
        j_divisibility_helper const & j_div,
        const unsigned int td_max, std::vector<uint32_t> &survivors,
        sublat_runtime_t sublat)
{
    search_survivors_in_line_sse2_choice<use_unsieve_patterns>(SS, bound, j,
            i0, i1, N, j_div, td_max, survivors, sublat);
}

void
search_survivors_in_line_sse2_oneside(unsigned char * const SS,
        const unsigned char bound,
        unsigned int j,
        int i0, int i1,
        int N,
        j_divisibility_helper const & j_div,
        const unsigned int td_max, std::vector<uint32_t> &survivors,
        sublat_runtime_t sublat)
{
    search_survivors_in_line_sse2_oneside_choice<use_unsieve_patterns>(SS,
            bound, j, i0, i1, N, j_div, td_max, survivors, sublat);
}

void
search_survivors_in_line_sse2_siqs(
        unsigned char * SS,
        unsigned char bound,
        unsigned int length,
        std::vector<uint16_t> &survivors)
{
    __m128i const B = _mm_xor_si128(_mm_set1_epi8(bound+1), sign_conversion);
    const unsigned int x_step = sizeof(__m128i);

    for (unsigned int x_start = 0; x_start < length; x_start += x_step)
    {
        /* Do bounds check using SSE pattern, set non-survivors in SS[0] array
           to 255 */
        auto * ptrS = (__m128i *)(SS + x_start);
        __m128i const s = *ptrS;
        __m128i m = _mm_cmpgt_epi8(B, _mm_xor_si128(s, sign_conversion));
        auto bitmask = (unsigned int) _mm_movemask_epi8(m);
        m = _mm_xor_si128(m, ff);
        *ptrS = _mm_or_si128(s, m);

        for (unsigned int x = x_start; UNLIKELY(bitmask != 0); ++x) {
            unsigned int const tz = ularith_ctz(bitmask);
            x += tz;
            bitmask >>= tz + 1u;

            survivors.push_back(x);
        }
    }
}
