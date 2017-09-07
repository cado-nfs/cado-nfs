#include "cado.h"
#include <stdlib.h>
#include <string.h>

#include "las-config.h"
#include "las-types.hpp"
#include "las-unsieve.hpp"
#include "ularith.h"
#include "las-norms.hpp"
#include "las-debug.hpp"
#include "gcd.h"
#include "memory.h"

static const int verify_gcd = 0; /* Enable slow but thorough test */

/* Set every stride-th byte, starting at index 0, to 255 in an array of
   stride unsieve_data::pattern_t's, and set all other bytes to 0. */
static void
minisieve(unsieve_data::pattern_t * const array, const size_t stride)
{
    memset (array, 0, stride * sizeof(unsieve_data::pattern_t));
    for (size_t i = 0; i < stride * sizeof(unsieve_data::pattern_t); i += stride)
      ((unsigned char *) array)[i] = 255;
}

unsieve_data::unsieve_data()
{
  entries = NULL;
  Jmax = 0;
}

unsieve_data::unsieve_data(int logI, int logA)
{
  unsigned long I = 1UL << logI;
  ASSERT_ALWAYS(logI >= 0);
  ASSERT_ALWAYS(logA >= logI);
  Jmax = 1UL << (logA - logI);
  /* Store largest prime factor of k in us.lpf[k], 0 for k=0, 1 for k=1 */
  entries = new entry[Jmax];
  entries[0] = entry(0,0,0);
  entries[1] = entry(1,1,0);
  for (unsigned int k = 2U; k < Jmax; k++)
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
      entries[k] = entry(p, c, (I / 2U) % p);
    }

    minisieve(pattern3, 3);
    minisieve(pattern5, 5);
    minisieve(pattern7, 7);
}

unsieve_data::unsieve_data(unsieve_data const & o) : Jmax(o.Jmax)
{
    entries = NULL;
    if (Jmax == 0) return;
    entries = new entry[Jmax];
    memcpy(entries, o.entries, Jmax * sizeof(entry));
    memcpy(pattern3, o.pattern3, sizeof(pattern3));
    memcpy(pattern5, o.pattern5, sizeof(pattern5));
    memcpy(pattern7, o.pattern7, sizeof(pattern7));
}

unsieve_data & unsieve_data::operator=(unsieve_data const & o)
{
    if (Jmax) delete[] entries;
    Jmax = o.Jmax;
    entries = new entry[Jmax];
    memcpy(entries, o.entries, Jmax * sizeof(entry));
    memcpy(pattern3, o.pattern3, sizeof(pattern3));
    memcpy(pattern5, o.pattern5, sizeof(pattern5));
    memcpy(pattern7, o.pattern7, sizeof(pattern7));
    return *this;
}

unsieve_data::~unsieve_data()
{
    if (Jmax == 0) return;
    delete[] entries;
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
          const unsigned int I, unsieve_data const & us)
{
  const unsigned int I_upt  = I / sizeof (unsieve_data::pattern_t);
  unsigned int i, pattern_idx;
  unsieve_data::pattern_t p0, p1, p2;
  unsieve_data::pattern_t * ul_line_start = (unsieve_data::pattern_t *) line_start;

  if (sizeof(unsieve_data::pattern_t) == 4) {
    /* -4^(-1) == 2 (mod 3) */
    pattern_idx = (2 * start_idx) % 3;
  } else if (sizeof(unsieve_data::pattern_t) == 8) {
    /* -8^(-1) == 1 (mod 3) */
    pattern_idx = start_idx;
  } else if (sizeof(unsieve_data::pattern_t) == 16) {
    /* -16^(-1) == 2 (mod 3) */
    pattern_idx = (2 * start_idx) % 3;
  } else
    abort();
  
  p0 = us.pattern3[pattern_idx];
  p1 = us.pattern3[(pattern_idx + 1) % 3];
  p2 = us.pattern3[(pattern_idx + 2) % 3];

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
          const unsigned int I, unsieve_data const & us)
{
  const unsigned int I_upt  = I / sizeof (unsieve_data::pattern_t);
  unsigned int i;
  unsieve_data::pattern_t p0, p1, p2, p3, p4;
  unsieve_data::pattern_t * ul_line_start = (unsieve_data::pattern_t *) line_start;
  size_t pattern_idx;

  if (sizeof(unsieve_data::pattern_t) == 4) {
    /* -4^(-1) == 1 (mod 5) */
    pattern_idx = start_idx;
  } else if (sizeof(unsieve_data::pattern_t) == 8) {
    /* -8^(-1) == 3 (mod 5) */
    pattern_idx = (3 * start_idx) % 5;
  } else if (sizeof(unsieve_data::pattern_t) == 16) {
    /* -16^(-1) == -1 (mod 5) */
    pattern_idx = (5 - start_idx) % 5;
  } else
    abort();
  
  p0 = us.pattern5[pattern_idx];
  p1 = us.pattern5[(pattern_idx + 1) % 5];
  p2 = us.pattern5[(pattern_idx + 2) % 5];
  p3 = us.pattern5[(pattern_idx + 3) % 5];
  p4 = us.pattern5[(pattern_idx + 4) % 5];

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
          const unsigned int I, unsieve_data const & us)
{
  const unsigned int I_upt  = I / sizeof (unsieve_data::pattern_t);
  unsigned int i;
  unsieve_data::pattern_t p0, p1, p2, p3, p4, p5, p6;
  unsieve_data::pattern_t * ul_line_start = (unsieve_data::pattern_t *) line_start;
  size_t pattern_idx;

  if (sizeof(unsieve_data::pattern_t) == 4) {
    /* -4^(-1) == 5 (mod 7) */
    pattern_idx = (5 * start_idx) % 7;
  } else if (sizeof(unsieve_data::pattern_t) == 8) {
    /* -8^(-1) == -1 (mod 7) */
    pattern_idx = (7 - start_idx) % 7;
  } else if (sizeof(unsieve_data::pattern_t) == 16) {
    /* -16^(-1) == 3 (mod 7) */
    pattern_idx = (3 * start_idx) % 7;
  } else
    abort();
  
  p0 = us.pattern7[pattern_idx];
  p1 = us.pattern7[(pattern_idx + 1) % 7];
  p2 = us.pattern7[(pattern_idx + 2) % 7];
  p3 = us.pattern7[(pattern_idx + 3) % 7];
  p4 = us.pattern7[(pattern_idx + 4) % 7];
  p5 = us.pattern7[(pattern_idx + 5) % 7];
  p6 = us.pattern7[(pattern_idx + 6) % 7];

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
                         const unsigned int I, unsieve_data const & us)
{
  unsigned int p, start_idx, c = j;

  if (j == 0)
    return;

  while (c % 2U == 0U) 
    c >>= 1;

  while (1)
    {
      p = us.entries[c].lpf; /* set p to largest prime factor of c */
      start_idx = us.entries[c].start;
      c = us.entries[c].cof;
      if (p < min_p)
        return;
      if (p <= 7)
        break;
      unsieve_one_prime (line_start, p, j, start_idx, I);
    }
  
  if (p == 7U)
    {
      unsieve_7(line_start, start_idx, I, us);
      p = us.entries[c].lpf;
      start_idx = us.entries[c].start;
      c = us.entries[c].cof;
    }

  if (p < min_p)
    return;
  
  if (p == 5U)
    {
      unsieve_5(line_start, start_idx, I, us);
      p = us.entries[c].lpf;
      start_idx = us.entries[c].start;
      c = us.entries[c].cof;
    }

  if (p < min_p)
    return;
  
  if (p == 3U)
    {
      unsieve_3(line_start, start_idx, I, us);
    }
  ASSERT_ALWAYS(c <= 1);
}


j_divisibility_helper::j_divisibility_helper(uint32_t J) : J(J)
{
  /* Store largest prime factor of k in j_div[k].p, for 1 < k < J,
     and store 0 for k=0, 1 for k=1 */
  entries = NULL;
  if (!J) return;
  ASSERT_ALWAYS(J >= 2);
  entries = new entry[J];
  entries[0].p   = 0U;
  entries[0].cof = 0U;
  entries[0].inv = 0U;
  entries[0].bound = 0U;
  entries[1].p   = 1U;
  entries[1].cof = 1U;
  entries[1].inv = 1U;
  entries[1].inv = UINT_MAX;
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
    entries[k].p = p;
    entries[k].cof = c;
    entries[k].inv = p == 2 ? 0 : (unsigned int)ularith_invmod(p);
    entries[k].bound = UINT_MAX / p;
  }
}

j_divisibility_helper::j_divisibility_helper(j_divisibility_helper const & o) : J(o.J)
{
    entries = NULL;
    if (J == 0) return;
    entries = new entry[J];
    memcpy(entries, o.entries, J * sizeof(entry));
}

j_divisibility_helper & j_divisibility_helper::operator=(j_divisibility_helper const & o)
{
    if (J) delete[] entries;
    J = o.J;
    entries = new entry[J];
    memcpy(entries, o.entries, J * sizeof(entry));
    return *this;
}


j_divisibility_helper::~j_divisibility_helper()
{
    if (!J) return;
    delete[] entries;
}

static inline int
sieve_info_test_lognorm (const unsigned char C1, const unsigned char C2,
                         const unsigned char S1, const unsigned char S2)
{
  return S1 <= C1 && S2 <= C2;
}

/* In SS[2][x_start] ... SS[2][x_start * 2^log_I - 1], look for survivors.
   We test divisibility of the resulting i value by the trial-divided primes.
   Return the number of survivors found. This function works for all j */
MAYBE_UNUSED static void
search_survivors_in_line1(unsigned char * const SS[2],
        const unsigned char bound[2], const unsigned int log_I,
        const unsigned int j, const int N MAYBE_UNUSED, j_divisibility_helper const & j_div,
        const unsigned int td_max, std::vector<uint32_t> &survivors)
{
    unsigned int div[6][2], nr_div;

    nr_div = extract_j_div(div, j, j_div, 3, td_max);
    ASSERT_ALWAYS(nr_div <= 6);

    for (int x = 0; x < (1 << log_I); x++) {
        if (!sieve_info_test_lognorm(bound[0], bound[1], SS[0][x], SS[1][x]))
        {
            SS[0][x] = 255;
            continue;
        }

        /* The very small prime used in the bound pattern, and unsieving larger
           primes have not identified this as gcd(i,j) > 1. It remains to check
           the trial-divided primes. */
        const unsigned int i = abs (x - (1 << (log_I - 1)));
        int divides = 0;
        switch (nr_div) {
            case 6: divides |= (i * div[5][0] <= div[5][1]);no_break();
            case 5: divides |= (i * div[4][0] <= div[4][1]);no_break();
            case 4: divides |= (i * div[3][0] <= div[3][1]);no_break();
            case 3: divides |= (i * div[2][0] <= div[2][1]);no_break();
            case 2: divides |= (i * div[1][0] <= div[1][1]);no_break();
            case 1: divides |= (i * div[0][0] <= div[0][1]);no_break();
            case 0: break;
        }

        if (divides) {
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
            survivors.push_back(x);
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
}


void
search_survivors_in_line(unsigned char * const SS[2], 
        const unsigned char bound[2], const unsigned int log_I,
        const unsigned int j, const int N, j_divisibility_helper const & j_div,
        const unsigned int td_max, unsieve_data const & us,
        std::vector<uint32_t> &survivors, sublat_t sublat)
{
    /* In line j = 0, only the coordinate (i, j) = (-1, 0) may survive */
    // FIXME: in sublat mode, this is broken!
    if (j == 0 && (!sublat.m)) {
        const size_t I = (size_t) 1 << log_I;
        const unsigned char s0 = SS[0][I / 2 - 1], s1 = SS[1][I / 2 - 1];
        memset(SS[0], 255, I);
        if (s0 <= bound[0] && s1 <= bound[1]) {
            SS[0][I / 2 - 1] = s0;
            SS[1][I / 2 - 1] = s1;
            survivors.push_back(I / 2 - 1);
        } else {
            return;
        }
    }

    // Naive version when we have sublattices, because unsieving is
    // harder. TODO: implement a fast version
    if (sublat.m) {
        for (int x = 0; x < (1 << log_I); x++) {
            if (!sieve_info_test_lognorm(bound[0], bound[1], SS[0][x], SS[1][x])) {
                SS[0][x] = 255;
                continue;
            }
            const unsigned int i = abs (int(sublat.m)*(x - (1 << (log_I - 1)))+int(sublat.i0));
            const unsigned int jj = sublat.m*j+sublat.j0;
            if ((((jj % 2) == 0) && ((i % 2) == 0)) ||
                    (bin_gcd_int64_safe (i, jj) != 1)) {
                SS[0][x] = 255;
            } else {
                survivors.push_back(x);
            }
        }
        return;
    }

    unsieve_not_coprime_line(SS[0], j, td_max + 1, 1U<<log_I, us);

#if defined(HAVE_SSE2)
    search_survivors_in_line_sse2(SS, bound, log_I, j, N, j_div, td_max,
            survivors);
#else
    search_survivors_in_line1(SS, bound, log_I, j, N, j_div, td_max,
            survivors);
#endif
}
