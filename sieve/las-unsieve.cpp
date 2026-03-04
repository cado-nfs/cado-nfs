#include "cado.h" // IWYU pragma: keep

/* This compilation units reacts to TRACK_CODE_PATH and uses macros
 * such as WHERE_AM_I_UPDATE.
 * This compilation unit _must_ produce different object files depending
 * on the value of TRACK_CODE_PATH.
 * The WHERE_AM_I_UPDATE macro itself is defined in las-where-am-i.hpp
 */

#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "fb-types.hpp"
#include "gcd.h"
#include "gpf.h"
#include "las-unsieve.hpp"
#include "macros.h"
#include "arith/ularith.h"
#ifdef TRACE_K
#include "las-where-am-i.hpp"
#include "las-output.hpp"
#include "verbose.h"
#endif

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
  ASSERT_ALWAYS(logI >= 0);
  ASSERT_ALWAYS(logA >= logI);
  Jmax = 1UL << (logA - logI);
  /* Store largest prime factor of k in us.lpf[k], 0 for k=0, 1 for k=1 */
  entries = new entry[Jmax];
  entries[0] = entry(0,0);
  entries[1] = entry(1,1);
  gpf_init(Jmax - 1);
  for (unsigned int k = 2U; k < Jmax; k++)
    {
      unsigned int p, c;
      p = gpf_get(k);
      c = k; do {c /= p;} while (c % p == 0);
      entries[k] = entry(p, c);
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
                         unsigned int j,
                         int i0, int i1,
                         unsigned int min_p,
                         unsieve_data const & us)
{
  unsigned int p, c=j;
  int start_idx;

  if (j == 0)
    return;

  while (c % 2U == 0U) 
    c >>= 1;

  while (1)
    {
      p = us.entries[c].lpf; /* set p to largest prime factor of c */
      if (p < min_p)
        return;
      start_idx = (-i0) % p; if (start_idx < 0) start_idx += p;
      c = us.entries[c].cof;
      if (p <= 7)
        break;
      unsieve_one_prime (line_start, p, j, start_idx, i1 - i0);
    }
  
  if (p == 7U)
    {
      unsieve_7(line_start, start_idx, i1 - i0, us);
      p = us.entries[c].lpf;
      start_idx = (-i0) % p; if (start_idx < 0) start_idx += p;
      c = us.entries[c].cof;
    }

  if (p < min_p)
    return;
  
  if (p == 5U)
    {
      unsieve_5(line_start, start_idx, i1 - i0, us);
      p = us.entries[c].lpf;
      start_idx = (-i0) % p; if (start_idx < 0) start_idx += p;
      c = us.entries[c].cof;
    }

  if (p < min_p)
    return;
  
  if (p == 3U)
      unsieve_3(line_start, start_idx, i1 - i0, us);

  ASSERT_ALWAYS(c <= 1);
}

j_divisibility_helper::j_divisibility_helper(uint32_t J)
{
    /* Store largest prime factor of k in j_div[k].p, for 1 < k < J,
       and store 0 for k=0, 1 for k=1 */
    ASSERT_ALWAYS(J >= 2);
    ASSERT_ALWAYS(!(J & (J-1)));        /* J must be a power of two */
    entries.reserve(J);
    entries.push_back({ 0U, 0U, 0U, 0U });
    entries.push_back({ 1U, 1U, 1U, UINT_MAX });
    for (unsigned int k = 2; k < J; k++) {
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
        unsigned int const inv = p == 2 ? 0U : (unsigned int)ularith_invmod(p);
        unsigned int const bound = UINT_MAX / p;
        entries.push_back({ p, c, inv, bound});
    }
}

template<std::size_t nsides>
static inline int
sieve_info_test_lognorm (
        const std::array<unsigned char, nsides> bound,
        std::array<unsigned char * const, nsides> SS,
        const unsigned int x);

template<>
inline int
sieve_info_test_lognorm (
        const std::array<unsigned char, 1> bound,
        std::array<unsigned char * const, 1> SS,
        const unsigned int x)
{
  return SS[0][x] <= bound[0];
}

template<>
inline int
sieve_info_test_lognorm (
        const std::array<unsigned char, 2> bound,
        std::array<unsigned char * const, 2> SS,
        const unsigned int x)
{
  return SS[0][x] <= bound[0] && SS[1][x] <= bound[1];
}

static inline int
sieve_info_test_lognorm (const unsigned char C1, const unsigned char C2,
                         const unsigned char S1, const unsigned char S2)
{
  return S1 <= C1 && S2 <= C2;
}

/* In SS[2][x_start] ... SS[2][x_start * 2^logI - 1], look for survivors.
   We test divisibility of the resulting i value by the trial-divided primes.
   Return the number of survivors found. This function works for all j */
MAYBE_UNUSED static void
search_survivors_in_line1(unsigned char * const SS[2],
        const unsigned char bound[2],
        unsigned int j, 
        int i0, int i1,
        int N MAYBE_UNUSED, j_divisibility_helper const & j_div,
        unsigned int td_max, std::vector<uint32_t> &survivors)
{
    unsigned int div[6][2], nr_div;

    nr_div = extract_j_div(div, j, j_div, 3, td_max);
    ASSERT_ALWAYS(nr_div <= 6);

    for (int x = 0; x < (i1 - i0); x++) {
        if (!sieve_info_test_lognorm(bound[0], bound[1], SS[0][x], SS[1][x]))
        {
            SS[0][x] = 255;
            continue;
        }

        /* The very small prime used in the bound pattern, and unsieving larger
           primes have not identified this as gcd(i,j) > 1. It remains to check
           the trial-divided primes. */
        const unsigned int i = abs (i0 + x);
        int divides = 0;
        switch (nr_div) {
            // coverity[unterminated_case]
            case 6: divides |= (i * div[5][0] <= div[5][1]);no_break();
            // coverity[unterminated_case]
            case 5: divides |= (i * div[4][0] <= div[4][1]);no_break();
            // coverity[unterminated_case]
            case 4: divides |= (i * div[3][0] <= div[3][1]);no_break();
            // coverity[unterminated_case]
            case 3: divides |= (i * div[2][0] <= div[2][1]);no_break();
            // coverity[unterminated_case]
            case 2: divides |= (i * div[1][0] <= div[1][1]);no_break();
            // coverity[unterminated_case]
            case 1: divides |= (i * div[0][0] <= div[0][1]);no_break();
            case 0: break;
        }

        if (divides) {
            if (verify_gcd)
                ASSERT_ALWAYS(bin_gcd_int64_safe (i, j) != 1);
  #ifdef TRACE_K
            if (trace_on_spot_Nx(N, x)) {
                verbose_fmt_print(TRACE_CHANNEL, 0, "# Slot [{}] in bucket {} has non coprime (i,j)=({},{})\n",
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
                verbose_fmt_print(TRACE_CHANNEL, 0, "# Slot [{}] in bucket {} is survivor with coprime (i,j)\n",
                        x, N);
            }
  #endif
        }
      }
}


/* same, but with one side only.
 *
 *
 * XXX I think the overall structure would need a cleanup anyway, but at
 * least here we can refactor somewhat.
 */

MAYBE_UNUSED static void
search_survivors_in_line1_oneside(unsigned char * Sf,
        const unsigned char bound,
        unsigned int j, 
        int i0, int i1,
        int N MAYBE_UNUSED, j_divisibility_helper const & j_div,
        unsigned int td_max, std::vector<uint32_t> &survivors)
{
    unsigned int div[6][2], nr_div;

    nr_div = extract_j_div(div, j, j_div, 3, td_max);
    ASSERT_ALWAYS(nr_div <= 6);

    for (int x = 0; x < (i1 - i0); x++) {
        if (Sf[x] > bound) {
            Sf[x] = 255;
            continue;
        }

        /* The very small prime used in the bound pattern, and unsieving larger
           primes have not identified this as gcd(i,j) > 1. It remains to check
           the trial-divided primes. */
        const unsigned int i = abs (i0 + x);
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
                verbose_fmt_print(TRACE_CHANNEL, 0, "# Slot [{}] in bucket {} has non coprime (i,j)=({},{})\n",
                        x, N, i, j);
            }
  #endif
            Sf[x] = 255;
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


/* linefragment is used when logI > LOG_BUCKET_REGION ; it is a positive
 * integer multiple of 2^LOG_BUCKET_REGION ; so the sub-line is actually
 * for:
 * -I/2 + linefragment <= i < -I/2 + MIN(2^LOG_BUCKET_REGION, I)
 */
void
search_survivors_in_line(unsigned char * const SS[2], 
        const unsigned char bound[2],
        unsigned int j,
        int i0, int i1,
        int N, j_divisibility_helper const & j_div,
        unsigned int td_max, unsieve_data const & us,
        std::vector<uint32_t> &survivors, sublat_t sublat)
{
    ASSERT_ALWAYS(SS[0] || SS[1]);
    unsigned char * Sf = SS[0];

    if (SS[0] && SS[1]) {

        /* In line j = 0, only the coordinate (i, j) = (-1, 0) may survive */
        // FIXME: in sublat mode, this is broken!
        if (j == 0 && (!sublat.m)) {
            if (i0 <= 0 && i1 > 0) {
                unsigned char const s0 = SS[0][1-i0];
                unsigned char const s1 = SS[1][1-i0];
                memset(SS[0], 255, i1 - i0);
                if (s0 <= bound[0] && s1 <= bound[1]) {
                    SS[0][1 - i0] = s0;
                    SS[1][1 - i0] = s1;
                    survivors.push_back(1 - i0);
                }
            } else {
                memset(SS[0], 255, i1 - i0);
            }
            return;
        }

        // Naive version when we have sublattices, because unsieving is
        // harder. TODO: implement a fast version
        if (sublat.m) {
            for (int x = 0; x < (i1 - i0); x++) {
                if (!sieve_info_test_lognorm(bound[0], bound[1], SS[0][x], SS[1][x])) {
                    SS[0][x] = 255;
                    continue;
                }
                const unsigned int i = abs(int(sublat.m)*(i0 + x)+int(sublat.i0));
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

        unsieve_not_coprime_line(Sf, j, i0, i1, td_max + 1, us);

#if defined(HAVE_SSE2)
        search_survivors_in_line_sse2(SS, bound, j, i0, i1, N, j_div, td_max,
                survivors);
#else
        search_survivors_in_line1(SS, bound, j, i0, i1, N, j_div, td_max,
                survivors);
#endif
    } else {
        unsigned char b = bound[0];
        if (!Sf) {
            Sf = SS[1];
            b = bound[1];
        }
        /* ok, here only Sf is non-null. We want the values below b. */

        /* In line j = 0, only the coordinate (i, j) = (-1, 0) may survive */
        // FIXME: in sublat mode, this is broken!
        if (j == 0 && (!sublat.m)) {
            if (i0 <= 0 && i1 > 0) {
                unsigned char const s = Sf[1-i0];
                memset(Sf, 255, i1 - i0);
                if (s <= b) {
                    Sf[1 - i0] = s;
                    survivors.push_back(1 - i0);
                }
            } else {
                memset(Sf, 255, i1 - i0);
            }
            return;
        }

        // Naive version when we have sublattices, because unsieving is
        // harder. TODO: implement a fast version
        if (sublat.m) {
            for (int x = 0; x < (i1 - i0); x++) {
                if (Sf[x] > b) {
                    Sf[x] = 255;
                    continue;
                }
                const unsigned int i = abs(int(sublat.m)*(i0 + x))+int(sublat.i0);
                const unsigned int jj = sublat.m*j+sublat.j0;
                if ((((jj % 2) == 0) && ((i % 2) == 0)) ||
                        (bin_gcd_int64_safe (i, jj) != 1)) {
                    Sf[x] = 255;
                } else {
                    survivors.push_back(x);
                }
            }
            return;
        }

        unsieve_not_coprime_line(Sf, j, i0, i1, td_max + 1, us);

#if defined(HAVE_SSE2)
        search_survivors_in_line_sse2_oneside(Sf, b, j, i0, i1, N, j_div, td_max,
                survivors);
#else
        search_survivors_in_line1_oneside(Sf, b, j, i0, i1, N, j_div, td_max,
                survivors);
#endif
    }
}

/* assume SS[i] != NULL for all 0 <= i < nsides. */
template<std::size_t nsides>
void search_survivors_in_line(
        std::array<unsigned char * const, nsides> SS,
        const std::array<unsigned char, nsides> bound,
        unsigned int length,
        std::vector<uint16_t> &survivors)
{
    for (unsigned int x = 0; x < length; ++x) {
        if (!sieve_info_test_lognorm(bound, SS, x)) {
            SS[0][x] = 255;
        } else {
            survivors.push_back(x);
        }
    }
}

template void search_survivors_in_line(
        std::array<unsigned char * const, 1> SS,
        const std::array<unsigned char, 1> bound,
        unsigned int length,
        std::vector<uint16_t> &survivors);
template void search_survivors_in_line(
        std::array<unsigned char * const, 2> SS,
        const std::array<unsigned char, 2> bound,
        unsigned int length,
        std::vector<uint16_t> &survivors);
