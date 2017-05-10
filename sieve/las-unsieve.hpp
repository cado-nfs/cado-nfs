#ifndef LAS_UNSIEVE_HPP_
#define LAS_UNSIEVE_HPP_

#include <stdint.h>
#include <vector>
#include "ularith.h"
#include "fb.hpp" // for sublat_t

#ifdef HAVE_SSE2
#include <xmmintrin.h>
#endif

struct unsieve_data {
    struct entry {
        unsigned int lpf, cof, start;
        entry(unsigned int lpf = 0,
                unsigned int cof = 0,
                unsigned int start = 0)
                : lpf(lpf), cof(cof), start(start) {}
    };
#ifdef HAVE_SSE2
typedef __m128i pattern_t;
// oh my...
// http://stackoverflow.com/questions/2804902/whats-the-difference-between-logical-sse-intrinsics
#define UNSIEVE_OR(a,b) do{(a) = _mm_or_si128((a),(b));}while(0)
#else
typedef unsigned long pattern_t;
#define UNSIEVE_OR(a,b) do{(a) |= (b);}while(0)
#endif
    /* entry[i].lpf is largest prime factor of i, for i < Jmax,
       cof is the cofactor i/lpf^k s.t. lfp^k || i,
       start is (I/2) % lpf */
    uint32_t Jmax;
    entry *entries;
    pattern_t pattern3[3];
    pattern_t pattern5[5];
    pattern_t pattern7[7];
    unsieve_data();
    unsieve_data(int logI, int logA);
    unsieve_data(unsieve_data const &);
    unsieve_data& operator=(unsieve_data const &);
    ~unsieve_data();
};

/* {{{ j_divisibility_helper
 * A structure for factoring the line-coordinate j into distinct odd primes
 * and performing fast divisibility tests of i by those primes.
 */
struct j_divisibility_helper {
    struct entry {
        unsigned int p,     /* The odd prime that divides this entry */
                     cof,   /* The cofactor after dividing out p as often as possible */
                     inv,   /* p^(-1) (mod 2^32) */
                     bound; /* (2^32-1) / p */
    };
    uint32_t J;
    entry *entries;
    entry& operator[](int i) { return entries[i]; }
    entry const & operator[](int i) const { return entries[i]; }
    j_divisibility_helper(uint32_t J = 0);
    j_divisibility_helper(j_divisibility_helper const &);
    j_divisibility_helper& operator=(j_divisibility_helper const &);
    ~j_divisibility_helper();
};


/* Helper function shared by las-unsieve.cpp and las-unsieve-sse2.cpp.
   Polluting global namespace is a little unfortunate here - #ifdef? */
static inline unsigned int 
extract_j_div(unsigned int (*div)[2], const unsigned int j, j_divisibility_helper const & j_div, 
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


j_divisibility_helper * init_j_div(uint32_t);
void clear_j_div(j_divisibility_helper *);
void search_survivors_in_line(unsigned char * const restrict[2],
        const unsigned char[2], unsigned int, unsigned int, int,
        j_divisibility_helper const &, unsigned int, unsieve_data const &,
        std::vector<uint32_t> &, sublat_t);
#ifdef HAVE_SSE2 
void search_survivors_in_line_sse2(unsigned char * const restrict[2],
        const unsigned char[2], unsigned int, unsigned int, int,
        j_divisibility_helper const &,
        unsigned int, std::vector<uint32_t> &);
#endif

#endif	/* LAS_UNSIEVE_HPP_ */
