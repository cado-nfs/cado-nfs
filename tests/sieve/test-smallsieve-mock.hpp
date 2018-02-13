#ifndef TEST_SMALLSIEVE_MOCK_HPP_
#define TEST_SMALLSIEVE_MOCK_HPP_

#include <stdint.h>
#include "macros.h"

#include "las-smallsieve.hpp"

#if 0
/* some mock datatypes */
#define WHERE_AM_I_UPDATE(a,b,c)  /**/
typedef void * where_am_I;

struct sublat_t {
    uint32_t m = 0;
    uint32_t i0 = 0;
    uint32_t j0 = 0;
};

typedef int64_t spos_t;

static inline void sieve_increase(unsigned char *S, const unsigned char logp, where_am_I& w MAYBE_UNUSED)
{
    *S += logp;
}
#endif

#ifdef LOG_BUCKET_REGION_IS_A_CONSTANT
#define LOG_BUCKET_REGION 16
#else
int LOG_BUCKET_REGION = 16;
#endif

/* this does the small sieve loop, only for nice primes.
 * In theory, we could do this stuff for all primes, but we'll be happy
 * enough if we can bench code for the main branch that goes fast enough
 */

typedef uint32_t fbprime_t;

#if 0
#define SSP_POW2        (1u<<0)
#define SSP_PROJ        (1u<<1)
#define SSP_DISCARD_SUBLAT     (1u<<2)
#define SSP_DISCARD_PROJ     (1u<<3)

struct ssp_t {
    /* ordinary primes. Equation is (i-r*j) = 0 mod p */
    /* p may be a prime power, and even a power of two */
    /* Note that we need three fields for the projective primes anyway,
     * so we may have three here. Yet, admittedly the "offset" field is
     * rather useless. (see ssp_init_oa.)  We can recompute it on the fly
     * when needed.
     */
    fbprime_t p;	// q if projective
    fbprime_t r;        // in [ 0, p [. g if projective
    // unused:
    fbprime_t offset;   // in [ 0, p [. U if projective
    uint16_t flags=0;
    public:
    unsigned char logp;
    ssp_t(fbprime_t p, fbprime_t r, fbprime_t offset = 0)
        : p(p), r(r), offset(offset) {
            // we don't care if we put garbage.
            logp = p/3;
        }
    ssp_t() = default;
    /* this code does not deal with proj primes anyway, but the accessors
     * below are real */
    inline bool is_proj() const { return false; }

    fbprime_t get_p() const {ASSERT(!is_proj()); return p;}
    fbprime_t get_r() const {ASSERT(!is_proj()); return r;}
    fbprime_t get_offset() const {ASSERT(!is_proj()); return offset;}
    fbprime_t get_q() const {ASSERT(is_proj()); ASSERT(p > 0); return p;}
    fbprime_t get_g() const {ASSERT(is_proj()); ASSERT(r > 0); return r;}
    fbprime_t get_U() const {ASSERT(is_proj()); return offset;}

    void set_q(const fbprime_t q) {ASSERT(is_proj()); p = q;}
    void set_g(const fbprime_t g) {ASSERT(is_proj()); ASSERT(g > 0); r = g;}
    void set_U(const fbprime_t U) {ASSERT(is_proj()); offset = U;}

    bool is_pow2() const {return (flags & SSP_POW2) != 0;}
    bool is_proj() const {return (flags & SSP_PROJ) != 0;}
    bool is_discarded_sublat() const {return (flags & SSP_DISCARD_SUBLAT) != 0;}
    bool is_discarded_proj() const {return (flags & SSP_DISCARD_PROJ) != 0;}
    bool is_discarded() const {return is_discarded_proj() || is_discarded_sublat();}
    bool is_nice() const {return !is_pow2() && !is_proj() && !is_discarded();}

    void set_pow2() {flags |= SSP_POW2;}
    void set_proj() {flags |= SSP_PROJ;}
    void set_discarded_sublat() {flags |= SSP_DISCARD_SUBLAT;}
    void set_discarded() {flags |= SSP_DISCARD_PROJ;}
    
    void print(FILE *) const;
};
#endif


#endif	/* TEST_SMALLSIEVE_MOCK_HPP_ */
