#ifndef LAS_SMALLSIEVE_HPP_
#define LAS_SMALLSIEVE_HPP_

#include <stdarg.h>
#include "fb-types.h"
#include "macros.h"

/* Structures for small sieve */

#define SSP_POW2        (1u<<0)
#define SSP_PROJ        (1u<<1)
#define SSP_DISCARD_SUBLAT     (1u<<14)
#define SSP_DISCARD_PROJ     (1u<<15)

class ssp_simple_t {
public:
    fbprime_t p;
    fbprime_t r;
    fbprime_t offset;
    unsigned char logp;

    void init(fbprime_t _p, fbprime_t _r, unsigned char _logp, unsigned int skip) {
        p = _p;
        r = _r;
        offset = (r * skip) % p;
        logp = _logp;
    }
    void print(FILE *) const;
};

class ssp_t {
    /* ordinary primes. Equation is (i-r*j) = 0 mod p */
    /* p may be a prime power, and even a power of two */
    /* Note that we need three fields for the projective primes anyway,
     * so we may have three here. Yet, admittedly the "offset" field is
     * rather useless. (see ssp_init_oa.)  We can recompute it on the fly
     * when needed.
     */
    fbprime_t p;	// q if projective
    fbprime_t r;        // in [ 0, p [. g if projective
    fbprime_t offset;   // in [ 0, p [. U if projective
    uint16_t flags;
public:
    unsigned char logp;

    fbprime_t get_p() const {ASSERT(!is_proj()); return p;}
    fbprime_t get_r() const {ASSERT(!is_proj()); return r;}
    fbprime_t get_offset() const {ASSERT(!is_proj()); return offset;}
    fbprime_t get_q() const {ASSERT(is_proj()); ASSERT(p > 0); return p;}
    fbprime_t get_g() const {ASSERT(is_proj()); ASSERT(r > 0); return r;}
    fbprime_t get_U() const {ASSERT(is_proj()); return offset;}
    void set_p(const fbprime_t _p) {p = _p;}
    void set_q(const fbprime_t q) {p = q;}
    void set_r(const fbprime_t _r) {r = _r;}
    void set_g(const fbprime_t g) {ASSERT(g > 0); r = g;}
    void set_offset(const fbprime_t _offset) {offset = _offset;}
    void set_U(const fbprime_t U) {offset = U;}

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
    
/* Initialization procedures for the ssp data */
    void init_affine(fbprime_t _p, fbprime_t _r, unsigned char _logp, unsigned int skip) {
        p = _p;
        r = _r;
        logp = _logp;
        offset = (r * skip) % p;
    }
    void init_proj(fbprime_t p, fbprime_t r, unsigned char _logp,
                   unsigned int skip MAYBE_UNUSED);
    void init(fbprime_t p, fbprime_t r, unsigned char logp, unsigned int skip,
              bool proj = false) {
        if (proj)
            init_proj(p, r, logp, skip);
        else
            init_affine(p, r, logp, skip);
    }

    void print(FILE *) const;
};

typedef struct {
    ssp_simple_t *ssps;
    ssp_t *ssp;
    int nb_ssps, nb_ssp;
} small_sieve_data_t;

/* Include this only now, as there's a cross dependency between the two
 * (our prototypes need sieve_info_t, which needs our datatypes...)
 */
#include "las-types.hpp"
#include "bucket.hpp"

extern void small_sieve_info(const char * what, int side, small_sieve_data_t * r);
extern int small_sieve_dump(FILE *, const char *, va_list);
extern void small_sieve_clear(small_sieve_data_t * ssd);
extern void small_sieve_extract_interval(small_sieve_data_t * r, small_sieve_data_t * s, int bounds[2]);
extern void small_sieve_init(small_sieve_data_t *ssd, unsigned int interleaving, const std::vector<fb_general_entry> *fb,
                      sieve_info const & si, int side);
extern int64_t * small_sieve_copy_start(int64_t * base, int bounds[2]);
extern int64_t * small_sieve_start(small_sieve_data_t *ssd, unsigned int N, sieve_info const & si);
/*
extern void small_sieve_skip_stride(small_sieve_data_t *ssd, int64_t * ssdpos, unsigned int N, unsigned int interleaving, sieve_info const & si);
*/
extern void sieve_small_bucket_region(unsigned char *S, int N,
			       small_sieve_data_t * ssd, int64_t * ssdpos, sieve_info& si,
                               int side,
                               int interleaving,
			       where_am_I& w MAYBE_UNUSED);

extern void
resieve_small_bucket_region (bucket_primes_t *BP, int N, unsigned char *S,
        small_sieve_data_t *ssd, int64_t * ssdpos,
        sieve_info const & si,
        int interleaving,
        where_am_I& w MAYBE_UNUSED);


#endif	/* LAS_SMALLSIEVE_HPP_ */
