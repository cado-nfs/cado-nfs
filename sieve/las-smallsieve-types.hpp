#ifndef LAS_SMALLSIEVE_TYPES_HPP_
#define LAS_SMALLSIEVE_TYPES_HPP_

#include <vector>
#include "fb-types.h"
#include "macros.h"

/* Structures for small sieve */

#define SSP_POW2        (1u<<0)
#define SSP_PROJ        (1u<<1)
#define SSP_DISCARD_SUBLAT     (1u<<2)
#define SSP_DISCARD_PROJ     (1u<<3)
#define SSP_PATTERN_SIEVED   (1u<<4)

/* spos_t is in las-forwardtypes.hpp */

/* Simple primes/roots. These are implicitly "nice", i.e., odd primes/powers
   with affine root.
   We also use them only for primes/roots which are sieved in the normal
   one-hit-at-a-time line sieving code. Pattern-sieved stuff is kept in
   ssp_t. */
class ssp_simple_t {
protected:
    /* Equation for ordinary primes is (i-r*j) = 0 mod p */
    fbprime_t p;
    fbprime_t r;
    fbprime_t offset;
public:
    unsigned char logp;

    ssp_simple_t() : p(0), r(0), offset(0), logp(0) {}
    ssp_simple_t(fbprime_t _p, fbprime_t _r, unsigned char _logp, unsigned int skip)
    : p(_p), r(_r), offset((r*skip)%p), logp(_logp)
    {}
    fbprime_t get_p() const {return p;}
    fbprime_t get_r() const {return r;}
    fbprime_t get_offset() const {return offset;}
    void set_p(const fbprime_t _p) {p = _p;}
    void set_r(const fbprime_t _r) {r = _r;}
    void set_offset(const fbprime_t _offset) {offset = _offset;}
    bool is_nice() const {return true;}
    void print(FILE *) const;
    bool operator<(ssp_simple_t const& x) const {
        return p < x.p;
    }
};

class ssp_t : public ssp_simple_t {
    /* use the remaining empty space in the struct so that we still have
     * the same size */
    uint8_t flags = 0;
    uint16_t rootp = 0;
    /* forbid comparison of ssp_t -- makes little sense I believe, as
     * it's a mixed bag. */
    bool operator<(ssp_simple_t const&) const { return false; }
public:

    /* Initialization procedures for the ssp data */
    /* Constructor for affine case */
    ssp_t() : ssp_simple_t(), flags(0) {}
    ssp_t(fbprime_t _p, fbprime_t _r, unsigned char _logp, unsigned int skip)
    : ssp_simple_t(_p, _r, _logp, skip)
    {}
    /* Constructor for affine or projective case */
    ssp_t(fbprime_t _p, fbprime_t _r, unsigned char _logp, unsigned int skip, bool proj);

    /* We could use the parent class' methods if we get rid of the ASSERT()s */
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
    bool is_pow() const { return rootp != 0; }
    bool is_pattern_sieved() const {return (flags & SSP_PATTERN_SIEVED) != 0;}

    void set_pow(unsigned int p) { rootp = p; }
    void set_pow2() {flags |= SSP_POW2; rootp=2;}
    void set_proj() {flags |= SSP_PROJ;}
    void set_discarded_sublat() {flags |= SSP_DISCARD_SUBLAT;}
    void set_discarded() {flags |= SSP_DISCARD_PROJ;}
    void set_pattern_sieved() {flags |= SSP_PATTERN_SIEVED;}
    
    void print(FILE *) const;
private:
    void init_proj(fbprime_t p, fbprime_t r, unsigned char _logp,
                   unsigned int skip MAYBE_UNUSED);
};

static_assert(sizeof(ssp_simple_t) == sizeof(ssp_t), "struct padding has been tampered with");

typedef struct {
    std::vector<ssp_simple_t> ssps;
    std::vector<ssp_t> ssp;
    /* These offsets, relative to the start of ssps, tell which of the ssps
       entries should be used for re-sieving */
    size_t resieve_start_offset, resieve_end_offset;
} small_sieve_data_t;

#endif
