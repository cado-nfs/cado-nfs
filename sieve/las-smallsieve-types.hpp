#ifndef LAS_SMALLSIEVE_TYPES_HPP_
#define LAS_SMALLSIEVE_TYPES_HPP_

#include <vector>
#include "las-forwardtypes.hpp"
#include "fb-types.h"
#include "fb.hpp"
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
    /* Equation for ordinary primes is (i-r*j) = 0 mod p.
     * We let L_p be the p-lattice.
     */
    fbprime_t p;
    fbprime_t r;
public:
    unsigned char logp;

    ssp_simple_t() : p(0), r(0), logp(0) {}
    ssp_simple_t(fbprime_t _p, fbprime_t _r, unsigned char _logp)
    : p(_p), r(_r), logp(_logp)
    {}
    fbprime_t get_p() const {return p;}
    fbprime_t get_r() const {return r;}
    void set_p(const fbprime_t _p) {p = _p;}
    void set_r(const fbprime_t _r) {r = _r;}
    bool is_nice() const {return true;}
    void print(FILE *) const;
    bool operator<(ssp_simple_t const& x) const {
        return p < x.p;
    }
};

class ssp_t : public ssp_simple_t {
    fbprime_t offset;   /* we used to have that in ssp_simple_t. Now it's
                           no longer here, so that ssp_t and ssp_simple_t
                           no longer have the same size. This used to be
                           a requirement that they do, but I think it's
                           not the case anymore.  */
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
    ssp_t(fbprime_t _p, fbprime_t _r, unsigned char _logp)
    : ssp_simple_t(_p, _r, _logp)
    {}
    /* Constructor for affine or projective case */
    ssp_t(fbprime_t _p, fbprime_t _r, unsigned char _logp, bool proj);

    /* We could use the parent class' methods if we get rid of the ASSERT()s */
    fbprime_t get_p() const {ASSERT(!is_proj()); return p;}
    fbprime_t get_r() const {ASSERT(!is_proj()); return r;}

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

// static_assert(sizeof(ssp_simple_t) == sizeof(ssp_t), "struct padding has been tampered with");

/* The prototypes in las-smallsieve.hpp should all be member functions
 * for this struct.  However we used to have an interface nesting problem
 * with sieve_info, perhaps it's gone now [FIXME, check]
 */
struct small_sieve_data_t {
    fb_factorbase::key_type fbK;
    std::vector<ssp_simple_t> ssps;
    std::vector<ssp_t> ssp;
    /* This counts the resieved primes in ssps */
    size_t resieve_end_offset;

    /* We have some vectors of small sieve positions prepared in
     * advanced (up to nb_buckets[1] of them). The ssdpos_many_next
     * area is for staging the next set of start positions, possible
     * while threads are using the positions in ssdpos_many.
     */
    std::vector<std::vector<spos_t>> ssdpos_many;
    std::vector<std::vector<spos_t>> ssdpos_many_next;
    size_t ssdpos_many_offset;

    /* This vectors are computed at the same time as ssd is
     * intialized, in small_sieve_start.  We use it to compute the
     * ssdpos fields for primes in ssps (easy ones).
     *
     * for logI <= logB:
     * c=offsets[i] is equal to ssd[i].offset when logI>logB.
     * In full generality, if logB = v + min(logI, logB), then c is
     * such that (c, 2^v) in L_p.
     *
     * for logI > logB:
     * c=offsets[i] is such that for the i-th p-lattice L_p, the
     * vector (2^min(logB, logI) + c, 0) is in L_p. We do not need it at
     * all when logI <= logB, so we may as well define it as (B+c,p)
     * in L_p.
     *
     */
    
    /* Note that at most one of these two vectors will be used anyway,
     * depnding on how logI and logB compare.
     */
    std::vector<fbprime_t> offsets;
};

#endif
