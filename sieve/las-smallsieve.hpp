#ifndef CADO_LAS_SMALLSIEVE_HPP
#define CADO_LAS_SMALLSIEVE_HPP

#include <cstddef>
#include <vector>

#include "fmt/base.h"

#include "fb-types.hpp"
#include "fb.hpp"
#include "las-forwardtypes.hpp"
#include "las-qlattice.hpp"
#include "macros.h"
#include "smallsieve.hpp"


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
    bool operator<(ssp_simple_t const& x) const {
        return p < x.p;
    }
    friend struct fmt::formatter<ssp_simple_t>;
};
namespace fmt {
    template<>
    struct formatter<ssp_simple_t> : public formatter<string_view> {
        auto format(ssp_simple_t const & a, format_context & ctx) const
            -> format_context::iterator;
    };
} /* namespace fmt */

class ssp_t : public ssp_simple_t {
    fbprime_t offset = 0;  /* we used to have that in ssp_simple_t. Now it's
                              no longer here, so that ssp_t and
                              ssp_simple_t no longer have the same size.
                              This used to be a requirement that they do,
                              but I think it's not the case anymore.  */
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
    ssp_t() = default;
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
    bool is_nice() const {return !is_pow2() && !is_proj();}
    bool is_pow() const { return rootp != 0; }
    bool is_pattern_sieved() const {return (flags & SSP_PATTERN_SIEVED) != 0;}

    void set_pow(unsigned int p) { rootp = p; }
    void set_pow2() {flags |= SSP_POW2; rootp=2;}
    void set_proj() {flags |= SSP_PROJ;}
    void set_pattern_sieved() {flags |= SSP_PATTERN_SIEVED;}

private:
    void init_proj(fbprime_t p, fbprime_t r, unsigned char _logp,
                   unsigned int skip MAYBE_UNUSED);
    friend struct fmt::formatter<ssp_t>;
};
namespace fmt {
    template<>
    struct formatter<ssp_t> : public formatter<string_view> {
        auto format(ssp_t const & a, format_context & ctx) const
            -> format_context::iterator;
    };
} /* namespace fmt */

class las_small_sieve_data : public small_sieve_data {
public:
    void small_sieve_init(
            std::vector<fb_entry_general> const & resieved,
            std::vector<fb_entry_general> const & rest,
            int logI,
            int side,
            fb_factorbase::key_type const & factorbaseK,
            qlattice_basis const & Q,
            double scale) final;

    void small_sieve_clear() final;

    void small_sieve_info(const char * what, int side) const final;

    void small_sieve_prepare_many_start_positions(
            unsigned int first_region_index,
            int nregions,
            int logI,
            sublat_t const & sl) final;

    void small_sieve_activate_many_start_positions() final;

    void sieve_small_bucket_region(
            unsigned char *S,
            unsigned int N,
            int bucket_relative_index,
            int logI,
            sublat_t const & sl,
            where_am_I & w) const final;

    void resieve_small_bucket_region(
            bucket_primes_t *BP,
            unsigned char *S,
            unsigned int N,
            int bucket_relative_index,
            int logI,
            sublat_t const & sl,
            where_am_I & w MAYBE_UNUSED) final;

protected:
    void small_sieve_start(
            std::vector<spos_t> & ssdpos,
            unsigned int first_region_index,
            int logI,
            sublat_t const & sl);

    void small_sieve_print_contents(const char * prefix) const;

private:
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

#endif	/* CADO_LAS_SMALLSIEVE_HPP */
