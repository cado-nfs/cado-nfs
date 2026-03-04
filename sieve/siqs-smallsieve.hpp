#ifndef CADO_SIQS_SMALLSIEVE_HPP
#define CADO_SIQS_SMALLSIEVE_HPP

#include <cstddef>
#include <memory>
#include <vector>

#include "fmt/base.h"

#include "fb-types.hpp"
#include "fb.hpp"
#include "las-forwardtypes.hpp"
#include "las-qlattice.hpp"
#include "smallsieve.hpp"

static_assert(std::is_same_v<siqs_pos_t, uint32_t>,
             "siqs-smallsieve.hpp assumes siqs_pos_t is uint32_t.");

class siqs_ssp_pdata_t {
    fbprime_t p;
    fbprime_t pmask; /* p-1 if p power of 2 else 0xff..ff */
    fbprime_t invq; /* 2^(word size)/q mod p */
    fbprime_t rq0; /* add(crt_data_modp)/2 mod p = add(crt_data_modq)/q mod p */
    /* crt_data_modp := [ 2*Rk/q mod p for Rk in Q.crt_data_modq ] */
    std::vector<fbprime_t> crt_data_modp;
    /* The only offset that can be used is -2^LOG_BUCKET_REGION. It is only
     * needed when logI > LOG_BUCKET_REGION. As it does not depend on the root
     * r, we store it here.
     */
    fbprime_t offset;

public:
    siqs_ssp_pdata_t(fb_entry_general const & e, siqs_special_q_data const & Q);

    fbprime_t get_p() const { return p; }
    fbprime_t get_pmask() const { return pmask; }
    fbprime_t get_offset() const { return offset; }
    std::vector<fbprime_t> const & CRT_data() const { return crt_data_modp; }
    /* Given a root rp modulo p, compute rp/q - r0 modulo p */
    fbroot_t transform_root(fbroot_t rp, redc_invp_t invp) const;
};


class siqs_ssp_simple_t {
public:
    siqs_ssp_simple_t(std::shared_ptr<siqs_ssp_pdata_t> const & pdata,
                      fbprime_t r,
                      unsigned char logp)
        : pdata(pdata), r(r), logp(logp)
    {}
    bool operator<(siqs_ssp_simple_t const & x) const {
        return get_p() < x.get_p();
    }
    /* Computes the position for the first bucket of the first region, i.e.,
     * the one with j=0 starting at i0=-I/2.
     * As r corresponds to the root for j = 0, it is just a matter of
     * substracting i0.
     * Returns r - i0 = r - (-I/2) = r + I/2
     */
    siqs_pos_t first_position_first_line(int logI) const {
        return ((r + (1u << (logI-1))) % pdata->get_p());
    }
    ///* Same but assumes p is a power of 2 */
    //siqs_pos_t first_position_first_line_power_of_two(int logI) const {
    //    return ((r + (1u << (logI-1))) & pdata->get_pmask());
    //}
    fbprime_t get_p() const { return pdata->get_p(); }
    fbprime_t get_r() const { return r; }
    fbprime_t get_offset() const { return pdata->get_offset(); }
    std::vector<fbprime_t> const & CRT_data() const { return pdata->CRT_data();}
    unsigned char get_logp() const { return logp; }

protected:
    std::shared_ptr<siqs_ssp_pdata_t> pdata;
    fbprime_t r;
    unsigned char logp;
};


class siqs_ssp_t : public siqs_ssp_simple_t {
public:
    siqs_ssp_t(std::shared_ptr<siqs_ssp_pdata_t> const & pdata,
               fbprime_t r,
               unsigned char logp);

    fbprime_t get_pmask() const { return pdata->get_pmask(); }
    bool is_pow2() const {return (flags & SSP_POW2) != 0;}
    bool is_nice() const {return !is_pow2();}
    bool is_pattern_sieved() const {return (flags & SSP_PATTERN_SIEVED) != 0;}

protected:
    uint8_t flags = 0u;
};

namespace fmt {
    template<>
    struct formatter<siqs_ssp_simple_t> : public formatter<string_view> {
        auto format(siqs_ssp_simple_t const & a, format_context & ctx) const
            -> format_context::iterator;
    };
    template<>
    struct formatter<siqs_ssp_t> : public formatter<string_view> {
        auto format(siqs_ssp_t const & a, format_context & ctx) const
            -> format_context::iterator;
    };
} /* namespace fmt */


class siqs_small_sieve_data : public small_sieve_data {
public:
    void small_sieve_init(
            std::vector<fb_entry_general> const & resieved,
            std::vector<fb_entry_general> const & rest,
            int logI,
            int side,
            fb_factorbase::key_type const & factorbaseK,
            siqs_special_q_data const & Q,
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

private:
    fb_factorbase::key_type fbK;
    std::vector<siqs_ssp_simple_t> ssps;
    std::vector<siqs_ssp_t> ssp;
    size_t resieve_end_offset;

    std::vector<std::vector<siqs_pos_t>> ssdpos_many;
    std::vector<std::vector<siqs_pos_t>> ssdpos_many_next;
};

#endif	/* CADO_SIQS_SMALLSIEVE_HPP */
