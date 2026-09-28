#ifndef CADO_SMALLSIEVE_HPP
#define CADO_SMALLSIEVE_HPP

#include <vector>
#include <stdexcept>

#include "fb.hpp"
#include "las-forwardtypes.hpp"
#include "las-qlattice.hpp"
#include "threadpool.hpp"
#include "fb-types.hpp"
#include "macros.h"
#include "las-sublat.hpp"

#define SSP_POW2            (1u<<0)
#define SSP_PROJ            (1u<<1)
#define SSP_PATTERN_SIEVED  (1u<<2)
#define SSP_SUBLAT_ROWWISE  (1u<<3)

class small_sieve_data {
public:
    /* Sum of the logs of the prime powers that divide the sublattice
     * modulus. Such a prime power divides either every norm of the
     * current sublattice class, or none of them, so sieving it would
     * merely add a constant to the whole sieve region. We leave it out
     * of the sieve and add the constant to the survivor bound instead.
     * Zero when there are no sublattices. */
    unsigned char constant_logp = 0;

    virtual void small_sieve_init(
            std::vector<fb_entry_general> const &,
            std::vector<fb_entry_general> const &,
            int,
            int,
            fb_factorbase::key_type const &,
            qlattice_basis const &,
            double)
    {
        throw std::runtime_error("small_sieve_init has not been implemented for"
                                 "siqs_special_q_data as special_q_data_class");
    }

    virtual void small_sieve_init(
            std::vector<fb_entry_general> const &,
            std::vector<fb_entry_general> const &,
            int,
            int,
            fb_factorbase::key_type const &,
            siqs_special_q_data const &,
            double)
    {
        throw std::runtime_error("small_sieve_init has not been implemented for"
                                 "siqs_special_q_data as special_q_data_class");
    }

    virtual void small_sieve_clear() = 0;

    virtual void small_sieve_info(const char * what, int side) const = 0;

    virtual void small_sieve_prepare_many_start_positions(
            unsigned int first_region_index,
            int nregions,
            int logI,
            sublat_runtime_t const & sl) = 0;

    virtual void small_sieve_prepare_many_start_positions(
            thread_pool &,
            task_group *,
            unsigned int first_region_index,
            int nregions,
            int logI,
            sublat_runtime_t const & sl)
    {
        // by default, fall back to the single-threaded code
        small_sieve_prepare_many_start_positions(first_region_index, nregions, logI, sl);
    }
    virtual void small_sieve_activate_many_start_positions() = 0;

    virtual void sieve_small_bucket_region(
            unsigned char *S,
            unsigned int N,
            int bucket_relative_index,
            int logI,
            sublat_runtime_t const & sl,
            where_am_I & w) const = 0;

    virtual void resieve_small_bucket_region(
            bucket_primes_t *BP,
            unsigned char *S,
            unsigned int N,
            int bucket_relative_index,
            int logI,
            sublat_runtime_t const & sl,
            where_am_I & w MAYBE_UNUSED) = 0;

    virtual ~small_sieve_data() = default;
};

#endif	/* CADO_SMALLSIEVE_HPP */
