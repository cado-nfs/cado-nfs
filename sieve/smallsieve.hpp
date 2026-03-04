#ifndef CADO_SMALLSIEVE_HPP
#define CADO_SMALLSIEVE_HPP

#include "fb.hpp"
#include "las-forwardtypes.hpp"
#include "las-qlattice.hpp"

#define SSP_POW2            (1u<<0)
#define SSP_PROJ            (1u<<1)
#define SSP_PATTERN_SIEVED  (1u<<2)

/* Do not compute start positions for more than this number of bucket
 * regions in advance. This defines the frequency of a synchronization
 * point, so it should not be too small. Typically one set of start
 * positions for one bucket region costs about 25k.
 *
 * This is capped to nb_buckets
 */
#define SMALL_SIEVE_START_POSITIONS_MAX_ADVANCE 1024


class small_sieve_data {
public:
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
                                 "qlattice_basis as special_q_data_class");
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
            sublat_t const & sl) = 0;

    virtual void small_sieve_activate_many_start_positions() = 0;

    virtual void sieve_small_bucket_region(
            unsigned char *S,
            unsigned int N,
            int bucket_relative_index,
            int logI,
            sublat_t const & sl,
            where_am_I & w) const = 0;

    virtual void resieve_small_bucket_region(
            bucket_primes_t *BP,
            unsigned char *S,
            unsigned int N,
            int bucket_relative_index,
            int logI,
            sublat_t const & sl,
            where_am_I & w MAYBE_UNUSED) = 0;

    virtual ~small_sieve_data() = default;
};

#endif	/* CADO_SMALLSIEVE_HPP */
