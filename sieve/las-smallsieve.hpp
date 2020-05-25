#ifndef LAS_SMALLSIEVE_HPP_
#define LAS_SMALLSIEVE_HPP_

#include <cstdarg>
#include <cstdio>
#include <vector>
#include "fb-types.h"   // sublat_t
#include "fb.hpp"
#include "las-forwardtypes.hpp"
#include "macros.h"
class bucket_primes_t;
struct qlattice_basis;
struct small_sieve_data_t;
struct where_am_I;

/* Do not compute start positions for more than this number of bucket
 * regions in advance. This defines the frequency of a synchronization
 * point, so it should not be too small. Typically one set of start
 * positions for one bucket region costs about 25k.
 *
 * This is capped to nb_buckets
 */
#define SMALL_SIEVE_START_POSITIONS_MAX_ADVANCE 1024

extern void small_sieve_info(const char * what, int side, small_sieve_data_t const & r);
extern int small_sieve_dump(FILE *, const char *, va_list);
extern void small_sieve_clear(small_sieve_data_t & ssd);

extern void small_sieve_init(small_sieve_data_t & ssd,
                      std::vector<fb_entry_general> const & resieved,
                      std::vector<fb_entry_general> const & rest,
                      int logI,
                      int side,
                      fb_factorbase::key_type const & fbK,
                      qlattice_basis const & Q,
                      double scale);

extern void small_sieve_start(
        std::vector<spos_t> & ssdpos,
        small_sieve_data_t & ssd,
        unsigned int first_region_index,
        int logI, sublat_t const & sl);

extern void small_sieve_prepare_many_start_positions(
        small_sieve_data_t & ssd,
        unsigned int first_region_index,
        int nregions,
        int logI, sublat_t const & sl);

extern void small_sieve_activate_many_start_positions(small_sieve_data_t & ssd);

extern void sieve_small_bucket_region(unsigned char *S, unsigned int N,
                               small_sieve_data_t const & ssd,
                               std::vector<spos_t> const & ssdpos,
                               int logI, sublat_t const & sl,
                               where_am_I & w);

extern void resieve_small_bucket_region (bucket_primes_t *BP,
        unsigned char *S,
        unsigned int N,
        small_sieve_data_t & ssd,
        std::vector<spos_t> const & ssdpos,
        int logI, sublat_t const & sl,
        where_am_I & w MAYBE_UNUSED);

#endif	/* LAS_SMALLSIEVE_HPP_ */
