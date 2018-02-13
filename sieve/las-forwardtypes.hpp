#ifndef LAS_FORWARDTYPES_HPP_
#define LAS_FORWARDTYPES_HPP_

/* These must be forward-declared, because various header files use them */
/* TODO: still true ? */

class bucket_primes_t;
struct sieve_info;
struct las_info;
struct where_am_I;;
// struct siever_config;;
#ifdef  DLP_DESCENT
// struct las_dlog_base;
struct descent_tree;
#endif

/* we use this type to store positions within the sieve arrays. This is
 * typically no more than 2*factor base bound, even though for the
 * largest deal of the processing needs, it is actually less than
 * LOG_BUCKET_REGION
 */
typedef int32_t spos_t;
typedef int64_t long_spos_t;

#endif
