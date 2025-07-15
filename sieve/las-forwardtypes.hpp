#ifndef CADO_LAS_FORWARDTYPES_HPP
#define CADO_LAS_FORWARDTYPES_HPP

#include <cstdint>

/* These must be forward-declared, because various header files use them */
/* TODO: still true ?
 *
 * Quite often, the forward declarations are just put as they are. It's
 * not totally clear to me which solution is better, between having this
 * "forwardtypes" header, and declaring them one by one when needed.
 */

class bucket_primes_t;
struct las_info;
struct where_am_I;;
// struct siever_config;;
// struct las_dlog_base;
struct special_q_task_collection;
class nfs_work;
class nfs_work_cofac;
class nfs_aux;

/* we use this type to store positions within the sieve arrays. This is
 * typically no more than 2*factor base bound, even though for the
 * largest deal of the processing needs, it is actually less than
 * LOG_BUCKET_REGION
 */
typedef int32_t spos_t;
typedef int64_t long_spos_t;

#endif
