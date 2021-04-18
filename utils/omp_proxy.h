#ifndef OMP_PROXY_H_
#define OMP_PROXY_H_

#include "macros.h"

// IWYU pragma: begin_exports
#ifdef HAVE_OPENMP
#include <omp.h>
#else
#ifdef  __GNUC__
#pragma GCC diagnostic ignored "-Wunknown-pragmas"
#endif
/* minimal stub */
#ifdef __cplusplus
extern "C" {
#endif
static inline int omp_get_max_threads()
{
  return 1;
}

static inline int omp_get_num_threads()
{
  return 1;
}

static inline void omp_set_num_threads(int n MAYBE_UNUSED)
{
}

static inline void omp_set_nested(int n MAYBE_UNUSED)
{
}

static inline int omp_get_thread_num()
{
  return 0;
}
#ifdef __cplusplus
}
#endif
#endif /* HAVE_OPENMP */
// IWYU pragma: end_exports


#endif	/* OMP_PROXY_H_ */
